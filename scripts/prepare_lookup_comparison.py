"""Prepare a controlled thermodynamic representation comparison, without altering baseline.

Lookup is a balanced interval search over cubic Hermite table coefficients.
The generated MOOSE parsed material differentiates THAT SAME piecewise energy.
No global polynomial refit and no independent chemical-potential interpolation.
"""
from pathlib import Path
import json, re, hashlib
import numpy as np
import pandas as pd
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import root

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT/'examples/ni_al_lookup'
OUT = ROOT/'results/ni_al_lookup'
CASE.mkdir(exist_ok=True); OUT.mkdir(exist_ok=True)
base = (ROOT/'examples/ni_al_coarsening/coarsening.i').read_text()
params = json.loads((ROOT/'results/ni_al_coarsening/model_parameters.json').read_text())
report = {'method':'piecewise cubic Hermite interpolation of reduced G and dG/dc, differentiated consistently', 'tables':{}}
expressions = {}
splines = []

def fmt(v): return format(float(v), '.17g')
def tree(x, coefficients, var, lo, hi):
    if hi-lo == 1:
        d = f'({var}-({fmt(x[lo])}))'
        a,b,c,e = coefficients[:,lo]
        return f'((({fmt(a)}*{d}+({fmt(b)}))*{d}+({fmt(c)}))*{d}+({fmt(e)}))'
    mid = (lo+hi)//2
    return f'if({var}<{fmt(x[mid])},{tree(x,coefficients,var,lo,mid)},{tree(x,coefficients,var,mid,hi)})'

for phase,var in [('gamma','cm'),('gamma_prime','cd')]:
    path = ROOT/f'thermodynamics/ni_al_1000K/{phase}_1000K.csv'
    df = pd.read_csv(path)
    x = df.x_Al.to_numpy()
    spline = CubicHermiteSpline(x,df.Gm_reduced_J_mol/1000,df.dGm_reduced_dx_J_mol/1000)
    splines.append(spline)
    expressions[phase] = tree(x,spline.c,var,0,len(x)-1)
    pd.DataFrame({'left_x_Al':x[:-1],'right_x_Al':x[1:], 'a3':spline.c[0], 'a2':spline.c[1], 'a1':spline.c[2], 'a0':spline.c[3]}).to_csv(CASE/f'{phase}_hermite_coefficients.csv',index=False)
    report['tables'][phase] = {'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'rows':len(df),'range':[float(x[0]),float(x[-1])], 'maximum_G_node_error_J_mol':float(np.max(np.abs(spline(x)*1000-df.Gm_reduced_J_mol))), 'maximum_mu_node_error_J_mol':float(np.max(np.abs(spline(x,1)*1000-df.dGm_reduced_dx_J_mol)))}

g,p = splines
def tangent(z):
    a,b=z; m=g(a,1)
    return [g(a,1)-p(b,1),g(a)-m*a-p(b)+m*b]
sol=root(tangent,[.124644,.237162],tol=1e-11)
assert sol.success and np.max(np.abs(tangent(sol.x)))<1e-10
report['common_tangent_x_Al'] = sol.x.tolist()
report['common_tangent_residual'] = tangent(sol.x)
report['curvatures_at_tangent_J_mol']=[float(g(sol.x[0],2)*1000),float(p(sol.x[1],2)*1000)]
report['quadratic_equilibrium_x_Al']=[.124644,.237162]
report['quadratic_curvatures_J_mol']=[347973.,424347.4]
report['controls'] = 'Same nodal eta,c IC, geometry, mesh, W,kappa,M,L, scales and final time. Replay baseline accepted time sequence.'
(OUT/'thermodynamic_validation.json').write_text(json.dumps(report,indent=2)+'\n')

kernels = '''[Kernels]
  [chemical_potential]
    type = NestedKKSSplitCHCRes
    variable = c
    global_cs = c
    w = mu
    all_etas = eta
    ca_names = 'cm cd'
    fa_name = fm
    coupled_variables = 'eta mu'
  []
  [conservation]
    type = CoupledTimeDerivative
    variable = mu
    v = c
  []
  [diffusion]
    type = SplitCHWRes
    variable = mu
    mob_name = M
  []
  [interface_time]
    type = TimeDerivative
    variable = eta
  []
  [bulk_f]
    type = NestedKKSACBulkF
    variable = eta
    global_cs = c
    ci_names = 'cm cd'
    fa_name = fm
    fb_name = fd
    mob_name = L
    w = 0.4614171612406061
    coupled_variables = c
  []
  [bulk_c]
    type = NestedKKSACBulkC
    variable = eta
    global_cs = c
    ci_names = 'cm cd'
    fa_name = fm
    mob_name = L
    coupled_variables = c
  []
  [interface_gradient]
    type = ACInterface
    variable = eta
    kappa_name = kappa_eta
    mob_name = L
  []
[]
'''
def materials(expr_g,expr_p):
    return f'''[Materials]
  [fm]
    type = DerivativeParsedMaterial
    property_name = fm
    expression = '{expr_g}'
    material_property_names = cm
    additional_derivative_symbols = cm
    derivative_order = 2
    compute = false
    enable_jit = true
  []
  [fd]
    type = DerivativeParsedMaterial
    property_name = fd
    expression = '{expr_p}'
    material_property_names = cd
    additional_derivative_symbols = cd
    derivative_order = 2
    compute = false
    enable_jit = true
  []
  [phase_concentrations]
    type = KKSPhaseConcentrationMaterial
    global_cs = c
    ci_names = 'cm cd'
    ci_IC = '0.124644 0.237162'
    fa_name = fm
    fb_name = fd
    h_name = h
    absolute_tolerance = 1e-12
    relative_tolerance = 1e-12
    min_iterations = 1
    max_iterations = 50
    damped_Newton = true
    conditions = table_domain
    damping_factor = 0.5
  []
  [table_domain]
    type = DerivativeParsedMaterial
    property_name = table_domain
    material_property_names = 'cm cd'
    expression = '(cm>=0.05)&(cm<=0.30)&(cd>=0.18)&(cd<=0.30)'
    compute = false
  []
  [phase_derivatives]
    type = KKSPhaseConcentrationDerivatives
    global_cs = c
    eta = eta
    ci_names = 'cm cd'
    fa_name = fm
    fb_name = fd
    h_name = h
  []
  [switch]
    type = SwitchingFunctionMaterial
    h_order = HIGH
    eta = eta
  []
  [barrier]
    type = BarrierFunctionMaterial
    g_order = SIMPLE
    eta = eta
  []
  [constants]
    type = GenericConstantMaterial
    prop_names = 'M L kappa_eta'
    prop_values = '0.002873786184560299 10 0.7646009503665432'
  []
[]
'''

times = pd.read_csv(ROOT/'results/ni_al_coarsening/diagnostics.csv').time.to_numpy()
pd.DataFrame({'time_star':times}).to_csv(CASE/'baseline_time_sequence.csv',index=False)
for name,eg,ep in [('lookup',expressions['gamma'],expressions['gamma_prime']),('quadratic_control','0.5*347.9729999999997*(cm-0.124644)^2','0.5*424.3473999999995*(cd-0.237162)^2')]:
    text = re.sub(r'\[Kernels\].*?(?=\[Materials\])',kernels,base,flags=re.S)
    text = re.sub(r'\[Materials\].*?(?=\[AuxVariables\])',lambda m:materials(eg,ep),text,flags=re.S)
    # Phase concentrations are material properties at quadrature points.
    for aux,var,prop in [('matrix_composition','c_gamma','cm'),('precipitate_composition','c_gamma_prime','cd')]:
        text = re.sub(r'  \['+aux+r'\].*?  \[\]',f'  [{aux}]\n    type = MaterialRealAux\n    variable = {var}\n    property = {prop}\n    execute_on = \'INITIAL TIMESTEP_END\'\n  []',text,flags=re.S)
    for var in ['c_gamma','c_gamma_prime']:
        text=text.replace(f'  [{var}]\n  []',f'  [{var}]\n    family = MONOMIAL\n    order = CONSTANT\n  []')
    text=text.replace('type = TotalFreeEnergy\n    variable = energy_density\n    f_name = F','type = KKSGlobalFreeEnergy\n    variable = energy_density\n    fa_name = fm\n    fb_name = fd\n    w = 0.4614171612406061')
    sequence=' '.join(fmt(t) for t in times)
    text=re.sub(r'  \[TimeStepper\].*?  \[\]',f"  [TimeStepper]\n    type = TimeSequenceStepper\n    time_sequence = '{sequence}'\n  []",text,flags=re.S)
    text=text.replace('results/ni_al_coarsening/coarsening',f'results/ni_al_lookup/{name}')
    text=text.replace('# Analytically eliminated KKS constraints for two local quadratic phase free energies.', '# Nested KKS; table-index lookup of Hermite coefficients, or quadratic numerical control.')
    (CASE/f'{name}.i').write_text(text)
print(json.dumps(report,indent=2))
