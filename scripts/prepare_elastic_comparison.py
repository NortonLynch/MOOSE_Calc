"""Add coherent elasticity to the completed lookup case without changing chemistry."""
from pathlib import Path
import hashlib, json, re, shutil

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT / 'examples/ni_al_elastic'
OUT = ROOT / 'results/ni_al_elastic'
CASE.mkdir(exist_ok=True)
OUT.mkdir(exist_ok=True)
source = ROOT / 'examples/ni_al_lookup/lookup.i'
text = source.read_text()
text = text.replace('# Ni-14 at.% Al, 1000 K. Demonstration with REFERENCE kinetics, not calibrated time.', '# Ni-14 at.% Al, 1000 K. Literature-model clock calibration: see time_calibration.json.')
text = text.replace('# Coordinates: 1 unit = 1 nm. Time: t*; reference t0 = 1 s.', '# Coordinates: 1 unit = 1 nm. Solver time remains t*; convert using scripts/calibrate_elastic_time.py.')

def append_block(text, name, body):
    pattern = r'(^\[' + re.escape(name) + r'\]\n)(.*?)(^\[\])'
    text, n = re.subn(pattern, lambda m: m[1]+m[2]+body+m[3], text,
                     count=1, flags=re.M|re.S)
    assert n == 1, name
    return text

# Dimensional stiffness / energy-density scale. Coordinates and displacements: nm.
f0 = 1000 / 7e-6
cg = [199,143,143,199,143,199,99,99,99]
cp = [207.8,148.8,148.8,207.8,148.8,207.8,102.3,102.3,102.3]
tensor = lambda values: ' '.join(format(v*1e9/f0,'.16g') for v in values)
text = text.replace('results/ni_al_lookup/lookup', 'results/ni_al_elastic/elastic')
text = text.replace('solve_type = NEWTON', 'solve_type = PJFNK')
text = text.replace('  solve_type = PJFNK', '  solve_type = PJFNK\n  automatic_scaling = true\n  compute_scaling_once = false',1)
text = text.replace('\n  l_max_its = 30', '\n  l_max_its = 100\n  l_tol = 1e-5')
text = append_block(text, 'Mesh', '''  [anchor]
    type = ExtraNodesetGenerator
    input = domain
    coord = '0 0'
    new_boundary = anchor
  []
''')
text = append_block(text, 'Variables', '''  [global_strain]
    family = SCALAR
    order = THIRD
  []
''')
text = text.replace("auto_direction = 'x y'", "auto_direction = 'x y'\n      variable = 'c eta mu disp_x disp_y'")
for direction in ('x','y'):
    text = append_block(text, 'BCs', f'''  [anchor_{direction}]
    type = DirichletBC
    variable = disp_{direction}
    boundary = anchor
    value = 0
  []
''')
text = append_block(text, 'Kernels', '''  [elastic_driving_force]
    type = AllenCahn
    variable = eta
    f_name = f_el
    mob_name = L
  []
''')
text = append_block(text, 'Materials', f'''  [elasticity_gamma]
    type = ComputeElasticityTensor
    base_name = gamma
    fill_method = symmetric9
    C_ijkl = '{tensor(cg)}'
  []
  [elasticity_gamma_prime]
    type = ComputeElasticityTensor
    base_name = gammaprime
    fill_method = symmetric9
    C_ijkl = '{tensor(cp)}'
  []
  [matrix_weight]
    type = DerivativeParsedMaterial
    property_name = one_minus_h
    coupled_variables = eta
    expression = '1-eta^3*(6*eta^2-15*eta+10)'
    derivative_order = 2
  []
  [elasticity_tensor]
    type = CompositeElasticityTensor
    coupled_variables = eta
    tensors = 'gamma gammaprime'
    weights = 'one_minus_h h'
  []
  [eigenstrain]
    type = ComputeVariableEigenstrain
    eigen_base = '0.0047 0.0047 0.0047 0 0 0'
    prefactor = h
    args = eta
    eigenstrain_name = eigenstrain
  []
  [stress]
    type = ComputeLinearElasticStress
  []
  [elastic_free_energy]
    type = ElasticEnergyMaterial
    property_name = f_el
    coupled_variables = eta
    derivative_order = 2
  []
  [global_strain_material]
    type = ComputeGlobalStrain
    scalar_global_strain = global_strain
    global_strain_uo = global_strain_uo
  []
''')
extra_fields = ['chemical_energy_density','elastic_energy_density','stress_xx_scaled','stress_yy_scaled','stress_xy_scaled','stress_zz_scaled']
for name in extra_fields:
    text = append_block(text, 'AuxVariables', f'''  [{name}]
    family = MONOMIAL
    order = CONSTANT
  []
''')
text = text.replace('type = KKSGlobalFreeEnergy\n    variable = energy_density', 'type = KKSGlobalFreeEnergy\n    variable = chemical_energy_density')
text = append_block(text, 'AuxKernels', '''  [elastic_energy]
    type = MaterialRealAux
    variable = elastic_energy_density
    property = f_el
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [combined_energy]
    type = ParsedAux
    variable = energy_density
    coupled_variables = 'chemical_energy_density elastic_energy_density'
    expression = 'chemical_energy_density+elastic_energy_density'
    execute_on = 'INITIAL TIMESTEP_END'
  []
''')
for name,i,j in [('xx',0,0),('yy',1,1),('xy',0,1),('zz',2,2)]:
    text = append_block(text,'AuxKernels',f'''  [stress_{name}]
    type = RankTwoAux
    variable = stress_{name}_scaled
    rank_two_tensor = stress
    index_i = {i}
    index_j = {j}
    execute_on = 'INITIAL TIMESTEP_END'
  []
''')
    text = append_block(text,'Postprocessors',f'''  [mean_stress_{name}_scaled]
    type = ElementAverageValue
    variable = stress_{name}_scaled
    execute_on = 'INITIAL TIMESTEP_END'
  []
''')
text = append_block(text, 'Postprocessors', '''  [elastic_energy]
    type = ElementIntegralVariablePostprocessor
    variable = elastic_energy_density
    execute_on = 'INITIAL TIMESTEP_END'
  []
''')
for i in range(3):
    text = append_block(text, 'Postprocessors', f'''  [global_strain_{i}]
    type = ScalarVariable
    variable = global_strain
    component = {i}
    execute_on = 'INITIAL TIMESTEP_END'
  []
''')
text += '''
[GlobalParams]
  displacements = 'disp_x disp_y'
[]
[Physics/SolidMechanics/QuasiStatic]
  [mechanics]
    strain = SMALL
    incremental = false
    add_variables = true
    eigenstrain_names = eigenstrain
    global_strain = global_strain
  []
[]
[ScalarKernels]
  [global_strain]
    type = GlobalStrain
    variable = global_strain
    global_strain_uo = global_strain_uo
  []
[]
[UserObjects]
  [global_strain_uo]
    type = GlobalStrainUserObject
    applied_stress_tensor = '0 0 0 0 0 0'
    execute_on = 'INITIAL LINEAR NONLINEAR'
  []
[]
'''
text = '# Lookup chemical thermodynamics + coherent cubic elasticity at 1000 K.\n# Plane strain; periodic fluctuations and zero in-plane average stress.\n' + text
(CASE/'elastic.i').write_text(text)
# A zero-misfit limit checks that adding mechanics recovers the original solution.
zero = text.replace("eigen_base = '0.0047 0.0047 0.0047 0 0 0'", "eigen_base = '0 0 0 0 0 0'").replace('results/ni_al_elastic/elastic','results/ni_al_elastic/zero_misfit')
(CASE/'zero_misfit.i').write_text(zero)
# Fixed phase fields for initial mechanical relaxation.
equilibrium = re.sub(r'^\[Kernels\].*?^\[\]\n', '', text, flags=re.M|re.S)
fields = re.search(r'^\[Variables\]\n(.*?)^\[\]',equilibrium,re.M|re.S)[1]
phase_fields = fields[:fields.index('  [global_strain]')]
equilibrium = equilibrium.replace(phase_fields,'',1)
equilibrium = append_block(equilibrium,'AuxVariables',phase_fields)
equilibrium = equilibrium.replace("variable = 'c eta mu disp_x disp_y'", "variable = 'disp_x disp_y'")
equilibrium = re.sub(r'^\[Executioner\].*?^\[\]', '''[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu mumps'
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-11
  nl_max_its = 50
  l_tol = 1e-6
  l_max_its = 100
[]''', equilibrium,flags=re.M|re.S)
equilibrium = equilibrium.replace('results/ni_al_elastic/elastic','results/ni_al_elastic/elastic_initial')
equilibrium = equilibrium.replace('time_step_interval = 5','time_step_interval = 1')
(CASE/'elastic_initial.i').write_text(equilibrium)
metadata = {
    'temperature_K':1000,'nominal_Al_atomic_fraction':0.14,
    'chemical_model':'unchanged lookup nested KKS, inherited from examples/ni_al_lookup/lookup.i',
    'input_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
    'elastic_data_source':'User attachment plasti.txt; values accepted as supplied, not independently calibrated',
    'gamma_c11_c12_c44_GPa':[199,143,99],
    'gamma_prime_c11_c12_c44_GPa':[207.8,148.8,102.3],
    'misfit':0.0047,'energy_density_scale_J_m3':f0,
    'GPa_to_dimensionless_factor':1e9/f0,
    'attachment_unit_correction':'199 GPa = 1242.0603 eV/nm^3, not 1.242 eV/nm^3; use the supplied GPa values',
    'boundary':'Periodic c, eta, mu, displacement fluctuations; one point removes translation; global in-plane strain enforces zero volume-average xx, yy, xy stress; epsilon_zz=0 (plane strain)',
    'orientation':'cubic [100] || x, [010] || y, [001] || z, common orientation',
    'energy':'C(eta)=(1-h)C_gamma+h C_gamma_prime; epsilon0=0.0047 h I; f_el=0.5 (epsilon-epsilon0):C:(epsilon-epsilon0)',
    'coupling':'AllenCahn derivative of f_el added to eta equation; stress equilibrium solved simultaneously; f_el has no explicit c dependence, so KKS chemical partition is unchanged at fixed eta',
    'time_mapping':'pending verification of user-provided diffusion literature; all simulations retain original dimensionless M and L',
}
(OUT/'model_parameters.json').write_text(json.dumps(metadata,indent=2)+'\n')
print(json.dumps(metadata,indent=2))
