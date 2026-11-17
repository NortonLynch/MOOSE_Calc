"""Map the unchanged model clock using published fcc Ni-Al mobilities.

Engstrom & Agren (1996), DOI 10.1515/ijmr-1996-870205;
parameters independently transcribed from Table 3 of Zhang et al. (2010),
DOI 10.3139/146.110428 (previous-assessment rows [16,22]).
This is a literature-model calibration, not a 1000 K experimental measurement.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.interpolate import CubicHermiteSpline

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results/ni_al_elastic'
T, R, l0, G0, Mstar = 1000., 8.31446261815324, 1e-9, 1000., 0.002873786184560299
meta = json.loads((ROOT/'thermodynamics/ni_al_1000K/metadata.json').read_text())
c = meta['validation']['equilibrium']['common_tangent_x_Al_from_interpolation'][0]
tab = pd.read_csv(ROOT/'thermodynamics/ni_al_1000K/gamma_1000K.csv')
g = CubicHermiteSpline(tab.x_Al, tab.Gm_reduced_J_mol, tab.dGm_reduced_dx_J_mol)
curvature = float(g(c, 2))

def evaluate(x, model='engstrom_1996'):
    if model == 'engstrom_1996':
        # U_i = -Q_i + R T ln(D0_i); all U coefficients in J/mol atoms.
        ua = x*(-142000-72.12*T)+(1-x)*(-284000-59.83*T)+x*(1-x)*(-41300-91.2*T)
        un = x*(-145900-64.26*T)+(1-x)*(-287000-69.8*T)+x*(1-x)*(-113000+65.5*T)
    elif model == 'zhang_2010':
        ua = x*(-123111.6-97.34*T)+(1-x)*(-268381-71.04*T)+x*(1-x)*(-308067.5+111.52*T)
        un = x*(-144600-64.85*T)+(1-x)*(-271377.6-81.79*T)+x*(1-x)*(-29571.8)
    else:
        raise ValueError(model)
    da, dn = np.exp(ua/(R*T)), np.exp(un/(R*T))
    phi = x*(1-x)*g(x, 2)/(R*T)
    return da, dn, phi, ((1-x)*da+x*dn)*phi

da, dn, phi, D = map(float, evaluate(c))
t0 = Mstar*(curvature/G0)*l0*l0/D
_, _, _, Dalt = map(float, evaluate(c, 'zhang_2010'))
result = {
    'status':'literature_model_calibration_not_direct_experimental_validation',
    'primary_model':'Engstrom & Agren 1996 fcc mobilities, reproduced in Zhang et al. 2010 Table 3 rows [16,22]',
    'primary_doi':'https://doi.org/10.1515/ijmr-1996-870205',
    'parameter_source_doi':'https://doi.org/10.3139/146.110428',
    'independent_parameter_crosscheck':'NIST C.E. Campbell, Re-assessment of gamma prime (L12-Ni3Al), July 9 2008, Disordered FCC table; same Engstrom-Agren parameters (exponential prefactors versus rounded linear T coefficients).',
    'independent_parameter_crosscheck_url':'https://materialsdata.nist.gov/bitstream/handle/11115/51/Re-assessment-l12-7-08.pdf?isAllowed=y&sequence=4',
    'parameter_source_pdf':'https://www.researchgate.net/profile/Lijun-Zhang-46/publication/273970625_Atomic_mobilities_and_diffusivities_in_the_fcc_L12_and_B2_phases_of_the_Ni-Al_system/links/5526789b0cf2628d5afdfd6b/Atomic-mobilities-and-diffusivities-in-the-fcc-L12-and-B2-phases-of-the-Ni-Al-system.pdf',
    'temperature_K':T,'calibration_gamma_Al_atomic_fraction':c,
    'thermodynamics':'Thermo-Calc Educational 2026b NIDEMO 2.0 exported gamma Hermite table; retained unchanged',
    'gamma_Gm_second_derivative_J_mol':curvature,
    'Al_tracer_diffusivity_m2_s':da,'Ni_tracer_diffusivity_m2_s':dn,
    'thermodynamic_factor':phi,'gamma_chemical_interdiffusivity_m2_s':D,
    'dimensionless_mobility':Mstar,'length_scale_m':l0,'molar_energy_scale_J_mol':G0,
    'time_scale_s_per_model_unit':t0,'end_time_s':60000*t0,'end_time_hours':60000*t0/3600,
    'formula_D':'Dchem = ((1-c)*DAl_tracer+c*DNi_tracer)*c*(1-c)*Gm_second_derivative/(R*T)',
    'formula_time':'t_seconds=t_star*Mstar*(Gm_second_derivative/G0)*length_scale_m^2/Dchem',
    'formula_rate':'K_nm3_s=K_nm3_per_model_time/time_scale_s_per_model_unit; K_m3_s=K_nm3_s*1e-27',
    'alternative_zhang_2010_Dchem_m2_s':Dalt,
    'alternative_zhang_2010_time_scale_s':Mstar*(curvature/G0)*l0*l0/Dalt,
    'scope':['1000 K evaluation of published high-temperature mobility assessment; extrapolation uncertainty is not quantified.',
             'Darken form with constant molar volume; vacancy-wind/correlation correction omitted as in the mobility framework.',
             'Mobilities and current NIDEMO thermodynamic factor are combined; this is not an exact reproduction of the original assessed thermodynamic database.',
             'Constant original Mstar and Lstar retained for comparison. This rescales the entire simulation, not a spatially composition-dependent mobility simulation.',
             'The original single constant mobility is used across both phases; gamma-prime ordered-phase mobility has not been independently calibrated. The time scale anchors matrix diffusion only.',
             'Elastic stress does not independently modify diffusivity in this model.',
             'Janssen temperature range is 1000-1300 degrees C, not 1000 K; no single measurement at 1000 K has been asserted.'],
}
OUT.mkdir(exist_ok=True)
(OUT/'time_calibration.json').write_text(json.dumps(result,indent=2)+'\n')
model_path=OUT/'model_parameters.json'
if model_path.exists():
    model=json.loads(model_path.read_text())
    model['time_mapping']='Literature-model calibration in time_calibration.json: Engstrom-Agren fcc mobilities + current gamma Hermite curvature; not direct experimental calibration.'
    model['numerical_solver']='PJFNK; full SMP preconditioning, MUMPS, automatic scaling recalculated every step; implicit Euler; original time sequence with cutback on failure; production 8 MPI ranks.'
    model_path.write_text(json.dumps(model,indent=2)+'\n')
xs=np.linspace(.10,.16,121)
aa,nn,pp,dd=evaluate(xs)
pd.DataFrame({'gamma_x_Al':xs,'Al_tracer_m2_s':aa,'Ni_tracer_m2_s':nn,'thermodynamic_factor':pp,'chemical_interdiffusivity_m2_s':dd}).to_csv(OUT/'diffusion_calibration_curve.csv',index=False)
print(json.dumps(result,indent=2))
