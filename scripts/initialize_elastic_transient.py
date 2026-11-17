"""Initialize mechanics at equilibrium while keeping original eta/c/mu fields."""
from pathlib import Path
import json,re
import pandas as pd
import netCDF4,numpy as np
ROOT=Path(__file__).resolve().parents[1]
CASE=ROOT/'examples/ni_al_elastic';OUT=ROOT/'results/ni_al_elastic'
with netCDF4.Dataset(OUT/'elastic_initial.e') as d:
    n=len(d['time_whole'])
    names=list(netCDF4.chartostring(d['name_nod_var'][:]))
    fields={s:np.array(d[f'vals_nod_var{names.index(s)+1}'][-1]) for s in ['c','eta']}
with netCDF4.Dataset(ROOT/'results/ni_al_lookup/lookup.e') as d:
    names=list(netCDF4.chartostring(d['name_nod_var'][:]))
    errors={s:float(np.max(np.abs(fields[s]-d[f'vals_nod_var{names.index(s)+1}'][0]))) for s in fields}
assert max(errors.values())<1e-13,errors
last=pd.read_csv(OUT/'elastic_initial.csv').iloc[-1]
assert max(abs(last[f'mean_stress_{s}_scaled']) for s in ['xx','yy','xy'])<1e-9
values=[float(last[f'global_strain_{i}']) for i in range(3)]
text=(CASE/'elastic.i').read_text()
text+='\n[ICs]\n  [global_strain_ic]\n    type = ScalarComponentIC\n    variable = global_strain\n    values = \''+' '.join(format(v,'.17g') for v in values)+"'\n  []\n[]\n"
var='[Variables]\n'
for direction in ['x','y']:
    var+=f'''  [disp_{direction}]
    [InitialCondition]
      type = SolutionIC
      solution_uo = relaxed_initial
      from_variable = disp_{direction}
    []
  []
'''
text=text.replace('[Variables]\n',var,1).replace('add_variables = true','add_variables = false')
text=text.replace('[UserObjects]\n',f'''[UserObjects]
  [relaxed_initial]
    type = SolutionUserObject
    mesh = ../../results/ni_al_elastic/elastic_initial.e
    timestep = {n}
    system_variables = 'disp_x disp_y'
  []
''',1)
(CASE/'elastic_relaxed.i').write_text(text)
result={'initial_eta_c_max_error_against_lookup':errors,'global_inplane_strain':values,'initial_elastic_energy_scaled_integral':float(last.elastic_energy),'initial_average_stress_scaled':{s:float(last[f'mean_stress_{s}_scaled']) for s in ['xx','yy','xy','zz']},'mechanical_relaxation_does_not_advance_diffusion_time':True}
(OUT/'initialization_check.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
