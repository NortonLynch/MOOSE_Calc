"""Independent controls for added mechanics and imported initial state."""
from pathlib import Path
import json
import numpy as np
import netCDF4

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/ni_al_elastic'

def nodal(path):
    with netCDF4.Dataset(path) as d:
        names=list(netCDF4.chartostring(d['name_nod_var'][:]))
        x,y=np.array(d['coordx'][:]),np.array(d['coordy'][:])
        order=np.lexsort((x,y))
        return np.array(d['time_whole'][:]),{s:np.array(d[f'vals_nod_var{i+1}'][:])[:,order] for i,s in enumerate(names)}

bt,b=nodal(ROOT/'results/ni_al_lookup/lookup.e')
zt,z=nodal(OUT/'zero_misfit.e')
errors={s:0. for s in ['c','eta','mu']}
matched=0
for i,t in enumerate(zt):
    j=int(np.argmin(abs(bt-t)))
    if abs(bt[j]-t)<1e-7:
        matched+=1
        for s in errors:errors[s]=max(errors[s],float(abs(z[s][i]-b[s][j]).max()))
assert matched>=4
assert max(errors.values())<1e-6,errors
_,ini=nodal(OUT/'elastic_initial.e')
_,run=nodal(OUT/'elastic.e')
initial_errors={s:float(abs(ini[s][-1]-run[s][0]).max()) for s in ['c','eta','mu','disp_x','disp_y']}
assert max(initial_errors.values())<1e-12,initial_errors
check={'zero_misfit_matching_saved_frames':matched,'zero_misfit_max_nodal_difference':errors,
       'steady_to_transient_initial_nodal_difference':initial_errors,'status':'PASS',
       'initial_auxiliary_output_note':'ComputeEigenstrainBase deliberately skips eigenstrain at step zero. Raw transient INITIAL stress/elastic-energy outputs therefore do not represent the imported relaxed state. Corrected analysis and seconds Exodus take ONLY those t=0 auxiliary/global mechanical quantities from elastic_initial.e final steady solution. Nodal c/eta/mu/disp fields are unchanged; their equality is checked above.',
       'source':'https://github.com/idaholab/moose/blob/master/modules/solid_mechanics/src/materials/ComputeEigenstrainBase.C'}
(OUT/'control_validation.json').write_text(json.dumps(check,indent=2)+'\n')
print(json.dumps(check,indent=2))
