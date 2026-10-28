from pathlib import Path
import json
import netCDF4
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def read(path):
    with netCDF4.Dataset(path) as d:
        names=list(netCDF4.chartostring(d['name_nod_var'][:]))
        order=np.lexsort((d['coordx'][:],d['coordy'][:]))
        return np.array(d['time_whole'][:]),{n:np.array(d[f'vals_nod_var{names.index(n)+1}'][:])[:,order] for n in ['c','eta','mu']}
t,a=read(ROOT/'results/ni_al_coarsening/coarsening.e')
s,b=read(ROOT/'results/ni_al_lookup/control_smoke.e')
report={}
for n in a:
    errors=[]
    for i,v in enumerate(s):
        j=np.argmin(abs(t-v))
        assert abs(t[j]-v)<1e-7
        errors.append(float(np.max(abs(a[n][j]-b[n][i]))))
    report[n+'_max_abs_difference']=max(errors)
    assert max(errors)<2e-7
report['frames_compared']=len(s)
report['last_time_star']=float(s[-1])
report['status']='PASS'
(ROOT/'results/ni_al_lookup/numerical_control.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
