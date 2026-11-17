"""Read-only mesh/output inventory for reliable Exodus conversion."""
from pathlib import Path
import json
import netCDF4,numpy as np
root=Path(__file__).resolve().parents[1]/'results/ni_al_elastic'
with netCDF4.Dataset(root/'elastic.e') as a,netCDF4.Dataset(root/'elastic_initial.e') as b:
    info={s:bool(np.array_equal(a[s][:],b[s][:])) for s in ['coordx','coordy','connect1']}
    info['global_names']=list(netCDF4.chartostring(a['name_glo_var'][:]))
    info['element_names']=list(netCDF4.chartostring(a['name_elem_var'][:]))
    info['string_dimension']=a['name_elem_var'].shape[1]
    print(json.dumps(info,indent=2))
