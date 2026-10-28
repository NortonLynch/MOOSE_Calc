"""Drop superseded Exodus records beyond the final solver time after recovery.

No field values are changed. Preserve the original file for audit.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import netCDF4

out = Path(__file__).resolve().parents[1] / 'results/ni_al_coarsening'
source = out / 'coarsening.e'
final_time = float(pd.read_csv(out/'coarsening.csv').time.iloc[-1])
with netCDF4.Dataset(source) as old:
    old.set_auto_chartostring(False)
    times = np.asarray(old['time_whole'][:])
    keep = np.flatnonzero(times <= final_time + 1e-8)
    assert len(keep)>0 and np.array_equal(keep,np.arange(len(keep)))
    assert np.isclose(times[keep[-1]],final_time)
    assert np.all(np.diff(times[keep])>0)
    if len(keep)==len(times):
        print('No superseded trailing records found.')
        raise SystemExit(0)
    temporary = out/'coarsening_trimmed.tmp.e'
    with netCDF4.Dataset(temporary,'w',format=old.data_model) as new:
        new.setncatts({key:old.getncattr(key) for key in old.ncattrs()})
        for name,dim in old.dimensions.items():
            new.createDimension(name,None if dim.isunlimited() else len(dim))
        for name,var in old.variables.items():
            kwargs={}
            if '_FillValue' in var.ncattrs():
                kwargs['fill_value']=var.getncattr('_FillValue')
            dst=new.createVariable(name,var.datatype,var.dimensions,**kwargs)
            dst.setncatts({key:var.getncattr(key) for key in var.ncattrs() if key!='_FillValue'})
            dst.set_auto_chartostring(False)
            if 'time_step' in var.dimensions:
                axis=var.dimensions.index('time_step')
                sl=[slice(None)]*var.ndim
                sl[axis]=slice(0,len(keep))
                dst[:]=var[tuple(sl)]
            else:
                dst[:]=var[:]
    removed=times[len(keep):].tolist()
backup=out/'coarsening_with_superseded_restart_tail.e'
assert not backup.exists(), 'Do not overwrite an earlier audit file'
source.rename(backup)
temporary.rename(source)
(out/'restart_tail_cleanup.json').write_text(json.dumps({'original_file':backup.name,'delivered_file':source.name,'final_solver_time':final_time,'kept_records':len(keep),'removed_superseded_times':removed,'operation':'Copy all original metadata and fields unchanged; remove only records after the final accepted solver time.'},indent=2)+'\n')
print(f'Kept {len(keep)} frames through {final_time}; removed {removed}')
