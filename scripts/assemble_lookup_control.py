"""Concatenate untouched saved prefix and freshly computed tail into native Exodus."""
from pathlib import Path
import argparse,json
import netCDF4,numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--output',default='results/ni_al_lookup');a=p.parse_args();OUT=ROOT/a.output
meta=json.loads((OUT/'matched_time_control.json').read_text());idx=meta['source_frame_index_zero_based']
target=OUT/'quadratic_matched.e';temp=OUT/'quadratic_assembled.tmp.e'
with netCDF4.Dataset(OUT/'quadratic_control.e') as old, netCDF4.Dataset(OUT/'quadratic_analytic_tail.e') as tail:
    old.set_auto_chartostring(False);tail.set_auto_chartostring(False)
    assert np.array_equal(old['coordx'][:],tail['coordx'][:]) and np.array_equal(old['coordy'][:],tail['coordy'][:])
    assert np.array_equal(old['connect1'][:],tail['connect1'][:])
    assert np.array_equal(old['name_nod_var'][:],tail['name_nod_var'][:])
    names=list(netCDF4.chartostring(old['name_nod_var'][:]))
    errors={n:float(abs(old[f'vals_nod_var{names.index(n)+1}'][idx]-tail[f'vals_nod_var{names.index(n)+1}'][0]).max()) for n in ['eta','c','mu']}
    assert max(errors.values())<1e-10,errors
    times=np.r_[old['time_whole'][:idx],tail['time_whole'][:]]
    assert np.all(np.diff(times)>0) and abs(times[-1]-60000)<1e-7
    with netCDF4.Dataset(temp,'w',format=old.data_model) as dst:
        dst.setncatts({k:old.getncattr(k) for k in old.ncattrs()})
        dst.setncattr('comparison_provenance','Exact prefix from quadratic_control.e and computed tail from quadratic_analytic_tail.e. No interpolation or field modification. See assembled_control.json.')
        for n,d in old.dimensions.items():dst.createDimension(n,None if d.isunlimited() else len(d))
        for n,var in old.variables.items():
            kw={'fill_value':var.getncattr('_FillValue')} if '_FillValue' in var.ncattrs() else {}
            new=dst.createVariable(n,var.datatype,var.dimensions,**kw)
            new.setncatts({k:var.getncattr(k) for k in var.ncattrs() if k!='_FillValue'});new.set_auto_chartostring(False)
            if 'time_step' in var.dimensions:
                axis=var.dimensions.index('time_step');sl=[slice(None)]*var.ndim;sl[axis]=slice(0,idx)
                values=np.concatenate([var[tuple(sl)],tail[n][:]],axis=axis)
                new[:]=values
                assert np.array_equal(new[:],values)
            else:new[:]=var[:]
if target.exists():
    backup=OUT/'quadratic_recovery_original_schedule.e';assert not backup.exists();target.rename(backup)
temp.rename(target)
qp=pd.read_csv(OUT/'quadratic_control.csv');tp=pd.read_csv(OUT/'quadratic_analytic_tail.csv')
combined=pd.concat([qp[qp.time<meta['checkpoint_time_star']-1e-7],tp],ignore_index=True)
reference=pd.read_csv(OUT/'lookup.csv').time.to_numpy()
assert all(min(abs(combined.time.to_numpy()-t))<1e-7 for t in reference)
extra=[float(t) for t in combined.time if min(abs(reference-t))>1e-7]
combined.to_csv(OUT/'quadratic_matched.csv',index=False,float_format='%.17g')
meta.update({'initial_field_max_abs_errors':errors,'saved_frames':len(times),'time_steps':len(combined)-1,'additional_control_times':extra,'exact_time_grid_match':not extra,'assembly':'Exact concatenation of original solver output fields; no resampling'})
(OUT/'assembled_control.json').write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps(meta,indent=2))
