"""Copy completed Exodus outputs to seconds/MPa units; preserve raw MOOSE files."""
from pathlib import Path
import json, shutil, hashlib
import numpy as np
import netCDF4

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/ni_al_elastic'
cal=json.loads((OUT/'time_calibration.json').read_text())
check=json.loads((OUT/'control_validation.json').read_text());assert check['status']=='PASS'
t0=cal['time_scale_s_per_model_unit'];f0=1000/7e-6

def names(d,kind):return list(netCDF4.chartostring(d[f'name_{kind}_var'][:]))
def rename(d,var,index,value):
    buf=np.zeros(d[var].shape[1],dtype='S1')
    b=value.encode();assert len(b)<=len(buf)
    buf[:len(b)]=np.frombuffer(b,dtype='S1');d[var][index,:]=buf

manifest={}
for key,src in [('no_elastic',ROOT/'results/ni_al_lookup/lookup.e'),('elastic',OUT/'elastic.e')]:
    dst=OUT/f'{key}_seconds.e'
    shutil.copy2(src,dst)
    with netCDF4.Dataset(dst,'r+') as d:
        ts=np.array(d['time_whole'][:]);assert abs(ts[-1]-60000)<1e-6
        en=names(d,'elem');gn=names(d,'glo')
        if key=='elastic':
            with netCDF4.Dataset(OUT/'elastic_initial.e') as ini:
                ien=names(ini,'elem');ign=names(ini,'glo')
                # Both meshes were generated with identical element IDs/order.
                assert np.array_equal(d['connect1'][:],ini['connect1'][:])
                assert np.array_equal(d['coordx'][:],ini['coordx'][:])
                for i,s in enumerate(en):
                    if s in ['energy_density','chemical_energy_density','elastic_energy_density'] or s.startswith('stress_'):
                        d[f'vals_elem_var{i+1}eb1'][0,:]=ini[f'vals_elem_var{ien.index(s)+1}eb1'][-1,:]
                for i,s in enumerate(gn):
                    if s in ['elastic_energy','total_energy'] or s.startswith('mean_stress_'):
                        d['vals_glo_var'][0,i]=ini['vals_glo_var'][-1,ign.index(s)]
        d['time_whole'][:]=ts*t0
        d['time_whole'].units='seconds'
        for coord in ['coordx','coordy']:d[coord].units='nm'
        for i,s in enumerate(en):
            v=d[f'vals_elem_var{i+1}eb1']
            if s.startswith('stress_'):
                v[:]=v[:]*f0/1e6;rename(d,'name_elem_var',i,s.replace('_scaled','_MPa'));v.units='MPa'
            elif s.endswith('energy_density'):
                v[:]=v[:]*f0;v.units='J/m^3'
        for i,s in enumerate(gn):
            if s.startswith('mean_stress_'):
                d['vals_glo_var'][:,i]=d['vals_glo_var'][:,i]*f0/1e6
                rename(d,'name_glo_var',i,s.replace('_scaled','_MPa'))
            elif s in ['total_energy','elastic_energy']:
                d['vals_glo_var'][:,i]=d['vals_glo_var'][:,i]*f0*1e-18
                rename(d,'name_glo_var',i,s+'_J_per_m')
            elif s in gn[:i]:
                # MOOSE emits both a ScalarVariable postprocessor and scalar output.
                rename(d,'name_glo_var',i,s+'_scalar')
        d.physical_time_basis='Literature-model calibration; see time_calibration.json, not direct experimental validation.'
        d.coordinate_unit='nm';d.energy_density_unit='J/m^3';d.stress_unit='MPa'
        d.initial_mechanics_note=check['initial_auxiliary_output_note'] if key=='elastic' else 'Original initial state.'
    with netCDF4.Dataset(dst) as d,netCDF4.Dataset(src) as raw:
        assert np.allclose(d['time_whole'][:],raw['time_whole'][:]*t0,rtol=0,atol=1e-10)
        nd=names(d,'nod');nr=names(raw,'nod')
        for s in ['c','eta','mu']:
            assert np.array_equal(d[f'vals_nod_var{nd.index(s)+1}'][:],raw[f'vals_nod_var{nr.index(s)+1}'][:])
        for v in d.variables.values():
            if v.dtype.kind in 'f':assert np.isfinite(v[:]).all(),v.name
    manifest[dst.name]={'bytes':dst.stat().st_size,'sha256':hashlib.sha256(dst.read_bytes()).hexdigest(),'source':str(src.relative_to(ROOT)),'time_unit':'s','coordinate_unit':'nm','stress_unit':'MPa','energy_density_unit':'J/m^3','integrated_energy_unit':'J/m of out-of-plane thickness','chemical_diffusion_potential_mu_unit':'1000 J/mol (unchanged)'}
(OUT/'physical_exodus_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps(manifest,indent=2))
