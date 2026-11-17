"""Analyze measured MOOSE fields; never impose an LSW growth law on the output."""
from pathlib import Path
import argparse, json, base64, shutil, hashlib
import numpy as np
import pandas as pd
import netCDF4
from scipy.ndimage import label
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/ni_al_elastic'
parser=argparse.ArgumentParser()
parser.add_argument('--partial',action='store_true')
parser.add_argument('--no-gif',action='store_true')
args=parser.parse_args()
cal=json.loads((OUT/'time_calibration.json').read_text())
t0=cal['time_scale_s_per_model_unit']
f0=1000/7e-6
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})

def read(path):
    with netCDF4.Dataset(path) as d:
        names=list(netCDF4.chartostring(d['name_nod_var'][:]))
        x,y=np.array(d['coordx'][:]),np.array(d['coordy'][:])
        order=np.lexsort((x,y)); n=len(np.unique(x))
        t=np.array(d['time_whole'][:]); assert np.all(np.diff(t)>0)
        nod={s:np.array(d[f'vals_nod_var{names.index(s)+1}'][:])[:,order].reshape(-1,n,n) for s in ['c','eta','mu']}
        enames=list(netCDF4.chartostring(d['name_elem_var'][:]))
        conn=np.array(d['connect1'][:])-1
        eo=np.lexsort((x[conn].mean(axis=1),y[conn].mean(axis=1)))
        wanted=['c_gamma','c_gamma_prime','stress_xx_scaled','stress_yy_scaled','stress_zz_scaled','elastic_energy_density']
        elem={s:np.array(d[f'vals_elem_var{enames.index(s)+1}eb1'][:])[:,eo].reshape(-1,n-1,n-1) for s in wanted if s in enames}
        for a in list(nod.values())+list(elem.values()):assert np.isfinite(a).all()
        return t,nod,elem

def history(t,nod):
    eta=nod['eta']; h=eta**3*(6*eta**2-15*eta+10)
    n=eta.shape[-1]; m=(n-1)//2; dx=96/(n-1)
    w=np.ones((n,n));w[[0,-1],:]*=.5;w[:,[0,-1]]*=.5
    slices=[(slice(0,m+1),slice(0,m+1)),(slice(0,m+1),slice(m,n)),(slice(m,n),slice(0,m+1)),(slice(m,n),slice(m,n))]
    rows=[]
    for k,tt in enumerate(t):
        row={'time_star':tt,'time_s':tt*t0,'time_h':tt*t0/3600,'mean_c':float((nod['c'][k]*w).sum()*dx**2/96**2),'gamma_prime_area_fraction':float((h[k]*w).sum()*dx**2/96**2)}
        radii=[]
        for i,(ys,xs) in enumerate(slices,1):
            q=h[k,ys,xs];qw=np.ones_like(q);qw[[0,-1],:]*=.5;qw[:,[0,-1]]*=.5
            r=float(np.sqrt(max(0,(q*qw).sum()*dx**2)/np.pi))
            present=int(eta[k,ys,xs].max()>=.5)
            row[f'R{i}_nm']=r;row[f'present{i}']=present
            if present:radii.append(r)
        # Four well-separated seeds stay away from domain edges; verify topology.
        count=label(eta[k]>=.5)[1]
        assert count==len(radii),f'Particle assignment invalid at {tt}: {count}, {len(radii)}'
        assert not (eta[k,0]>=.5).any() and not (eta[k,-1]>=.5).any()
        assert not (eta[k,:,0]>=.5).any() and not (eta[k,:,-1]>=.5).any()
        row['count']=count
        row['mean_R_nm']=float(np.mean(radii)) if radii else np.nan
        row['mean_R_cubed_nm3']=row['mean_R_nm']**3
        row['mean_of_R_cubes_nm3']=float(np.mean(np.array(radii)**3)) if radii else np.nan
        rows.append(row)
    return pd.DataFrame(rows)

def fit(tr,lo,hi,name):
    use=tr[(tr.time_star>=lo)&(tr.time_star<=hi)]
    if len(use)<4:return None
    x=use.time_s.to_numpy();y=use.mean_R_cubed_nm3.to_numpy()
    w=np.r_[np.diff(x)[0],x[2:]-x[:-2],np.diff(x)[-1]]/2
    xx=x-np.average(x,weights=w); yy=y-np.average(y,weights=w)
    slope=float(np.dot(w*xx,yy)/np.dot(w*xx,xx))
    intercept=float(np.average(y-slope*x,weights=w))
    pred=intercept+slope*x; sse=float(np.dot(w,(y-pred)**2));sst=float(np.dot(w,yy**2))
    return {'window':name,'lo_star':float(use.time_star.iloc[0]),'hi_star':float(use.time_star.iloc[-1]),'samples':len(use),'K_nm3_s':slope,'K_m3_s':slope*1e-27,'intercept_nm3':intercept,'R_squared_time_weighted':1-sse/sst if sst>0 else None,'RMSE_nm3':float(np.sqrt(sse/w.sum())),'particle_counts':sorted(map(int,use['count'].unique())),'crosses_particle_loss':bool(use['count'].nunique()>1),'interpretation':'descriptive effective slope; not a validated 3D LSW rate constant'}

def shape_p4(eta):
    """Directional eta=0.5 radii, using the weighted centroid of the P4 quadrant."""
    x=np.linspace(0,96,eta.shape[0]);mid=(len(x)-1)//2
    yy,xx=np.meshgrid(x[mid:],x[mid:],indexing='ij')
    q=eta[mid:,mid:];h=q**3*(6*q**2-15*q+10)
    w=np.ones_like(h);w[[0,-1],:]*=.5;w[:,[0,-1]]*=.5;h=h*w
    cx=float((h*xx).sum()/h.sum());cy=float((h*yy).sum()/h.sum())
    interp=RegularGridInterpolator((x,x),eta)
    radii=[]
    for theta in np.arange(8)*np.pi/4:
        fn=lambda r:float(interp([[cy+r*np.sin(theta),cx+r*np.cos(theta)]])[0])-.5
        radii.append(brentq(fn,0,22))
    axis=float(np.mean(radii[::2]));diagonal=float(np.mean(radii[1::2]))
    return {'centroid_nm':[cx,cy],'axis_radius_nm':axis,'diagonal_radius_nm':diagonal,'diagonal_axis_ratio':diagonal/axis,'definition':'eta=0.5 directional radii; diagonal mean / x-y-axis mean; 1 for a circle; P4 only'}

runs={};metrics={'time_calibration':cal,'fits':{},'status':'PARTIAL' if args.partial else 'COMPLETE','radius_definition':'R_i=sqrt(integral_over_seed_quadrant h(eta) dA/pi), present if max eta >=0.5; number average over surviving particles, then cubed; 2D area-equivalent radius','regression':'time-weighted least squares (trapezoidal weights) prevents adaptive output clusters from dominating; no forced zero intercept'}
for key,path,title,color in [('no_elastic',ROOT/'results/ni_al_lookup/lookup.e','Without elasticity','#277da8'),('elastic',OUT/'elastic.e','With elasticity','#cf5b36')]:
    t,nod,elem=read(path);tr=history(t,nod);post=pd.read_csv(path.with_suffix('.csv'))
    if not args.partial:assert abs(t[-1]-60000)<1e-6 and abs(post.time.iloc[-1]-60000)<1e-6
    if key=='elastic':
        initial=pd.read_csv(OUT/'elastic_initial.csv').iloc[-1]
        for s in ['elastic_energy','total_energy','mean_stress_xx_scaled','mean_stress_yy_scaled','mean_stress_xy_scaled','mean_stress_zz_scaled']:
            post.loc[0,s]=initial[s]
        _,_,ie=read(OUT/'elastic_initial.e')
        for s in ie:
            if s.startswith('stress_') or s=='elastic_energy_density':elem[s][0]=ie[s][-1]
    tr.to_csv(OUT/f'{key}_particle_history.csv',index=False)
    post['time_s']=post.time*t0;post['time_h']=post.time*t0/3600
    post.to_csv(OUT/f'{key}_diagnostics_physical.csv',index=False)
    m={'saved_frames':len(t),'accepted_steps':len(post)-1,'end_time_star':float(t[-1]),'end_time_h':float(t[-1]*t0/3600),'initial_mean_R_nm':float(tr.mean_R_nm.iloc[0]),'final_mean_R_nm':float(tr.mean_R_nm.iloc[-1]),'initial_mean_R_cubed_nm3':float(tr.mean_R_cubed_nm3.iloc[0]),'final_mean_R_cubed_nm3':float(tr.mean_R_cubed_nm3.iloc[-1]),'final_count':int(tr['count'].iloc[-1]),'radii_final_nm':tr[[f'R{j}_nm' for j in range(1,5)]].iloc[-1].tolist(),'max_mass_error':float(abs(post.mean_c-.14).max()),'max_nodal_integrated_mass_error':float(abs(tr.mean_c-.14).max()),'initial_energy_integral_scaled':float(post.total_energy.iloc[0]),'final_energy_integral_scaled':float(post.total_energy.iloc[-1]),'maximum_energy_increase_scaled':float(max(0,np.diff(post.total_energy).max())),'dissolution_brackets_star':{}}
    assert m['max_mass_error']<1e-8 and m['max_nodal_integrated_mass_error']<1e-8
    assert m['maximum_energy_increase_scaled']<1e-5,m
    if key=='elastic':
        m['max_average_inplane_stress_MPa']=float(abs(post[['mean_stress_xx_scaled','mean_stress_yy_scaled','mean_stress_xy_scaled']]).to_numpy().max()*f0/1e6)
        assert m['max_average_inplane_stress_MPa']<1e-3,m
    center_eta=(nod['eta'][:,:-1,:-1]+nod['eta'][:,1:,:-1]+nod['eta'][:,:-1,1:]+nod['eta'][:,1:,1:])/4
    core=pd.DataFrame({'time_star':t,'time_s':t*t0,'time_h':t*t0/3600})
    for s in ['c_gamma','c_gamma_prime']:
        mask=center_eta<.01 if s=='c_gamma' else center_eta>.99
        core[s+'_core_mean']=[float(elem[s][i][mask[i]].mean()) if mask[i].any() else np.nan for i in range(len(t))]
        m[s+'_final_core_mean']=float(core[s+'_core_mean'].iloc[-1])
        m[s+'_output_range']=[float(elem[s].min()),float(elem[s].max())]
        bounds=(.05,.30) if s=='c_gamma' else (.18,.30)
        assert bounds[0]<=m[s+'_output_range'][0]<=m[s+'_output_range'][1]<=bounds[1]
    m['final_area_fraction']=float(tr.gamma_prime_area_fraction.iloc[-1])
    m['P4_final_shape']=shape_p4(nod['eta'][-1])
    core.to_csv(OUT/f'{key}_core_compositions.csv',index=False)
    for j in range(1,5):
        lost=np.flatnonzero(tr[f'present{j}'].to_numpy()==0)
        m['dissolution_brackets_star'][str(j)]=[float(t[lost[0]-1]),float(t[lost[0]])] if len(lost) else None
    metrics[key]=m
    windows=[(5000,30000,'5000-30000'),(10000,40000,'10000-40000'),(50000,60000,'50000-60000'),(5000,60000,'5000-60000')]
    # Also report every sufficiently sampled constant-count interval, without cherry-picking.
    boundaries=np.r_[0,np.flatnonzero(np.diff(tr['count'])!=0)+1,len(tr)]
    for a,b in zip(boundaries[:-1],boundaries[1:]):
        windows.append((max(5000,t[a]),t[b-1],f'constant_N{int(tr["count"].iloc[a])}_{a}_{b-1}'))
    fits=[fit(tr,lo,hi,name) for lo,hi,name in windows if t[-1]>=hi]
    metrics['fits'][key]=[f for f in fits if f is not None]
    runs[key]=(t,nod,elem,tr,post,title,color)
assert np.array_equal(runs['no_elastic'][1]['eta'][0],runs['elastic'][1]['eta'][0])
assert np.array_equal(runs['no_elastic'][1]['c'][0],runs['elastic'][1]['c'][0])
end=min(runs['elastic'][0][-1],60000)
bt=runs['no_elastic'][4].time.to_numpy();et=runs['elastic'][4].time.to_numpy()
metrics['integration_time_grid_comparison']={
    'identical_accepted_grid':bool(len(bt)==len(et) and np.allclose(bt,et,rtol=0,atol=1e-7)),
    'same_prescribed_schedule':'Both start from the original baseline time sequence. Nonlinear failures may subdivide accepted steps, so this does not guarantee identical accepted grids near dissolution.',
    'identical_through_30000':bool(len(bt[bt<=30000])==len(et[et<=30000]) and np.allclose(bt[bt<=30000],et[et<=30000],rtol=0,atol=1e-7)) if end>=30000 else None,
}

fig,axs=plt.subplots(2,3,figsize=(11.4,7.4),layout='constrained')
for row,(key,run) in enumerate(runs.items()):
    t,nod,_,tr,_,title,color=run
    for col,target in enumerate([0,end/2,end]):
        k=int(np.argmin(abs(t-target)));ax=axs[row,col]
        im=ax.imshow(nod['eta'][k],origin='lower',extent=[0,96,0,96],vmin=0,vmax=1,cmap='Blues')
        ax.contour(np.linspace(0,96,97),np.linspace(0,96,97),nod['eta'][k],levels=[.5],colors='k',linewidths=.7)
        ax.set_title(f'{title}\nt = {t[k]*t0/3600:.3f} h; N = {tr["count"].iloc[k]}')
        ax.set(xlabel='x (nm)',ylabel='y (nm)')
fig.colorbar(im,ax=axs,label=r'Order parameter $\eta$',shrink=.8)
fig.suptitle('Ni-14 at.% Al, 1000 K | identical initial state and chemical free energy')
fig.savefig(OUT/'morphology_comparison.png');plt.close(fig)

fig,axs=plt.subplots(2,1,figsize=(10,7),sharex=True,layout='constrained',gridspec_kw={'height_ratios':[3,1]})
for key,run in runs.items():
    _,_,_,tr,_,title,color=run;use=tr[tr.time_star<=end+1e-8]
    axs[0].plot(use.time_h,use.mean_R_cubed_nm3,label=title,color=color,lw=1.8)
    axs[1].step(use.time_h,use['count'],where='post',color=color)
    for f in metrics['fits'][key]:
        if f['window']=='5000-60000' and end>=60000:
            xx=np.array([f['lo_star'],f['hi_star']])*t0
            axs[0].plot(xx/3600,f['intercept_nm3']+f['K_nm3_s']*xx,'--',color=color,lw=1,label=f'Whole-window fit: K={f["K_nm3_s"]:.4g} nm³/s, R²={f["R_squared_time_weighted"]:.3f}')
axs[0].set(ylabel=r'$\langle R\rangle^3$ (nm$^3$)',title='Cube of the number-average radius (surviving particles)')
axs[0].legend(fontsize=8);axs[1].set(xlabel='Time (h), literature-model calibration',ylabel='Particle count',yticks=[1,2,3,4])
for ax in axs:ax.grid(alpha=.2)
fig.savefig(OUT/'mean_radius_cubed_vs_time.png');plt.close(fig)

fig,axs=plt.subplots(1,2,figsize=(11,4.3),sharey=True,layout='constrained')
for ax,(key,run) in zip(axs,runs.items()):
    _,_,_,tr,_,title,color=run;use=tr[tr.time_star<=end+1e-8]
    for j in range(1,5):ax.plot(use.time_h,use[f'R{j}_nm'],label=f'P{j}')
    ax.set(xlabel='Time (h)',ylabel='Area-equivalent radius (nm)',title=title);ax.grid(alpha=.2);ax.legend()
fig.savefig(OUT/'particle_radii_comparison.png');plt.close(fig)

fig,axs=plt.subplots(1,3,figsize=(12,3.8),layout='constrained')
for key,run in runs.items():
    _,_,_,tr,post,title,color=run;p=post[post.time<=end+1e-8];u=tr[tr.time_star<=end+1e-8]
    axs[0].plot(p.time_h,(p.mean_c-.14)*1e10,color=color,label=title)
    axs[1].plot(p.time_h,p.total_energy*f0*1e-18,color=color)
    axs[2].plot(u.time_h,u.gamma_prime_area_fraction*100,color=color)
axs[0].set(ylabel=r'$(\bar c-0.14)\times 10^{10}$',title='Al conservation');axs[0].legend(fontsize=8)
axs[1].set(ylabel='Energy per out-of-plane length (J/m)',title='Total free energy')
axs[2].set(ylabel="Gamma-prime area fraction (%)",title='Phase fraction')
for ax in axs:ax.set_xlabel('Time (h)');ax.grid(alpha=.2)
fig.savefig(OUT/'conservation_energy_fraction.png');plt.close(fig)

_,nod,elem,_,_,_,_=runs['elastic']
fig,axs=plt.subplots(1,2,figsize=(9.5,4.2),layout='constrained')
lim=max(np.max(abs(elem[s][-1]))*f0/1e6 for s in ['stress_xx_scaled','stress_yy_scaled'])
for ax,s,title in zip(axs,['stress_xx_scaled','stress_yy_scaled'],[r'$\sigma_{xx}$',r'$\sigma_{yy}$']):
    im=ax.imshow(elem[s][-1]*f0/1e6,origin='lower',extent=[0,96,0,96],cmap='RdBu_r',vmin=-lim,vmax=lim)
    ax.contour(np.linspace(0,96,97),np.linspace(0,96,97),nod['eta'][-1],levels=[.5],colors='k',linewidths=.6)
    ax.set(xlabel='x (nm)',ylabel='y (nm)',title=title)
fig.colorbar(im,ax=axs,label='Stress (MPa)');fig.suptitle(f'Elastic case at {end*t0/3600:.3f} h | plane strain, zero mean in-plane stress')
fig.savefig(OUT/'elastic_stress_final.png');plt.close(fig)

fits=[]
for key,ff in metrics['fits'].items():
    for f in ff:
        f['lo_s']=f['lo_star']*t0;f['hi_s']=f['hi_star']*t0
        f['lo_h']=f['lo_s']/3600;f['hi_h']=f['hi_s']/3600
        f['K_nm3_per_model_time']=f['K_nm3_s']*t0
        fits.append({'case':key,**f})
pd.DataFrame(fits).to_csv(OUT/'coarsening_rate_fits.csv',index=False)
(OUT/'comparison_metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')
print(json.dumps({k:metrics[k] for k in ['status','no_elastic','elastic','fits']},indent=2))

if not args.partial and not args.no_gif:
    fig,axs=plt.subplots(1,2,figsize=(8,4),layout='constrained')
    images=[]
    for ax,(key,run) in zip(axs,runs.items()):
        images.append(ax.imshow(run[1]['eta'][0],origin='lower',extent=[0,96,0,96],vmin=0,vmax=1,cmap='Blues'))
        ax.set(title=run[5],xlabel='x (nm)',ylabel='y (nm)')
    suptitle=fig.suptitle('')
    targets=np.linspace(0,end,75)
    def update(i):
        for im,run in zip(images,runs.values()):im.set_data(run[1]['eta'][np.argmin(abs(run[0]-targets[i]))])
        suptitle.set_text(f'Ni-14 at.% Al, 1000 K | t = {targets[i]*t0/3600:.3f} h')
        return images+[suptitle]
    FuncAnimation(fig,update,frames=len(targets),interval=110).save(OUT/'elastic_comparison.gif',writer=PillowWriter(fps=9));plt.close(fig)

# Report is generated separately from these measured metrics.
