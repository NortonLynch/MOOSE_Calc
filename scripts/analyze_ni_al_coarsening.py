"""Read MOOSE Exodus results, track particles, check conservation, and render.

The spatial coordinates are nm. t* is dimensionless reference-model time.
No trajectory or field is synthesized: all figures use the solver output.
"""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
import netCDF4
from scipy.ndimage import label
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

parser = argparse.ArgumentParser()
parser.add_argument('--stem', default='coarsening')
parser.add_argument('--partial', action='store_true')
parser.add_argument('--no-gif', action='store_true')
args = parser.parse_args()
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/ni_al_coarsening'
params = json.loads((OUT/'model_parameters.json').read_text())
with netCDF4.Dataset(OUT/f'{args.stem}.e') as ds:
    names = list(netCDF4.chartostring(ds['name_nod_var'][:]))
    x, y = np.asarray(ds['coordx'][:]), np.asarray(ds['coordy'][:])
    order = np.lexsort((x,y))
    nx, ny = len(np.unique(x)), len(np.unique(y))
    assert nx*ny == len(x)
    times = np.asarray(ds['time_whole'][:])
    fields = {}
    for name in ['eta','c','mu','gamma_prime_fraction','c_gamma','c_gamma_prime']:
        values = ds[f'vals_nod_var{names.index(name)+1}'][:]
        assert not np.ma.getmaskarray(values).any(), f'Incomplete field: {name}'
        fields[name] = np.asarray(values)[:,order].reshape((-1,ny,nx))
        assert np.isfinite(fields[name]).all(), name
    X,Y = x[order].reshape(ny,nx),y[order].reshape(ny,nx)

assert np.all(np.diff(times) > 0), 'Nonmonotonic or duplicate Exodus times'
raw_post = pd.read_csv(OUT/f'{args.stem}.csv')
# Recovery appends to CSV: discard superseded rows after each rewind to a checkpoint.
# Exodus itself overwrites those abandoned frames, and must remain strictly monotonic.
accepted = []
previous_raw_time = -np.inf
for row in raw_post.itertuples(index=False, name=None):
    if row[0] < previous_raw_time-1e-7:
        # Each recovery checkpoint in this run is also a saved Exodus frame
        # (checkpoint every 100 steps, Exodus every 5). The first new CSV row
        # is the NEXT step, so also remove abandoned rows between checkpoint
        # time and that new step, not just rows later than the new step.
        index = np.searchsorted(times, row[0]-1e-7)-1
        assert index >= 0
        checkpoint_time = times[index]
        while accepted and accepted[-1][0] > checkpoint_time+1e-7:
            accepted.pop()
    while accepted and accepted[-1][0] >= row[0]-1e-7:
        accepted.pop()
    accepted.append(row)
    previous_raw_time = row[0]
post = pd.DataFrame(accepted, columns=raw_post.columns)
assert np.all(np.diff(post.time) > 0)
eta,c = fields['eta'],fields['c']
h = eta**3*(6*eta**2-15*eta+10)
dx,dy = float(X[0,1]-X[0,0]),float(Y[1,0]-Y[0,0])
weights = np.ones((ny,nx))
weights[[0,-1],:] *= .5
weights[:,[0,-1]] *= .5
area = (X.max()-X.min())*(Y.max()-Y.min())
mean_c = (c*weights).sum(axis=(1,2))*dx*dy/area
mass_drift = float(np.max(np.abs(mean_c-0.14)))
assert mass_drift < 1e-7, mass_drift
assert float(np.max(np.abs(post.mean_c-.14))) < 1e-7
assert eta.min() > -0.03 and eta.max() < 1.03
assert c.min() > 0 and c.max() < 1

# The four well-separated particles stay inside their initial Voronoi quadrants.
# Integrating each quadrant by the Q1 trapezoidal weights partitions the domain.
mid_y,mid_x = (ny-1)//2,(nx-1)//2
slices = [(slice(0,mid_y+1),slice(0,mid_x+1)),
          (slice(0,mid_y+1),slice(mid_x,nx)),
          (slice(mid_y,ny),slice(0,mid_x+1)),
          (slice(mid_y,ny),slice(mid_x,nx))]
rows = []
for i,t in enumerate(times):
    row = {'time_star':float(t),'reference_time_s':float(t*params['scales']['time_s']),
           'mean_x_Al':float(mean_c[i]),'gamma_prime_area_fraction':float((h[i]*weights).sum()*dx*dy/area)}
    _, n_components = label(eta[i]>=.5)
    row['particle_count_eta_half'] = int(n_components)
    for j,(ys,xs) in enumerate(slices,1):
        hs = h[i,ys,xs]
        qw = np.ones_like(hs)
        qw[[0,-1],:]*=.5
        qw[:,[0,-1]]*=.5
        a = float((hs*qw).sum()*dx*dy)
        present = bool(np.max(eta[i,ys,xs]) >= .5)
        row[f'p{j}_area_nm2'] = a
        row[f'p{j}_equivalent_radius_nm'] = np.sqrt(max(0,a)/np.pi)
        row[f'p{j}_present'] = int(present)
        row[f'quadrant{j}_Al_integral_nm2'] = float((c[i,ys,xs]*qw).sum()*dx*dy)
    rows.append(row)
track = pd.DataFrame(rows)
prefix = '' if args.stem=='coarsening' else args.stem+'_'
post.to_csv(OUT/f'{prefix}diagnostics.csv',index=False,float_format='%.14g')
track.to_csv(OUT/f'{prefix}particle_history.csv',index=False,float_format='%.12g')
initial_r = [float(track[f'p{j}_equivalent_radius_nm'].iloc[0]) for j in range(1,5)]
final_r = [float(track[f'p{j}_equivalent_radius_nm'].iloc[-1]) for j in range(1,5)]
de = np.diff(post.total_energy)
energy_increase = float(max(0.,de.max(initial=0)))
report = {
    'status':'partial' if args.partial else 'complete',
    'exodus_file':f'{args.stem}.e','readback':'PASS','saved_frames':len(times),
    'nonlinear_time_steps':len(post)-1,'last_time_star':float(times[-1]),
    'superseded_restart_CSV_rows_removed':len(raw_post)-len(post),
    'maximum_mean_Al_deviation_from_0p14':mass_drift,
    'initial_particle_count':int(track.particle_count_eta_half.iloc[0]),
    'final_particle_count':int(track.particle_count_eta_half.iloc[-1]),
    'initial_equivalent_radii_nm':initial_r,'final_equivalent_radii_nm':final_r,
    'total_free_energy_initial_dimensionless':float(post.total_energy.iloc[0]),
    'total_free_energy_final_dimensionless':float(post.total_energy.iloc[-1]),
    'largest_positive_energy_step_dimensionless':energy_increase,
    'eta_range':[float(eta.min()),float(eta.max())],
    'x_Al_range':[float(c.min()),float(c.max())],
    'matrix_phase_x_Al_range':[float(fields['c_gamma'].min()),float(fields['c_gamma'].max())],
    'precipitate_phase_x_Al_range':[float(fields['c_gamma_prime'].min()),float(fields['c_gamma_prime'].max())],
    'radius_measure':'sqrt(integral of h(eta) in the original particle quadrant / pi)',
    'first_dissolution_time_star':{f'p{j}':(float(track.loc[track[f'p{j}_present']==0,'time_star'].iloc[0]) if (track[f'p{j}_present']==0).any() else None) for j in range(1,5)},
    'quadrant_Al_transfer_nm2':[float(track[f'quadrant{j}_Al_integral_nm2'].iloc[-1]-track[f'quadrant{j}_Al_integral_nm2'].iloc[0]) for j in range(1,5)],
    'scope':'Reference-parameter 2D demonstration, not an experimentally calibrated aging-time prediction',
}
if not args.partial:
    assert np.isclose(times[-1],post.time.iloc[-1]), 'Last field output missing'
    assert all(np.min(np.abs(post.time.to_numpy()-t))<1e-6 for t in times), 'Field time absent from accepted solver history'
    assert int(track.particle_count_eta_half.iloc[0]) == 4
    assert int(track.particle_count_eta_half.iloc[-1]) < 4, 'Continue until at least one particle dissolves'
    assert min(final_r)<min(initial_r)-1
    assert max(final_r)>max(initial_r)+1
    assert energy_increase < 1e-5*post.total_energy.iloc[0]
    report['coarsening_observed']='PASS: small-particle dissolution and large-particle growth'
(OUT/f'{prefix}validation.json').write_text(json.dumps(report,indent=2)+'\n')

plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
chosen = [0]+[int(np.argmin(np.abs(times-times[-1]*f))) for f in [.33,.66,1.]]
fig,axs = plt.subplots(2,4,figsize=(13,6.6),layout='constrained')
for col,k in enumerate(chosen):
    im1=axs[0,col].imshow(eta[k],origin='lower',extent=[0,X.max(),0,Y.max()],cmap='magma',vmin=0,vmax=1)
    im2=axs[1,col].imshow(c[k],origin='lower',extent=[0,X.max(),0,Y.max()],cmap='viridis',vmin=.124,vmax=.24)
    axs[0,col].set_title(f't* = {times[k]:,.0f}')
    for row in [0,1]:
        axs[row,col].set(xlabel='x (nm)',ylabel='y (nm)' if col==0 else '')
    if col==0:
        for j,(a,b,_) in enumerate(params['initial_particles_center_x_y_radius_nm'],1):
            axs[0,col].text(a,b,str(j),ha='center',va='center',color='black',weight='bold')
fig.colorbar(im1,ax=list(axs[0]),label=r'$\eta$: gamma-prime order parameter',shrink=.86)
fig.colorbar(im2,ax=list(axs[1]),label='Al mole fraction',shrink=.86)
fig.suptitle('Ni–14 at.% Al, 1000 K | MOOSE | reference-parameter 2D coarsening')
fig.savefig(OUT/f'{prefix}evolution.png',dpi=160)
plt.close(fig)
fig,axs=plt.subplots(1,3,figsize=(13,3.7),layout='constrained')
for j in range(1,5):
    axs[0].plot(track.time_star,track[f'p{j}_equivalent_radius_nm'],label=f'Particle {j}')
axs[0].set(xlabel='Time t* (reference model)',ylabel='Equivalent radius (nm)',title='Particle-size evolution')
axs[0].legend(fontsize=9)
axs[1].plot(post.time,post.mean_c-.14)
axs[1].set(xlabel='Time t*',ylabel='Mean Al fraction − 0.14',title='Solute conservation',ylim=(-max(1e-12,1.2*mass_drift),max(1e-12,1.2*mass_drift)))
axs[1].ticklabel_format(axis='y',style='sci',scilimits=(0,0))
axs[2].plot(post.time,post.total_energy)
axs[2].set(xlabel='Time t*',ylabel='Total reduced free energy (scaled)',title='Free-energy decay')
fig.savefig(OUT/f'{prefix}diagnostics.png',dpi=160)
plt.close(fig)

first_loss = np.flatnonzero(track.particle_count_eta_half.to_numpy()<4)
show_t = times[first_loss[0]]*0.8 if len(first_loss) else times[-1]*0.8
k = int(np.argmin(np.abs(times-show_t)))
mu = fields['mu'][k]
dmu_y,dmu_x = np.gradient(mu,dy,dx)
norm = np.hypot(dmu_x,dmu_y)
qx = np.divide(-dmu_x,norm,out=np.zeros_like(norm),where=norm>1e-14)
qy = np.divide(-dmu_y,norm,out=np.zeros_like(norm),where=norm>1e-14)
qx[eta[k]>.2]=np.nan
qy[eta[k]>.2]=np.nan
fig,ax=plt.subplots(figsize=(6.4,5.2),layout='constrained')
img=ax.imshow(mu*1000,origin='lower',extent=[0,X.max(),0,Y.max()],cmap='viridis')
ax.contour(X,Y,eta[k],levels=[.5],colors='white',linewidths=1.2)
sample=(slice(3,-3,7),slice(3,-3,7))
ax.quiver(X[sample],Y[sample],qx[sample],qy[sample],color='white',scale=30,width=.003)
ax.set(xlabel='x (nm)',ylabel='y (nm)',title=f'Diffusion-potential field at t*={times[k]:,.0f}\nArrows: direction of solute flux; white line: eta=0.5')
fig.colorbar(img,ax=ax,label='Reduced diffusion potential (J/mol)')
fig.savefig(OUT/f'{prefix}diffusion_evidence.png',dpi=160)
plt.close(fig)

if not args.no_gif:
    fig,ax=plt.subplots(figsize=(5.5,5),layout='constrained')
    im=ax.imshow(eta[0],origin='lower',extent=[0,X.max(),0,Y.max()],cmap='magma',vmin=0,vmax=1)
    ax.set(xlabel='x (nm)',ylabel='y (nm)')
    title=ax.set_title('')
    fig.colorbar(im,ax=ax,label=r'$\eta$')
    # Equal t* spacing in the preview; choose nearest actual saved frame, no field interpolation.
    indices=np.unique([np.argmin(abs(times-t)) for t in np.linspace(0,times[-1],100)])
    def update(k):
        im.set_data(eta[k])
        title.set_text(f'1000 K | t*={times[k]:,.0f}\nReference-parameter 2D coarsening')
        return im,title
    anim=FuncAnimation(fig,update,frames=indices,interval=100)
    anim.save(OUT/'coarsening_preview.gif',writer=PillowWriter(fps=10),dpi=110)
    plt.close(fig)
print(json.dumps(report,indent=2))
