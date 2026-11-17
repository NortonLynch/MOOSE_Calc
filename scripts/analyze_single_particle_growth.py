"""Track the same final-largest particle in existing MOOSE fields, not a population mean."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import netCDF4
import contourpy
from scipy.ndimage import label
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path as PolygonPath

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'results/ni_al_elastic'
OUT = BASE / 'single_particle'
OUT.mkdir(exist_ok=True)
cal = json.loads((BASE / 'time_calibration.json').read_text())
old = json.loads((BASE / 'comparison_metrics.json').read_text())
assert old['status'] == 'COMPLETE'
t0 = cal['time_scale_s_per_model_unit']
factor = np.pi**1.5
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False,
                     'axes.spines.right': False, 'savefig.dpi': 180})
cases = [('no_elastic', ROOT / 'results/ni_al_lookup/lookup.e', 'Without elasticity', '#277da8'),
         ('elastic', BASE / 'elastic.e', 'With elasticity', '#cf5b36')]
windows = [(0, 60000, '0-60000'), (5000, 60000, '5000-60000'),
           (5000, 30000, '5000-30000'), (5000, 15000, '5000-15000'),
           (15000, 30000, '15000-30000'), (30000, 45000, '30000-45000'),
           (45000, 60000, '45000-60000')]


def integrate(q, dx, dy):
    w = np.ones(q.shape[-2:])
    w[[0, -1], :] *= .5
    w[:, [0, -1]] *= .5
    return (q * w).sum(axis=(-2, -1)) * dx * dy


def regress(tr, lo, hi, name, metric):
    use = tr[(tr.time_star >= lo) & (tr.time_star <= hi)]
    if len(use) < 4:
        return None
    x = use.time_s.to_numpy()
    y = use[metric].to_numpy()
    w = np.r_[np.diff(x)[0], x[2:] - x[:-2], np.diff(x)[-1]] / 2
    xm = np.average(x, weights=w)
    ym = np.average(y, weights=w)
    slope = np.dot(w * (x - xm), y - ym) / np.dot(w * (x - xm), x - xm)
    intercept = ym - slope * xm
    residual = y - (intercept + slope*x)
    sse = np.dot(w, residual**2)
    sst = np.dot(w, (y - ym)**2)
    # Curvature diagnostic only; deterministic time series has no independent noise model.
    z = (x - xm) / np.ptp(x)
    design = np.column_stack([np.ones_like(z), z, z*z])
    q = np.linalg.lstsq(design*np.sqrt(w)[:, None], y*np.sqrt(w), rcond=None)[0]
    qres = y - design @ q
    return dict(window=name, metric=metric, requested_lo_star=lo, requested_hi_star=hi,
                lo_star=float(use.time_star.iloc[0]), hi_star=float(use.time_star.iloc[-1]),
                lo_s=float(x[0]), hi_s=float(x[-1]), lo_h=float(x[0]/3600), hi_h=float(x[-1]/3600),
                samples=len(use), K_nm3_s=float(slope), K_m3_s=float(slope*1e-27),
                intercept_nm3=float(intercept), R2=float(1-sse/sst),
                RMSE_nm3=float(np.sqrt(sse/w.sum())),
                RMSE_percent_of_observed_range=float(100*np.sqrt(sse/w.sum())/np.ptp(y)),
                max_abs_residual_nm3=float(abs(residual).max()),
                quadratic_R2=float(1-np.dot(w, qres*qres)/sst),
                quadratic_RMSE_nm3=float(np.sqrt(np.dot(w, qres*qres)/w.sum())),
                population_counts=sorted(map(int, use.population_count.unique())),
                fitted_particle='P4 fixed identity; not the maximum reselected at each time')


runs = {}
fits = []
metrics = {'status': 'COMPLETE', 'tracked_particle': 'P4', 'initial_center_nm': [72, 72],
           'area_definition': 'A_h = trapezoidal integral of h(eta) in the fixed P4 quadrant, h=eta^3(6eta^2-15eta+10); same convention as the previous analysis',
           'contour_crosscheck': 'Area of the closed eta=0.5 contour containing the fixed P4 center; polygon shoelace area',
           'radius_definition': 'r_eq=sqrt(A/pi); r_eq^3=A^(3/2)/pi^(3/2)',
           'pi_power_3_2': float(factor), 'time_scale_s': t0, 'time_calibration': cal,
           'fitting': 'time-weighted least squares with trapezoidal time weights and free intercept; R2 is descriptive, not a test of universal LSW kinetics',
           'initial_relaxation_cutoff_star': 5000, 'cases': {}}
for key, source, title, color in cases:
    with netCDF4.Dataset(source) as d:
        names = list(netCDF4.chartostring(d['name_nod_var'][:]))
        x, y = np.asarray(d['coordx'][:]), np.asarray(d['coordy'][:])
        ux, uy = np.unique(x), np.unique(y)
        order = np.lexsort((x, y))
        eta = np.asarray(d[f'vals_nod_var{names.index("eta")+1}'][:])[:, order].reshape(-1, len(uy), len(ux))
        times = np.asarray(d['time_whole'][:])
    assert np.isfinite(eta).all() and np.all(np.diff(times) > 0)
    assert times[0] == 0 and abs(times[-1]-60000) < 1e-7
    assert np.allclose(np.diff(ux), 1) and np.allclose(np.diff(uy), 1)
    midx, midy = len(ux)//2, len(uy)//2
    quadrants = [(slice(0, midy+1), slice(0, midx+1)), (slice(0, midy+1), slice(midx, None)),
                 (slice(midy, None), slice(0, midx+1)), (slice(midy, None), slice(midx, None))]
    h = eta**3 * (6*eta**2-15*eta+10)
    all_areas = np.stack([integrate(h[:, ys, xs], 1, 1) for ys, xs in quadrants], axis=1)
    final_largest_id = int(np.argmax(all_areas[-1]))+1
    assert final_largest_id == 4, 'Track identity must be reconsidered if final-largest particle changes.'
    # Identity remains P4 even if it were not the instantaneous maximum in an intermediate frame.
    area = all_areas[:, 3]
    radius = np.sqrt(area/np.pi)
    contour_areas = []
    for z in eta:
        paths = contourpy.contour_generator(x=ux, y=uy, z=z).lines(.5)
        polygons = [p for p in paths if PolygonPath(p).contains_point((72, 72))]
        assert len(polygons) == 1
        p = polygons[0]
        assert np.allclose(p[0], p[-1]), 'P4 contour must be closed.'
        assert (p[:, 0] > 48).all() and (p[:, 1] > 48).all()
        assert (p[:, 0] < 96).all() and (p[:, 1] < 96).all()
        contour_areas.append(.5*abs(np.dot(p[:-1, 0], p[1:, 1])-np.dot(p[1:, 0], p[:-1, 1])))
    contour_areas = np.array(contour_areas)
    previous = pd.read_csv(BASE/f'{key}_particle_history.csv')
    assert np.allclose(previous.time_star, times, rtol=0, atol=1e-8)
    error = float(np.max(abs(previous.R4_nm-radius)))
    assert error < 1e-11, error
    counts = np.array([label(z >= .5)[1] for z in eta])
    assert np.array_equal(counts, previous['count'])
    assert (previous.present4 == 1).all()
    tr = pd.DataFrame({'time_star': times, 'time_s': times*t0, 'time_h': times*t0/3600,
                       'particle_id': 4, 'area_nm2': area, 'area_3_2_nm3': area**1.5,
                       'radius_nm': radius, 'radius_cubed_nm3': radius**3,
                       'contour_area_nm2': contour_areas, 'contour_area_3_2_nm3': contour_areas**1.5,
                       'contour_radius_nm': np.sqrt(contour_areas/np.pi),
                       'contour_radius_cubed_nm3': (contour_areas/np.pi)**1.5,
                       'population_count': counts})
    assert np.allclose(tr.area_3_2_nm3/factor, tr.radius_cubed_nm3, rtol=1e-14)
    tr.to_csv(OUT/f'{key}_P4_history.csv', index=False)
    local_windows = list(windows)
    boundaries = np.r_[0, np.flatnonzero(np.diff(counts) != 0)+1, len(times)]
    for a, b in zip(boundaries[:-1], boundaries[1:]):
        local_windows.append((max(5000, times[a]), times[b-1], f'population_N{counts[a]}_{a}_{b-1}'))
    for lo, hi, name in local_windows:
        for metric in ['radius_cubed_nm3', 'area_3_2_nm3', 'contour_radius_cubed_nm3']:
            f = regress(tr, lo, hi, name, metric)
            if f:
                fits.append({'case': key, **f})
    metrics['cases'][key] = {'source': str(source.relative_to(ROOT)), 'saved_frames': len(times),
        'selected_final_largest_id': final_largest_id, 'P4_is_largest_all_saved_frames': bool(np.all(np.argmax(all_areas, axis=1) == 3)),
        'max_radius_difference_from_previous_extraction_nm': error,
        'initial_radius_nm': float(radius[0]), 'final_radius_nm': float(radius[-1]),
        'initial_area_nm2': float(area[0]), 'final_area_nm2': float(area[-1]),
        'initial_r3_nm3': float(radius[0]**3), 'final_r3_nm3': float(radius[-1]**3),
        'initial_A32_nm3': float(area[0]**1.5), 'final_A32_nm3': float(area[-1]**1.5),
        'contour_max_radius_relative_difference_percent': float(np.max(abs(np.sqrt(contour_areas/np.pi)/radius-1))*100),
        'contour_all_frames_valid': True,
        'dissolution_events_star': old[key]['dissolution_brackets_star']}
    runs[key] = (tr, eta, title, color)
fit_table = pd.DataFrame(fits)
fit_table.to_csv(OUT/'P4_growth_fits.csv', index=False)


def selected_fit(key, window, metric='radius_cubed_nm3'):
    return next(f for f in fits if f['case'] == key and f['window'] == window and f['metric'] == metric)


for key in runs:
    for window in ['5000-30000', '5000-60000']:
        rf = selected_fit(key, window)
        af = selected_fit(key, window, 'area_3_2_nm3')
        assert np.isclose(af['K_nm3_s']/rf['K_nm3_s'], factor, rtol=1e-12)
        assert np.isclose(af['R2'], rf['R2'], rtol=1e-12)
    metrics['cases'][key]['primary_fit'] = selected_fit(key, '5000-30000')
    metrics['cases'][key]['whole_fit'] = selected_fit(key, '5000-60000')
    metrics['cases'][key]['contour_primary_fit'] = selected_fit(key, '5000-30000', 'contour_radius_cubed_nm3')
metrics['primary_K_elastic_to_no_elastic_ratio'] = selected_fit('elastic', '5000-30000')['K_nm3_s']/selected_fit('no_elastic', '5000-30000')['K_nm3_s']
metrics['validation'] = 'PASS: final particle identity, closed isolated contour in every frame, original-radius agreement, area-radius identity, slope scaling identity, time conversion, finite fields'
(OUT/'metrics.json').write_text(json.dumps(metrics, indent=2)+'\n')


def save(fig, name):
    fig.savefig(OUT/f'{name}.png')
    fig.savefig(OUT/f'{name}.svg')
    plt.close(fig)


for metric, name, ylabel in [('radius_cubed_nm3', 'P4_radius_cubed_vs_time', r'$r_{P4}^3$ (nm$^3$)'),
                             ('area_3_2_nm3', 'P4_area_3_2_vs_time', r'$A_{P4}^{3/2}$ (nm$^3$)')]:
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.8), layout='constrained', sharex=True)
    for ax, (key, (tr, eta, title, color)) in zip(axs, runs.items()):
        ax.plot(tr.time_h, tr[metric], color=color, lw=2, label='Fixed particle P4')
        for window, style, fitcolor, labeltext in [('5000-30000', '--', '#1b6b41', 'Early fit'), ('5000-60000', ':', '#333333', 'Whole-window fit')]:
            f = selected_fit(key, window, metric)
            xx = np.array([f['lo_s'], f['hi_s']])
            ax.plot(xx/3600, f['intercept_nm3']+f['K_nm3_s']*xx, style, color=fitcolor, lw=1.5,
                    label=f'{labeltext}: K={f["K_nm3_s"]:.4g} nm³/s; R²={f["R2"]:.4f}')
        for j, bracket in old[key]['dissolution_brackets_star'].items():
            if bracket:
                ev = np.mean(bracket)*t0/3600
                ax.axvline(ev, color='#777777', alpha=.35, lw=.8)
                ax.text(ev, .52, f' P{j} dissolves', transform=ax.get_xaxis_transform(), rotation=90, va='center', fontsize=8, color='#666666')
        ax.set(title=title, xlabel='Physical time (h), literature-model calibration', ylabel=ylabel)
        ax.grid(alpha=.2); ax.legend(fontsize=8, loc='upper left')
    fig.suptitle('Ni-14 at.% Al, 1000 K | same P4 tracked at every frame; no population averaging')
    save(fig, name)

fig, axs = plt.subplots(2, 2, figsize=(12, 6.8), layout='constrained')
for col, (key, (tr, eta, title, color)) in enumerate(runs.items()):
    for window, fitcolor in [('5000-30000', '#1b6b41'), ('5000-60000', '#333333')]:
        f = selected_fit(key, window)
        u = tr[(tr.time_star >= f['lo_star']) & (tr.time_star <= f['hi_star'])]
        residual = u.radius_cubed_nm3-(f['intercept_nm3']+f['K_nm3_s']*u.time_s)
        axs[0, col].plot(u.time_h, residual, color=fitcolor, label=f't*={window}')
    axs[0, col].axhline(0, color='gray', lw=.7)
    axs[0, col].set(title=title, ylabel=r'Linear-fit residual in $r^3$ (nm$^3$)')
    axs[0, col].legend(fontsize=8)
    for lo, hi, window in windows[3:]:
        f = selected_fit(key, window)
        axs[1, col].hlines(f['K_nm3_s'], f['lo_h'], f['hi_h'], color=color, lw=4)
        axs[1, col].plot((f['lo_h']+f['hi_h'])/2, f['K_nm3_s'], 'o', color=color)
    axs[1, col].set(ylabel=r'Window slope $K_{P4}$ (nm$^3$/s)')
    for ax in axs[:, col]:
        ax.set_xlabel('Physical time (h)'); ax.grid(alpha=.2)
save(fig, 'P4_linearity_diagnostics')

fig, axs = plt.subplots(1, 3, figsize=(11.4, 4), layout='constrained')
for ax, (key, idx, title) in zip(axs, [('no_elastic', 0, 'Identical initial state'), ('no_elastic', -1, 'Without elasticity: final'), ('elastic', -1, 'With elasticity: final')]):
    tr, eta, _, _ = runs[key]
    ax.imshow(eta[idx], origin='lower', extent=[0, 96, 0, 96], vmin=0, vmax=1, cmap='Blues')
    ax.contour(ux, uy, eta[idx], levels=[.5], colors='k', linewidths=.65)
    ax.annotate('P4: tracked', xy=(72, 72), xytext=(8, 91), color='#b13e25', fontsize=11,
                arrowprops={'arrowstyle': '->', 'color': '#b13e25', 'lw': 1.5})
    ax.set(title=title, xlabel='x (nm)', ylabel='y (nm)')
save(fig, 'P4_particle_identity')

fig, axs = plt.subplots(1, 2, figsize=(11, 4.2), layout='constrained')
for ax, (key, (tr, eta, title, color)) in zip(axs, runs.items()):
    ax.plot(tr.time_h, tr.radius_cubed_nm3, color=color, label='Diffuse area integral')
    ax.plot(tr.time_h, tr.contour_radius_cubed_nm3, '--', color='#333333', label=r'$\eta=0.5$ contour area')
    ax.set(title=title, xlabel='Physical time (h)', ylabel=r'P4 equivalent-radius cube (nm$^3$)')
    ax.grid(alpha=.2); ax.legend(fontsize=9)
save(fig, 'P4_area_definition_check')
print(json.dumps({'validation': metrics['validation'], 'cases': metrics['cases'],
                  'primary_K_elastic_to_no_elastic_ratio': metrics['primary_K_elastic_to_no_elastic_ratio']}, indent=2))
