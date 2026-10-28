"""Publish measured diameters at accepted solver steps, using existing result CSVs."""
from pathlib import Path
import argparse
import hashlib
import json

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', type=Path, default=Path(__file__).resolve().parents[1]
                        / 'results/ni_al_coarsening')
    args = parser.parse_args()
    out = args.results
    history = pd.read_csv(out / 'particle_history.csv')
    diagnostics = pd.read_csv(out / 'diagnostics.csv')
    times = diagnostics.time.to_numpy()
    saved = history.time_star.to_numpy()
    if not (np.isfinite(times).all() and np.isfinite(saved).all()
            and np.all(np.diff(times) > 0) and np.all(np.diff(saved) > 0)):
        raise ValueError('Both input time axes must be finite and strictly increasing.')
    # diagnostics has the initial state, then one row per accepted step after
    # removal of abandoned restart branches. Its row index is the step number.
    right = np.clip(np.searchsorted(times, saved), 0, len(times) - 1)
    left = np.maximum(right - 1, 0)
    steps = np.where(abs(times[left] - saved) <= abs(times[right] - saved), left, right)
    error = abs(times[steps] - saved)
    if np.max(error) > 1e-6 or not np.all(np.diff(steps) > 0):
        raise ValueError('Saved frames do not uniquely match accepted solver steps.')
    table = pd.DataFrame({'saved_frame': np.arange(len(history)),
                          'accepted_step': steps, 'time_star': saved,
                          'particle_count': history.particle_count_eta_half})
    for j in range(1, 5):
        diameter = 2 * history[f'p{j}_equivalent_radius_nm'].to_numpy()
        if not (np.isfinite(diameter).all() and np.all(diameter >= 0)):
            raise ValueError(f'Invalid diameter for particle {j}')
        table[f'p{j}_diameter_nm'] = diameter
        table[f'p{j}_present'] = history[f'p{j}_present']
    table.to_csv(out / 'particle_diameter_vs_step.csv', index=False, float_format='%.12g')

    plt.rcParams.update({'font.size': 11, 'axes.spines.top': False,
                         'axes.spines.right': False})
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), layout='constrained')
    for j, color in enumerate(['#D55E00', '#0072B2', '#009E73', '#CC79A7'], 1):
        for ax, x in zip(axes, [steps, saved]):
            ax.plot(x, table[f'p{j}_diameter_nm'], color=color, lw=2, label=f'Particle {j}')
    axes[0].set(xlabel='Accepted solver time step', ylabel='Equivalent diameter (nm)',
                title='Diameter vs. solver step')
    axes[1].set(xlabel=r'Model time $t^*$', ylabel='Equivalent diameter (nm)',
                title='Diameter vs. model time')
    for ax in axes:
        ax.grid(alpha=.2)
        ax.set_ylim(bottom=0)
        ax.legend(fontsize=9, loc='lower left')
    fig.suptitle('Ni-14 at.% Al, 1000 K | 2D reference-parameter coarsening')
    fig.savefig(out / 'particle_diameter_vs_step.png', dpi=180)
    plt.close(fig)
    report = {
        'diameter_definition': '2*sqrt(integral h(eta) dA in original quadrant / pi)',
        'accepted_step_definition': 'row index in diagnostics.csv; initial state = 0',
        'saved_frames': len(table), 'final_accepted_step': int(steps[-1]),
        'maximum_time_match_error': float(error.max()),
        'initial_diameters_nm': [float(table[f'p{j}_diameter_nm'].iloc[0]) for j in range(1, 5)],
        'final_diameters_nm': [float(table[f'p{j}_diameter_nm'].iloc[-1]) for j in range(1, 5)],
        'first_absent_saved_step': {f'p{j}': (int(table.loc[table[f'p{j}_present'] == 0,
                                                        'accepted_step'].iloc[0])
                                           if (table[f'p{j}_present'] == 0).any() else None)
                                    for j in range(1, 5)},
        'input_sha256': {name: hashlib.sha256((out / name).read_bytes()).hexdigest()
                         for name in ['particle_history.csv', 'diagnostics.csv']},
        'note': 'Diameters exist at saved frames only, not at every accepted step. '
                'Model time is not experimentally calibrated aging time.'
    }
    (out / 'diameter_summary.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
