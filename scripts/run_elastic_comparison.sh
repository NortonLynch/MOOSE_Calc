#!/usr/bin/env bash
set -eo pipefail
cd "$(dirname "$0")/.."
if [[ -f results/ni_al_elastic/elastic.e ]]; then
  echo 'results/ni_al_elastic/elastic.e already exists. Use a separate project copy for a fresh run; completed results are preserved.' >&2
  exit 2
fi
[[ -f results/ni_al_lookup/lookup.e ]] || { echo 'Restore the completed lookup reference fields first (scripts/download_lookup_results.py).' >&2; exit 2; }
mkdir -p logs results/ni_al_elastic
source "$HOME/miniforge3-moose/etc/profile.d/conda.sh"
conda activate moose-post
set -u
python scripts/prepare_elastic_comparison.py
bash scripts/run_moose.sh 8 -i examples/ni_al_elastic/elastic_initial.i > logs/elastic_initial.log 2>&1
python scripts/initialize_elastic_transient.py
bash scripts/run_moose.sh 4 -i examples/ni_al_elastic/zero_misfit.i Executioner/num_steps=20 > logs/zero_misfit.log 2>&1
bash scripts/run_moose.sh 8 -i examples/ni_al_elastic/elastic_relaxed.i > logs/elastic.log 2>&1
python scripts/calibrate_elastic_time.py
python scripts/check_elastic_controls.py
python scripts/analyze_elastic_comparison.py
python scripts/export_elastic_physical.py
python scripts/write_elastic_report.py
python scripts/analyze_single_particle_growth.py
python scripts/write_single_particle_report.py
echo 'Completed: results/ni_al_elastic/REPORT.html'
echo 'Fixed-particle kinetics: results/ni_al_elastic/single_particle/REPORT.html'
