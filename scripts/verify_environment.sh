#!/usr/bin/env bash
set -eo pipefail
project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$project_root"
mkdir -p results/environment_smoke logs
bash scripts/run_moose.sh 1 -i examples/environment_smoke/ch_smoke.i \
  Outputs/file_base=results/environment_smoke/serial \
  2>&1 | tee logs/smoke-serial.log
bash scripts/run_moose.sh 2 -i examples/environment_smoke/ch_smoke.i \
  Outputs/file_base=results/environment_smoke/mpi2 \
  2>&1 | tee logs/smoke-mpi2.log
bash scripts/run_moose.sh 2 -i examples/environment_smoke/kks_smoke.i \
  2>&1 | tee logs/smoke-kks.log
source "$HOME/miniforge3-moose/etc/profile.d/conda.sh"
conda activate moose-post
python scripts/verify_smoke.py
