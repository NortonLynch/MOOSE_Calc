#!/usr/bin/env bash
set -eo pipefail
project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$HOME/miniforge3-moose/etc/profile.d/conda.sh"
conda activate moose
set -u
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
ranks="${1:-1}"
shift || true
[[ "$ranks" =~ ^[1-9][0-9]*$ ]] || { echo 'MPI ranks must be a positive integer.' >&2; exit 2; }
cd "$project_root"
if [[ $# -eq 0 ]]; then
  set -- --help
fi
if [[ "$ranks" -eq 1 ]]; then
  exec moose-opt "$@"
else
  exec mpiexec -n "$ranks" moose-opt "$@"
fi
