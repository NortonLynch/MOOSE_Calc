#!/usr/bin/env bash
set -eo pipefail
project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$project_root"
mkdir -p results/ni_al_coarsening logs
bash scripts/run_moose.sh "${1:-4}" -i examples/ni_al_coarsening/coarsening.i "${@:2}" > logs/ni_al_coarsening.log 2>&1
