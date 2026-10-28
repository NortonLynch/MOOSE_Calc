#!/usr/bin/env bash
set -eo pipefail
cd "$(dirname "$0")/.."
bash scripts/run_moose.sh 4 -i examples/ni_al_coarsening/coarsening.i \
  --recover results/ni_al_coarsening/coarsening_cp/0900 \
  Executioner/num_steps=902 Executioner/dtmin=400 Executioner/dtmax=400 \
  Executioner/nl_max_its=100 \
  "Executioner/petsc_options=-snes_monitor -snes_converged_reason" \
  Outputs/file_base=results/ni_al_coarsening/benchmark \
  Outputs/fields/execute_on=NONE Outputs/restart/time_step_interval=99999999 Outputs/csv=false \
  > logs/ni_al_coarsening_benchmark.log 2>&1
