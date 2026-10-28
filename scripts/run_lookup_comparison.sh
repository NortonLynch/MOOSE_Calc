#!/usr/bin/env bash
# Run in WSL/Linux. Original coarsening results are never overwritten.
set -eo pipefail
project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$project_root"
tag="${1:-rerun_$(date +%Y%m%d_%H%M%S)}"
[[ "$tag" =~ ^[A-Za-z0-9_-]+$ ]] || { echo 'Use letters, digits, _ or - for the run name.' >&2; exit 2; }
out="results/ni_al_lookup_$tag"
[[ ! -e "$out" ]] || { echo "Output already exists: $out" >&2; exit 2; }
source "$HOME/miniforge3-moose/etc/profile.d/conda.sh"
conda activate moose-post
python scripts/prepare_lookup_comparison.py
mkdir -p "$out" logs
cp results/ni_al_lookup/thermodynamic_validation.json "$out/"
for name in quadratic_control lookup; do
  bash scripts/run_moose.sh "${MOOSE_RANKS:-4}" -i "examples/ni_al_lookup/$name.i" \
    "Outputs/file_base=$out/$name" > "$out/$name.log" 2>&1
done
python scripts/prepare_lookup_tail.py --output "$out"
bash scripts/run_moose.sh "${MOOSE_RANKS:-4}" -i examples/ni_al_lookup/quadratic_matched_tail.i \
  > "$out/quadratic_matched_tail.log" 2>&1
python scripts/assemble_lookup_control.py --output "$out"
python scripts/analyze_lookup_comparison.py --output "$out"
echo "Report: $out/REPORT.html"
