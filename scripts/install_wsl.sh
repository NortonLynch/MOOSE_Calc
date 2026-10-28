#!/usr/bin/env bash
# Install the recorded linux-64 solver environment without rewriting its lock file.
set -eo pipefail
project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
conda_root="$HOME/miniforge3-moose"
release=26.7.2-0
installer="$project_root/downloads/Miniforge3-Linux-x86_64.sh"
installer_sha256=281b0ac7d550802efc81af633225a5e6116d29ae72f3ab4eae7168c3931a4c05
[[ "$(uname -s)" == Linux && "$(uname -m)" == x86_64 ]] || {
  echo 'This lock file requires Linux x86_64 (WSL2 / Ubuntu on Windows).' >&2
  exit 2
}
if [[ ! -x "$conda_root/bin/conda" ]]; then
  mkdir -p "$project_root/downloads"
  if [[ ! -f "$installer" ]]; then
    curl --fail --location --retry 3 \
      "https://github.com/conda-forge/miniforge/releases/download/$release/Miniforge3-Linux-x86_64.sh" \
      --output "$installer.part"
    mv "$installer.part" "$installer"
  fi
  printf '%s  %s\n' "$installer_sha256" "$installer" | sha256sum --check -
  bash "$installer" -b -p "$conda_root"
fi
source "$conda_root/etc/profile.d/conda.sh"
if [[ ! -f "$conda_root/envs/moose/conda-meta/history" ]]; then
  conda create --yes --name moose --file "$project_root/environment/moose-linux-64.lock.txt"
fi
conda activate moose
set -u
# Do not silently accept another solver version, or replace an existing environment.
python - <<'PY'
import json
import os
from pathlib import Path
expected = {'moose': '2026.10.14', 'moose-petsc': '3.25.4',
            'moose-libmesh': '2026.08.18_0772c3d', 'mpich': '5.0.1'}
installed = {}
for path in (Path(os.environ['CONDA_PREFIX']) / 'conda-meta').glob('*.json'):
    package = json.loads(path.read_text())
    installed[package['name']] = package['version']
for name, version in expected.items():
    if installed.get(name) != version:
        raise SystemExit(f'{name}: expected {version}, found {installed.get(name)}. '
                         'Existing environment was preserved; use a clean WSL user for exact deployment.')
print('Solver package versions match the recorded environment.')
PY
mkdir -p "$project_root/results/environment_smoke" "$project_root/logs"
cd "$project_root"
moose-opt --version
moose-opt -i examples/environment_smoke/ch_smoke.i --check-input
echo 'MOOSE deployment verified. Run: bash scripts/run_moose.sh 2 -i examples/environment_smoke/ch_smoke.i'
