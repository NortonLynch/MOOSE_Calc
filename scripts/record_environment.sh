#!/usr/bin/env bash
set -eo pipefail
project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$HOME/miniforge3-moose/etc/profile.d/conda.sh"
conda activate moose
{
  date --iso-8601=seconds
  uname -a
  cat /etc/os-release
  free -h
  command -v moose-opt
  moose-opt --version
  mpiexec --version
  conda list
} > "$project_root/environment/installed-moose.txt"
conda activate moose-post
{
  python --version
  conda list
} > "$project_root/environment/installed-postprocess.txt"
conda list -n moose --explicit --md5 > "$project_root/environment/moose-linux-64.lock.txt"
conda list -n moose-post --explicit --md5 > "$project_root/environment/postprocess-linux-64.lock.txt"
