#!/usr/bin/env bash
set -euo pipefail

SKIP_DEP=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    -skip-dep)
      SKIP_DEP=1
      shift
      ;;
    -h|--help)
      echo "Usage: $0 [-skip-dep]"
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      echo "Usage: $0 [-skip-dep]"
      exit 2
      ;;
  esac
done

# Pick conda frontend if available
if command -v mamba >/dev/null 2>&1; then
  CONDA_CMD="mamba"
elif command -v conda >/dev/null 2>&1; then
  CONDA_CMD="conda"
else
  CONDA_CMD=""
fi

if [[ "$SKIP_DEP" -eq 0 ]]; then
  if [[ -z "$CONDA_CMD" ]]; then
    echo "No conda/mamba found. Using pip for requirements."
    pip install -r requirements.txt
  else
    echo "Using $CONDA_CMD for requirements."
    "$CONDA_CMD" install -c astropy -c conda-forge --file requirements.txt
  fi
else
  echo "Skipping dependency install (-skip-dep)."
fi

# Your package install + import smoke test
pip install --no-build-isolation --no-deps .
mkdir -p tmp
cd tmp
python - <<'PY'
import jetset
from jetset.test_data_helper import test_SEDs
from jetset.ebl_data import *
from jetset.Spectral_Templates_Repo import *
from jetset.jetkernel import mathkernel
from jetset.jetkernel import jetkernel
print(jetkernel.__file__)
import jetset
print(jetset.__file__)
print(jetset.__version__)

PY
cd ..
code=$?
if [ $code -ne 0 ]; then
  printf "\n\33[31mError while installing jetset.\33[0m\n"
  exit -63
fi
