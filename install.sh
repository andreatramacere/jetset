#!/usr/bin/env bash
set -euo pipefail

SKIP_DEP=0
IN_CONDA_ENV=0
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

while [[ $# -gt 0 ]]; do
  case "$1" in
    -skip-dep)
      SKIP_DEP=1
      shift
      ;;
    -h|--help)
      echo "Usage: $0 [-skip-dep]"
      echo "  -skip-dep    Skip dependency installation."
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      echo "Usage: $0 [-skip-dep]"
      exit 2
      ;;
  esac
done

# Work from repository root (where requirements.txt lives)
cd "$SCRIPT_DIR"

# Detect active conda-style env (conda/mamba/micromamba)
if [[ -n "${CONDA_PREFIX:-}" && -d "${CONDA_PREFIX}/conda-meta" ]]; then
  IN_CONDA_ENV=1
fi

# Pick conda frontend if available
if command -v micromamba >/dev/null 2>&1; then
  CONDA_CMD="micromamba"
elif command -v mamba >/dev/null 2>&1; then
  CONDA_CMD="mamba"
elif command -v conda >/dev/null 2>&1; then
  CONDA_CMD="conda"
else
  CONDA_CMD=""
fi

sync_conda_pins() {
  local pin_file tmp_specs tmp_pin

  pin_file="${CONDA_PREFIX}/conda-meta/pinned"
  tmp_specs="$(mktemp)"
  tmp_pin="$(mktemp)"

  # Keep only conda-compatible specs from requirements.txt and normalize spaces.
  awk '
    /^[[:space:]]*#/ || /^[[:space:]]*$/ { next }
    {
      gsub(/[[:space:]]+/, "", $0)
      if ($0 ~ /^-/) next
      if ($0 ~ /:\/\//) next
      if ($0 ~ /@/) next
      print $0
    }
  ' requirements.txt > "$tmp_specs"

  if [[ ! -s "$tmp_specs" ]]; then
    echo "No conda-compatible constraints found in requirements.txt; skipping pin sync."
    rm -f "$tmp_specs" "$tmp_pin"
    return
  fi

  {
    if [[ -f "$pin_file" ]]; then
      awk '
        BEGIN { skip=0 }
        /^# >>> jetset constraints >>>$/ { skip=1; next }
        /^# <<< jetset constraints <<</ { skip=0; next }
        skip==0 { print }
      ' "$pin_file"
    fi
    echo "# >>> jetset constraints >>>"
    cat "$tmp_specs"
    echo "# <<< jetset constraints <<<"
  } > "$tmp_pin"

  mv "$tmp_pin" "$pin_file"
  rm -f "$tmp_specs"
  echo "Synced JetSeT constraints to: $pin_file"
}

if [[ "$SKIP_DEP" -eq 0 ]]; then
  if [[ "$IN_CONDA_ENV" -eq 1 && -n "$CONDA_CMD" ]]; then
    echo "Detected active conda env: ${CONDA_PREFIX}"
    sync_conda_pins
    echo "Using $CONDA_CMD for requirements."
    "$CONDA_CMD" install --yes -c astropy -c conda-forge --file requirements.txt
  elif [[ -z "$CONDA_CMD" ]]; then
    echo "No conda/mamba found. Using pip for requirements."
    python -m pip install -r requirements.txt
  else
    echo "conda/mamba found but no active conda env. Using pip for requirements."
    python -m pip install -r requirements.txt
  fi
else
  echo "Skipping dependency install (-skip-dep)."
fi

# Install package from source with pip (deps already handled above)
PIP_INSTALL_ARGS=(--no-deps)
if [[ "$SKIP_DEP" -eq 1 ]]; then
  # In skip-dep mode, rely on already-installed build tools and avoid networked
  # build-isolation environments created from pyproject.toml requirements.
  PIP_INSTALL_ARGS+=(--no-build-isolation)
fi
python -m pip install "${PIP_INSTALL_ARGS[@]}" .
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
print(jetset.__file__)
print(jetset.__version__)

PY
cd ..
code=$?
if [ $code -ne 0 ]; then
  printf "\n\33[31mError while installing jetset.\33[0m\n"
  exit -63
fi
