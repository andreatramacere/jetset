#!/bin/bash
set -ex

echo ">>>>>>>>>>>>>>>>>>> $PWD"
echo "Using Python: $PYTHON"
$PYTHON --version

# Install into the conda-build prefix
$PYTHON -m pip install . --no-deps -vv