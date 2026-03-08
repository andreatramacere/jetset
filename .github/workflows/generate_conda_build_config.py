import os
import yaml

# Get the Python version from the matrix (passed as an env var)
python_version = os.environ.get("MATRIX_PYTHON_VERSION")
if not python_version:
    raise RuntimeError("MATRIX_PYTHON_VERSION not set! Pass it from the workflow matrix.")

# Strip quotes and normalize
python_version = python_version.strip()

# Optional: add more versions if you want to build multiple at once
# Here we just write the single version from the matrix
config = {
    "python": [python_version],
    "numpy": ["2.4.2"],  # pin to numpy >= 2
}

# Write the conda_build_config.yaml file
output_path = ".github/conda-pipeline/github/conda_build_config.yaml"
os.makedirs(os.path.dirname(output_path), exist_ok=True)

with open(output_path, "w") as f:
    yaml.dump(config, f, sort_keys=False)

print(f"[INFO] Wrote {output_path} with Python={python_version} and NumPy >=2")
