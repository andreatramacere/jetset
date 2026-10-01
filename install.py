#!/usr/bin/env python3
"""Install JetSeT from source using the current Python environment."""

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


SCRIPT_DIR = Path(__file__).resolve().parent
PIN_START = "# >>> jetset constraints >>>"
PIN_END = "# <<< jetset constraints <<<"


def sync_conda_pins(prefix):
    """Replace JetSeT's marked constraints, preserving other Conda pins."""
    specs = []
    for line in (SCRIPT_DIR / "requirements.txt").read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        spec = "".join(line.split())
        # Match the original installer's filtering of pip-only entries.
        if spec.startswith("-") or "://" in spec or "@" in spec:
            continue
        specs.append(spec)

    if not specs:
        print("No conda-compatible constraints found in requirements.txt; skipping pin sync.")
        return

    pin_file = prefix / "conda-meta" / "pinned"
    retained = []
    skip = False
    if pin_file.is_file():
        for line in pin_file.read_text().splitlines():
            if line == PIN_START:
                skip = True
            elif line.startswith(PIN_END):
                skip = False
            elif not skip:
                retained.append(line)

    content = "\n".join(retained + [PIN_START] + specs + [PIN_END]) + "\n"
    # Write beside the destination so replacement is atomic on its filesystem.
    temporary_path = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", dir=pin_file.parent,
                                         delete=False) as temporary:
            temporary_path = Path(temporary.name)
            temporary.write(content)
        temporary_path.replace(pin_file)
    finally:
        if temporary_path is not None and temporary_path.exists():
            temporary_path.unlink()
    print("Synced JetSeT constraints to: {}".format(pin_file))


def run(command, cwd=SCRIPT_DIR):
    """Run a command and stop installation if it fails."""
    subprocess.run(command, cwd=cwd, check=True)


def install(skip_dependencies=False):
    """Install dependencies and JetSeT, then check the installed package."""
    pip = [sys.executable, "-m", "pip"]
    if skip_dependencies:
        print("Skipping dependency install (-skip-dep).", flush=True)
    else:
        prefix_value = os.environ.get("CONDA_PREFIX")
        prefix = Path(prefix_value) if prefix_value else None
        active_conda = prefix is not None and (prefix / "conda-meta").is_dir()
        conda = next((path for name in ("micromamba", "mamba", "conda")
                      for path in [shutil.which(name)] if path), None)
        if active_conda and conda:
            print("Detected active conda env: {}".format(prefix), flush=True)
            sync_conda_pins(prefix)
            print("Using {} for requirements.".format(conda), flush=True)
            run([conda, "install", "--yes", "-c", "astropy", "-c", "conda-forge",
                 "--file", "requirements.txt"])
        else:
            print("Using pip for requirements (no active Conda environment or frontend).",
                  flush=True)
            run(pip + ["install", "-r", "requirements.txt"])

    arguments = ["install", "--no-deps"]
    if skip_dependencies:
        # Use existing build tools instead of downloading an isolated build environment.
        arguments.append("--no-build-isolation")
    run(pip + arguments + ["."])

    # Check from a subdirectory so the source tree does not shadow the installed package.
    check_directory = SCRIPT_DIR / "tmp"
    check_directory.mkdir(exist_ok=True)
    run([sys.executable, "-c", """
import jetset
from jetset.test_data_helper import test_SEDs
from jetset.ebl_data import *
from jetset.Spectral_Templates_Repo import *
from jetset.jetkernel import mathkernel
from jetset.jetkernel import jetkernel
print(jetkernel.__file__)
print(jetset.__file__)
print(jetset.__version__)
"""], cwd=check_directory)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-skip-dep", "--skip-dep", action="store_true",
                        help="Skip dependencies and use already-installed build tools.")
    args = parser.parse_args()
    try:
        install(skip_dependencies=args.skip_dep)
    except (OSError, subprocess.CalledProcessError) as error:
        print("Error while installing JetSeT: {}".format(error), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
