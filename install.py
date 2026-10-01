#!/usr/bin/env python3
"""Install JetSeT from source using the current Python environment."""

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


SCRIPT_DIR = Path(__file__).resolve().parent


def atomic_write(path, content):
    """Replace a file only after its new contents have been written."""
    temporary_path = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8",
                                         dir=path.parent, delete=False) as temporary:
            temporary_path = Path(temporary.name)
            temporary.write(content)
        temporary_path.replace(path)
    finally:
        if temporary_path is not None and temporary_path.exists():
            temporary_path.unlink()


def sync_conda_pins(prefix):
    """Track our constraints in separate metadata, never in pin comments."""
    specs = []
    for line in (SCRIPT_DIR / "requirements.txt").read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        spec = "".join(line.split())
        # Match the original installer's filtering of pip-only entries.
        if spec.startswith("-") or "://" in spec or "@" in spec:
            continue
        if spec not in specs:
            specs.append(spec)

    pin_file = prefix / "conda-meta" / "pinned"
    metadata_file = prefix / ".jetset" / "conda-pins.json"
    original = pin_file.read_text() if pin_file.is_file() else ""
    managed = []
    if metadata_file.is_file():
        metadata = json.loads(metadata_file.read_text())
        if (not isinstance(metadata, dict) or metadata.get("version") != 1
                or not isinstance(metadata.get("managed_specs"), list)
                or not all(isinstance(spec, str) for spec in metadata["managed_specs"])):
            raise ValueError("Invalid JetSeT pin metadata: {}".format(metadata_file))
        # Recover if a previous run stopped between the two file replacements.
        if "pending" in metadata:
            pending = metadata["pending"]
            if (not isinstance(pending, dict)
                    or not isinstance(pending.get("before"), str)
                    or not isinstance(pending.get("after"), str)
                    or not isinstance(pending.get("managed_specs"), list)
                    or not all(isinstance(spec, str) for spec in pending["managed_specs"])):
                raise ValueError("Invalid pending JetSeT pin metadata")
            if original == pending["after"]:
                metadata["managed_specs"] = pending["managed_specs"]
            elif original != pending["before"]:
                raise ValueError("Pins changed during an interrupted JetSeT update; "
                                 "inspect {} before retrying".format(metadata_file))
        managed = metadata["managed_specs"]

    retained = original.splitlines()
    # Remove one occurrence for each pin we added, working from the end where
    # we append our pins. Identical pre-existing user pins are never adopted.
    for spec in reversed(managed):
        for index in range(len(retained) - 1, -1, -1):
            if retained[index] == spec:
                del retained[index]
                break

    existing = {"".join(line.split()) for line in retained}
    added = [spec for spec in specs if spec not in existing]
    lines = retained + added
    content = "\n".join(lines) + ("\n" if lines else "")

    metadata_file.parent.mkdir(parents=True, exist_ok=True)
    # Journal the transition so a failed write cannot lose pin ownership.
    pending = {"version": 1, "managed_specs": managed,
               "pending": {"before": original, "after": content,
                           "managed_specs": added}}
    atomic_write(metadata_file, json.dumps(pending, indent=2) + "\n")
    atomic_write(pin_file, content)
    atomic_write(metadata_file, json.dumps(
        {"version": 1, "managed_specs": added}, indent=2) + "\n")
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
            command = [conda, "install", "--yes", "--prefix", str(prefix),
                       "-c", "astropy", "-c", "conda-forge",
                       "--file", "requirements.txt"]
            # Solve against the real environment and its existing pins before
            # touching either the pinned file or our ownership metadata.
            print("Checking JetSeT requirements against existing pins and installed "
                  "Conda packages (dry run). The plan may include package changes.",
                  flush=True)
            try:
                run(command + ["--dry-run"])
            except subprocess.CalledProcessError:
                print("Dependency check failed; pins and JetSeT metadata were not changed.",
                      file=sys.stderr)
                raise
            sync_conda_pins(prefix)
            print("Using {} for requirements.".format(conda), flush=True)
            run(command)
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
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        print("Error while installing JetSeT: {}".format(error), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
