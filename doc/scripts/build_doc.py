#!/usr/bin/env python3

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import contextmanager
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
DOC_ROOT = SCRIPT_DIR.parent
NOTEBOOKS_ROOT = DOC_ROOT / "documentation_notebooks" / "notebooks"
API_ROOT = DOC_ROOT / "api"
DEST_ROOT = DOC_ROOT / "user_guide" / "documentation_notebooks_rst"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build notebook .rst files. "
            "Supports the documented workflow flags: -c (clean), -e (execute), -b (sphinx build)."
        )
    )
    parser.add_argument(
        "-e",
        "--execute",
        action="store_true",
        help="Execute notebooks in-place before .rst conversion.",
    )
    parser.add_argument(
        "-c",
        "--clean",
        action="store_true",
        help=(
            "Clean generated notebook .rst/.png files, then convert notebooks."
        ),
    )
    parser.add_argument(
        "-ncv",
        "--notebook-convert",
        action="store_true",
        help=(
            "Convert Notebooks"
        ),
    )
    parser.add_argument(
        "-b",
        "--build",
        action="store_true",
        help="Run sphinx-build after conversion and update steps.",
    )
    parser.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=None,
        help=(
            "Parallel workers for notebook conversion. "
            "Default: 1 when --execute is used, otherwise min(cpu_count, 10)."
        ),
    )
    parser.add_argument(
        "--sphinx-jobs",
        type=int,
        default=10,
        help="Parallel jobs passed to sphinx-build when --build is enabled (default: 10).",
    )
    parser.add_argument(
        "--clean-only", action="store_true",
        help="Only clean generated notebook .rst/.png files.",
    )
    parser.add_argument(
        "--copy-only", action="store_true",
        help="Only copy existing notebook .rst/.png files and shared PNG images.",
    )
    parser.add_argument(
        "dir_name",
        nargs="?",
        default="",
        help="Optional subdirectory under documentation_notebooks/notebooks.",
    )
    args = parser.parse_args()
    if args.clean_only or args.copy_only:
        if (args.clean_only and args.copy_only) or any(
            (args.clean, args.notebook_convert, args.execute, args.build)
        ):
            parser.error("--clean-only and --copy-only cannot be combined with other workflow actions")
    return args


def run_cmd(cmd: list[str], env: dict[str, str] | None = None) -> None:
    print("$", " ".join(cmd))
    subprocess.run(cmd, cwd=str(DOC_ROOT), check=True, env=env)


def validate_search_dir(dir_name: str) -> Path:
    if dir_name:
        search_dir = NOTEBOOKS_ROOT / dir_name
    else:
        search_dir = NOTEBOOKS_ROOT

    search_dir = search_dir.resolve()
    if not search_dir.is_relative_to(NOTEBOOKS_ROOT.resolve()):
        raise ValueError("Notebook directory must be inside documentation_notebooks/notebooks")
    if not search_dir.is_dir():
        raise FileNotFoundError(f"Directory not found: {search_dir}")

    return search_dir


def is_hidden(path: Path) -> bool:
    return any(part.startswith(".") for part in path.parts)


def collect_notebooks(search_dir: Path) -> list[Path]:
    notebooks: list[Path] = []
    for nb_path in sorted(search_dir.rglob("*.ipynb")):
        rel_path = nb_path.relative_to(DOC_ROOT)
        if ".ipynb_checkpoints" in rel_path.parts:
            continue
        if is_hidden(rel_path):
            continue
        notebooks.append(nb_path)
    return notebooks


def clean_notebook_files(search_dir: Path) -> None:
    destination = DEST_ROOT / "notebooks" / search_dir.relative_to(NOTEBOOKS_ROOT.resolve())
    for root in (search_dir, destination):
        if root.is_dir():
            for path in root.rglob("*"):
                if path.is_file() and path.suffix in (".rst", ".png"):
                    path.unlink()


def copy_generated_files(search_dir: Path) -> None:
    source_root = (DOC_ROOT / "documentation_notebooks").resolve()
    DEST_ROOT.mkdir(parents=True, exist_ok=True)
    for root, suffixes in (
        (search_dir, (".rst", ".png")),
        (source_root / "slides", (".png",)),
        (source_root / "images", (".png",)),
    ):
        if not root.is_dir():
            continue
        for path in root.rglob("*"):
            if path.is_file() and path.suffix in suffixes:
                destination = DEST_ROOT / path.relative_to(source_root)
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(path, destination)


def clean_api_dir_keep_gitkeep() -> None:
    if not API_ROOT.exists():
        return
    if not API_ROOT.is_dir():
        raise RuntimeError(f"Expected directory, found: {API_ROOT}")

    for entry in API_ROOT.iterdir():
        if entry.name == ".gitkeep":
            continue
        if entry.is_dir():
            shutil.rmtree(entry)
        else:
            entry.unlink()


def _truncate_ultranest_stream(lines: list[str], limit: int = 20) -> list[str]:
    out: list[str] = []
    seen_marker = False
    i = 0
    while i < len(lines):
        line = lines[i]
        if not seen_marker and line == "    ====== ultranest script ========":
            seen_marker = True
            block: list[str] = [line]
            i += 1

            while i < len(lines) and len(block) < limit:
                block.append(lines[i])
                i += 1

            out.extend(block)

            skipped = False
            while i < len(lines) and (lines[i].startswith("    ") or lines[i].strip() == ""):
                skipped = True
                i += 1

            if skipped:
                out.append(f"    ... [output truncated: first {limit} lines shown] ...")
                out.append("")
            continue

        out.append(line)
        i += 1
    return out


def postprocess_rst(rst_path: Path) -> None:
    if not rst_path.is_file():
        return

    text = rst_path.read_text(encoding="utf-8")
    original = text

    marker = ".. parsed-literal::\n\n    ====== ultranest script ========"
    replacement = ".. code-block:: text\n\n    ====== ultranest script ========"
    text = text.replace(marker, replacement, 1)

    lines = text.splitlines()
    lines = _truncate_ultranest_stream(lines, limit=20)
    text = "\n".join(lines)
    if original.endswith("\n"):
        text += "\n"

    if text != original:
        rst_path.write_text(text, encoding="utf-8")


def convert_one_notebook(nb_path: Path, execute: bool) -> None:
    rel_nb = nb_path.relative_to(DOC_ROOT)
    print(rel_nb)

    if execute:
        run_cmd(
            [
                "jupyter",
                "nbconvert",
                "--to",
                "notebook",
                "--execute",
                "--inplace",
                str(rel_nb),
            ]
        )

    run_cmd(["jupyter", "nbconvert", str(rel_nb), "--to", "rst"])
    postprocess_rst(nb_path.with_suffix(".rst"))


def convert_notebooks(notebooks: list[Path], execute: bool, jobs: int | None) -> None:
    if jobs is None:
        if execute:
            jobs = 1
        else:
            jobs = max(1, min(os.cpu_count() or 1, 10))
    if jobs < 1:
        raise ValueError("--jobs must be >= 1")

    print(f"Converting {len(notebooks)} notebook(s) with jobs={jobs}, execute={execute}")
    if jobs == 1:
        for nb_path in notebooks:
            convert_one_notebook(nb_path, execute=execute)
        return

    failures: list[tuple[Path, Exception]] = []
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = {
            pool.submit(convert_one_notebook, nb_path, execute): nb_path for nb_path in notebooks
        }
        for future in as_completed(futures):
            nb_path = futures[future]
            try:
                future.result()
            except Exception as exc:
                failures.append((nb_path, exc))

    if failures:
        for nb_path, exc in failures:
            print(f"[ERROR] {nb_path.relative_to(DOC_ROOT)} -> {exc}", file=sys.stderr)
        raise RuntimeError(f"{len(failures)} notebook conversion task(s) failed")


def run_workflow(args: argparse.Namespace) -> None:
    search_dir = validate_search_dir(args.dir_name)
    if args.clean_only:
        clean_notebook_files(search_dir)
        return
    if args.copy_only:
        copy_generated_files(search_dir)
        return
    notebooks = collect_notebooks(search_dir)
    if not notebooks:
        print(f"No notebooks found in {search_dir}")
        return

    clean_api_dir_keep_gitkeep()
    run_cmd([sys.executable, "make_apidoc_and_uml_graphs.py"])
    if args.clean:
        clean_notebook_files(search_dir)

    if args.notebook_convert or args.clean:
        convert_notebooks(notebooks, execute=args.execute, jobs=args.jobs)

    copy_generated_files(search_dir)

    if args.build:
        if args.sphinx_jobs < 1:
            raise ValueError("--sphinx-jobs must be >= 1")
        if shutil_which("sphinx-build") is None:
            raise RuntimeError("sphinx-build command not found")

        with readthedocs_env_true():
            run_cmd([
                "sphinx-build",
                "-j",
                str(args.sphinx_jobs),
                "-b",
                "html",
                "./",
                "build",
            ])


def shutil_which(command: str) -> str | None:
    from shutil import which

    return which(command)


@contextmanager
def readthedocs_env_true():
    previous = os.environ.get("READTHEDOCS")
    os.environ["READTHEDOCS"] = "True"
    try:
        yield
    finally:
        if previous is None:
            os.environ.pop("READTHEDOCS", None)
        else:
            os.environ["READTHEDOCS"] = previous


def main() -> int:
    args = parse_args()
    try:
        run_workflow(args)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
