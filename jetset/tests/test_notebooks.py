import importlib.util
import os
import shutil
from pathlib import Path

import pytest

from .base_class import TestBase

nbformat = pytest.importorskip("nbformat")
nbclient = pytest.importorskip("nbclient")
NotebookClient = nbclient.NotebookClient


NOTEBOOK_ROOT_RELATIVE = Path("doc") / "documentation_notebooks" / "notebooks"

# Coverage-target notebook subset based on testing/.coverage analysis.
# These notebooks are the ones that can exercise currently uncovered modules
# (for example: gammapy_plugin, mcmc_ultranest, sherpa_plugin,
# jet_radio_component, poly_fit_tools).
COVERAGE_TARGET_NOTEBOOKS = {
    Path("gammapy_plugin") / "gammapy_plugin.ipynb",
    Path("jet_radio_component") / "Jet_example_radio_component.ipynb",
    Path("model_fit_with_ultranest") / "Jet_example_model_fit_only_ultranest_plain.ipynb",
    Path("phen_constr") / "SSC_th_bkg.ipynb",
    Path("sherpa_plugin") / "sherpa-plugin-sherpa-interface.ipynb",
}


def _find_notebook_root():
    search_starts = [Path.cwd().resolve(), Path(__file__).resolve().parent]
    for start in search_starts:
        for parent in (start, *start.parents):
            candidate = parent / NOTEBOOK_ROOT_RELATIVE
            if candidate.is_dir():
                return candidate
    return None


def _is_hidden_notebook(path):
    return any(part.startswith(".") for part in path.parts)


def _module_available(module_name):
    return importlib.util.find_spec(module_name) is not None


def _collect_notebooks(notebook_root):
    notebooks = []
    for nb_path in sorted(notebook_root.rglob("*.ipynb")):
        rel_path = nb_path.relative_to(notebook_root)
        if _is_hidden_notebook(rel_path):
            continue
        if ".ipynb_checkpoints" in rel_path.parts:
            continue
        if rel_path not in COVERAGE_TARGET_NOTEBOOKS:
            continue
        notebooks.append(nb_path)
    return notebooks


NOTEBOOK_ROOT = _find_notebook_root()
NOTEBOOK_PATHS = [] if NOTEBOOK_ROOT is None else _collect_notebooks(NOTEBOOK_ROOT)
NOTEBOOK_IDS = (
    []
    if NOTEBOOK_ROOT is None
    else [str(path.relative_to(NOTEBOOK_ROOT)) for path in NOTEBOOK_PATHS]
)


def _skip_reason_for_notebook(notebook_path):
    rel_path = notebook_path.relative_to(NOTEBOOK_ROOT)
    parts = rel_path.parts

    if "sherpa_plugin" in parts and not _module_available("sherpa"):
        return "requires optional dependency 'sherpa'"

    if "gammapy_plugin" in parts and not _module_available("gammapy"):
        return "requires optional dependency 'gammapy'"

    if "model_fit_with_ultranest" in parts:
        if not _module_available("ultranest"):
            return "requires optional dependency 'ultranest>=4.0'"
        if not _module_available("h5py"):
            return "requires optional dependency 'h5py'"
        if rel_path.name.endswith("_openmpi.ipynb") and shutil.which("mpirun") is None:
            return "requires 'mpirun' for OpenMPI notebook"

    if (
        "load_data" in parts
        and rel_path.name == "Jet_example_load_data.ipynb"
        and not _module_available("sedbuilder")
    ):
        return "requires optional dependency 'ssdc-sedbuilder'"

    return None


def _execute_notebook(notebook_path, timeout=1800):
    os.environ.setdefault("MPLBACKEND", "Agg")

    with notebook_path.open("r", encoding="utf-8") as file_handle:
        notebook = nbformat.read(file_handle, as_version=4)

    client = NotebookClient(
        notebook,
        timeout=timeout,
        kernel_name="python3",
        resources={"metadata": {"path": str(notebook_path.parent)}},
    )
    client.execute()


class TestNotebooks(TestBase):
    def integration_suite(self, timeout=1800):
        if NOTEBOOK_ROOT is None:
            pytest.skip(
                "unable to locate doc/documentation_notebooks/notebooks from this environment"
            )

        for notebook_path in NOTEBOOK_PATHS:
            skip_reason = _skip_reason_for_notebook(notebook_path)
            if skip_reason:
                continue
            _execute_notebook(notebook_path, timeout=timeout)


@pytest.mark.integration
@pytest.mark.slow
@pytest.mark.parametrize("notebook_path", NOTEBOOK_PATHS, ids=NOTEBOOK_IDS)
def test_execute_notebook(notebook_path):
    if NOTEBOOK_ROOT is None:
        pytest.skip("unable to locate documentation notebooks directory")

    skip_reason = _skip_reason_for_notebook(notebook_path)
    if skip_reason:
        pytest.skip(f"{notebook_path.name}: {skip_reason}")

    _execute_notebook(notebook_path)
