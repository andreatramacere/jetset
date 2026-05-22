import importlib.util
import inspect

import pytest


def _module_available(module_name):
    return importlib.util.find_spec(module_name) is not None


def _iter_test_methods(instance):
    for name, method in inspect.getmembers(instance, predicate=inspect.ismethod):
        if name.startswith("test_"):
            yield name, method


def _method_kwargs(method):
    kwargs = {}
    signature = inspect.signature(method)
    if "plot" in signature.parameters:
        kwargs["plot"] = False
    return kwargs


@pytest.mark.integration
@pytest.mark.slow
def test_coverage_integration_and_notebooks():
    """Run integration tests and notebook tests from a single entry point."""
    try:
        from .test_integration import TestIntegration
        from .test_notebooks import TestNotebooks
    except Exception as exc:
        pytest.skip(f"unable to import integration/notebook suites: {exc}")

    integration = TestIntegration()
    executed_integration_tests = []

    ultranest_ready = _module_available("ultranest") and _module_available("h5py")

    for name, method in _iter_test_methods(integration):
        if name == "test_ultranest" and not ultranest_ready:
            continue
        method(**_method_kwargs(method))
        executed_integration_tests.append(name)

    assert executed_integration_tests, "no integration test method was executed"

    notebooks = TestNotebooks()
    notebooks.integration_suite(timeout=1800)
