import pytest

def pytest_addoption(parser):
    parser.addoption(
        "--workFlow",
        action="store",
        metavar="NAME",
        help="run validation tests for only this workflow (nanopolish|medaka)",
    ),
    parser.addoption(
        "--numValidations",
        action="store",
        default=-1,
        help="the number of validation test datasets to run (default= -1 (all))")

def pytest_configure(config):
    # register an additional marker
    config.addinivalue_line(
        "markers", "env(name): mark test to run only on named environment"
    )

def pytest_runtest_setup(item):
    envnames = [mark.args[0] for mark in item.iter_markers(name="env")]
    if envnames:
        if item.config.getoption("--workFlow") not in envnames:
            pytest.skip("test requires env in {!r}".format(envnames))

@pytest.fixture(scope='session')
def numValidations(request):
    value = request.config.option.numValidations
    if value is None:
        value = -1
    return int(value)