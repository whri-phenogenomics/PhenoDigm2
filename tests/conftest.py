""" Configuration file to share fixtures between modules

@author: Diego Pava
"""
import pytest
from pd2 import tools as pd2tools


# Creates a mock config and patches the pd2tools.getPD2dirs function.
@pytest.fixture
def mock_pd2dirs(monkeypatch, tmp_path):
    """General fixture to mock the config object in PhenoDigm2 and getPD2dirs.

    Creates mock paths for the necessary directories for most functions
    through `config`.
    Uses pytest's monkeypatch to replace the functionality of `pd2tools.getPD2dirs`.
    It returns the mock temporary paths instead of creating new ones,
    these tmp_paths are cleaned after each test.

    Returns: a dictionary containing the temporary paths.
    """
    # Mock config paths
    rootdir = tmp_path / "mock_rootdir"
    datadir = tmp_path / "mock_datadir"
    resourcesdir = tmp_path / "mock_resourcesdir"
    dbdir = tmp_path / "mock_dbdir"

    # Mock get PD2dirs function
    def mock_get_pd2dirs(*args, **kwargs):
        return (rootdir, datadir, resourcesdir, dbdir)

    monkeypatch.setattr(pd2tools, "getPD2dirs", mock_get_pd2dirs)

    return {
        "rootdir": rootdir,
        "resourcesdir": resourcesdir,
        "datadir": datadir,
        "dbdir": dbdir,
    }
