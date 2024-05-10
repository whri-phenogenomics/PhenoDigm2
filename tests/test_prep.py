""" Test suite for the prep module
    Tests the behaviour of phenodigm2.py download 

    @author: Diego Pava
"""
import pytest
import os
import pd2.prep as pd2prep
from pd2.prep import fetchFromURL, fetchUsingXmlQuery, runDirPrep, runDownloads
import gzip
import json
from xml.etree.ElementTree import Element, SubElement, tostring
from xml.dom import minidom
import shutil


# Test for fetchFromURL function

# Define a mock class that behaves like a response to request.get.
# Each function should mock the methods used by response in the fetchFromURL function.
class MockResponse:
    def __init__(self, content):
        self.content = content

    # Mock the behaviour of response.iter_content()
    def iter_content(self, *args, **kwargs):
        yield self.content

    # Mock requests.close() behaviour
    def close(self):
        pass


# Pytest fixture to mock the behaviour of requests.get
# Mock the response as a bytes literal (since fetchFromUrl downloads as a binary stream).
@pytest.fixture
def mock_get(monkeypatch):
    def mock_get_function(*args, **kwargs):
        return MockResponse(b"fake_response")

    # Monkeypatch requests.get to behave like mock_get
    monkeypatch.setattr("requests.get", mock_get_function)

# Test the normal functioning of fetchFromURL
# Parametrize to test compressed and uncompressed files


@pytest.mark.parametrize(
    "filename,tocompress", [("test_file.xml", False), ("test_file.txt", True)]
)
def test_fetch_from_ulr(filename, tocompress, tmp_path, mock_get):
    # Call fetchFromURL
    fetchFromURL("http://fakephenodigmwebiste.com", filename, tmp_path)

    # If it's a compressed test case, adjust filename and open with gzip
    open_func = open
    filelog = filename + ".log"
    if tocompress:
        filename += ".gz"
        open_func = gzip.open

    # Test the file and its log were written
    assert (tmp_path / filename).exists()
    assert (tmp_path / filelog).exists()

    # Test the file is readable and has the expected content
    with open_func(os.path.join(tmp_path, filename), "rb") as f:
        assert f.read() == b"fake_response"


# Tests for fetchUsingXMLQuery

# Test fetchUsingXMLQuery on an empty querybase
def test_fetch_using_xml_query_empty_querybase(tmp_path, mock_get):
    fetchUsingXmlQuery(
        urlpath="http://fakephenodigmwebsite.com",
        querybase="",
        filename="test_file.csv",
        outdir=tmp_path,
    )
    # Check contents of zipped file created.
    with gzip.open(tmp_path / "test_file.csv.gz", "rb") as f:
        assert f.read() == b"fake_response"


# Fixture to create a mock xml file for query testing
@pytest.fixture
def mock_xml_file(tmp_path):
    query_xml_file = tmp_path / "test_query.xml"
    with open(query_xml_file, "w") as file:
        file.write('<Test name = "test_gene_ensembl">')

    # The file is yielded to the test and deleted after.
    yield query_xml_file


# Test fetchUsingXMLQuery when there is a query. 
def test_fetch_using_xml_query_non_empty_querybase(tmp_path, mock_get, mock_xml_file):
    fetchUsingXmlQuery(
        "http://fakephenodigmwebsite.com",
        mock_xml_file.stem,
        "test_file.csv",
        querydir=tmp_path,
        outdir=tmp_path,
    )
    with gzip.open(tmp_path / "test_file.csv.gz", "rb") as f:
        assert f.read() == b"fake_response"


# Tests for runDirPrep

# Test the creation of datadir and dbdir
# The `mock_pd2dirs` fixture mocks a config object for tests.
def test_run_dir_prep_existing_resources_dir(mock_pd2dirs):
    # Create a temporary resources folder
    mock_pd2dirs["resourcesdir"].mkdir()
    # Now the resources directory does exist, so runDirPrep should not raise an Exception
    runDirPrep(mock_pd2dirs)
    assert mock_pd2dirs["datadir"].exists()
    assert mock_pd2dirs["dbdir"].exists()


# Test raising exception when resources is missing
def test_run_dir_prep_missing_resources_dir(mock_pd2dirs):
    with pytest.raises(Exception, match="Missing resources directory"):
        runDirPrep(mock_pd2dirs)


# Test runDownloads

# Fixture to mock `dependencies.json` and `catalog.xml`
@pytest.fixture
def mock_run_downloads_dependencies(monkeypatch, mock_pd2dirs):
    """ Fixture to write mock `dependecies.json` and `catalog.xml` files.

        The fixture takes mock directories from the `mock_pd2dirs` fixture,
        writes the `dependencies.json` and `catalog.xml` files, monkeypatches
        the functions used within `runDownloads` to isolate their behaviour and
        returns the necessary directories.

        Returns: The 4 directories from mock_pd2dirs with necessary download files.

        Notes: The functions used within `runDownloads` are monkeypatched to test
        the r`unDownloads` function irrespective of the its internal functions.
        This is done so if the functionality of the internal methods change, the
        test for `runDownloads` still assesses the function independently.

    """

    # Call the mock_dp2dirs
    rootdir = mock_pd2dirs["rootdir"]
    datadir = mock_pd2dirs["datadir"]
    resourcesdir = mock_pd2dirs["resourcesdir"]
    dbdir = mock_pd2dirs["dbdir"]

    # Create only the mock directories that are used
    datadir.mkdir()
    resourcesdir.mkdir()

    # Create a mock dependencies.json
    dependencies = {
        "id": {
            "url": "http://fakeurl.com",
            "targetdir": "mock_annotations",
            "query": ["query"],
            "filename": ["test_file.xml"],
        }
    }
    # Write mock dependencies file to resourcedir
    dependencies_file = resourcesdir / "dependencies.json"
    with open(dependencies_file, "w") as f:
        json.dump(dependencies, f)

    # Create mock targetdir from dependencies.json
    # runDownloads does this automatically, the test does not.
    (datadir / dependencies["id"]["targetdir"]).mkdir()

    # Create a mock catalog.xml
    root = Element("root")
    uri_element = SubElement(root, "uri")
    uri_element.set("name", "fake_test_url")
    uri_element.set("uri", "test_catalog_file")

    # Write mock catalog file to resourcedir
    catalog_file = resourcesdir / "catalog.xml"
    with open(catalog_file, "w") as f:
        f.write(minidom.parseString(tostring(root)).toprettyxml(indent="\t"))

    # Monkeypatch the functions in runDownloads

    # Monkeypatch fetchUsingXmlQuery to write a simplified file
    def mock_fetchUsingXmlQuery(url, querybase, filename, querydir=None, outdir=None):
        with open(os.path.join(outdir, filename), "wb") as f:
            f.write(b"fake_response")

    monkeypatch.setattr(pd2prep, "fetchUsingXmlQuery", mock_fetchUsingXmlQuery)

    # Monkeypatch fetchFromURL to write a simplified zipped file. 
    def mock_fetchFromURL(url, filename, outdir):
        with open(os.path.join(outdir, filename), "wb") as f:
            f.write(b"fake_response")

    monkeypatch.setattr(pd2prep, "fetchFromURL", mock_fetchFromURL)

    # Monkeypatch shutil.copy to move and write the necessary files.
    def mock_shutil_copy(src, dst):
        filename = src.split("/")[-1]
        if filename in ["hp-importer.owl", "mp-importer.owl", "catalog.xml"]:
            with open(datadir / filename, "w") as f:
                f.write("Mock file")

    monkeypatch.setattr(shutil, "copy", mock_shutil_copy)

    return rootdir, datadir, resourcesdir, dbdir

# TODO: write tests for the OMIM exception
# TODO: write tests to check the secret API key


# Test to check contents of dependencies.json are downloaded correctly
def test_run_downloads(mock_run_downloads_dependencies, mock_pd2dirs):
    rootdir, datadir, resourcesdir, dbdir = mock_run_downloads_dependencies
    runDownloads(mock_pd2dirs)

    # Assert the file from dependencies.json exists in correct dir
    assert (datadir / "mock_annotations" / "test_file.xml").exists()

    # Asser the file from catalog.xml exists
    assert (datadir / "test_catalog_file").exists()

    # Assert importers and catalog exist in the correct directory
    assert (datadir / "hp-importer.owl").exists()
    assert (datadir / "mp-importer.owl").exists()
    assert (datadir / "catalog.xml").exists()
