"""
Tests for contents of pd2.tools.py
(functions with generic tools)

@author: Diego Pava
"""

import pytest
from unittest.mock import patch, MagicMock
import os
from datetime import datetime
import sqlite3
from pd2.tools import (
    readHeader,
    time,
    matches,
    log,
    getPD2dirs,
    getPD2coredir,
    getDBfile,
    getDBconn,
    runProcess,
)


# Test for readHeader
# Defines a fixture that creates a tmp file that readHeader should read
@pytest.fixture
def mock_header_file(tmp_path):
    headers = "Header_1,Header_2,Header_3"
    path = tmp_path / "mock_header"

    with open(path, "w") as f:
        f.write(headers)

    return path


# Test readHeader produces expected content from mock_header_file
def test_read_Header(mock_header_file):
    assert readHeader(mock_header_file) == ["Header_1", "Header_2", "Header_3"]


# Monkeypatch datetime.now method to return the same time during testing.
@pytest.fixture
def mock_datetime_now(monkeypatch):
    class TestDatetime:
        @classmethod
        def now(cls):
            return datetime(2023, 6, 29)

    monkeypatch.setattr("pd2.tools.datetime", TestDatetime)


# Test the time function returns expected time.
def test_time(mock_datetime_now):
    assert time() == "2023-06-29 00:00:00"


# TODO?: Test for matches --> seems unused by the software
# Test cases for matches function
@pytest.mark.parametrize(
    "pattern, string, expected_result",
    [
        (r"\d+", "Pheno2023", True),  # pattern matches the string
        (r"\d+", "Pheno", False),  # pattern doesn't match the string
        (r"\d+", "", False),  # string is empty
        ("", "Pheno", True),  # pattern is empty
    ],
)
def test_matches(pattern, string, expected_result):
    assert matches(pattern, string) == expected_result


# Test the log function
# mock_datetime_now is passed to have a static time in the log output.
# This is because log calls time()
def test_log(capsys, mock_datetime_now):
    # Set arguments for log()
    test_message = "This is a test message"
    indent = 4

    # Call log function
    log(test_message, indent)

    # Capture log's output
    output = capsys.readouterr()

    # Assert the log message has expceted content and structure
    assert output.out == f"[2023-06-29 00:00:00] {' ' * indent}{test_message}\n"


# Class to mock the config argument
class MockConfig:
    def __init__(self, db):
        self.db = db


# Test for pd2Dirs
def test_get_pd2Dirs(tmp_path):
    # Set path for a temporary database
    dbpath = tmp_path / "test_db"

    # Create a mock config object
    mock_config = MockConfig(dbpath)

    # Call the function with the mock config
    _, datadir, resourcesdir, procdir = getPD2dirs(mock_config)

    # Assert paths are what is expected
    assert datadir == os.path.join(dbpath, "data_raw")
    assert resourcesdir == os.path.join(dbpath, "resources")
    assert procdir == os.path.join(dbpath, "data_processed")


# Test for getPD2coredir
# Database names to test as parameters
@pytest.mark.parametrize("dbname", ["test_db/", "test_db"])
def test_get_PD2coredir(dbname, tmp_path):
    # Set a config object with database path
    mock_config = MockConfig(dbname)

    # Add the solr_cores_dir attribute and assign it to the tmp_path
    mock_config.solr_cores_dir = tmp_path / "test_solr_cores_dir"

    # Call function with the mock config
    coredir, confdir, datadir = getPD2coredir(mock_config)

    # Assert trailing "/" was removed when naming the coredir
    assert not coredir.endswith("/")
    # Assert paths are what is expected
    assert coredir == os.path.join(mock_config.solr_cores_dir, "phenodigm2_test_db")
    assert confdir == os.path.join(coredir, "conf")
    assert datadir == os.path.join(coredir, "data")


# Test for getDBfile
def test_get_DBfile(tmp_path):
    dbpath = tmp_path / "test_db"
    mock_config = MockConfig(dbpath)

    dbfile = getDBfile(mock_config)

    assert dbfile == os.path.join(dbpath, "phenodigm2-test_db.sqlite")


# Basic getDBconn test to ensure right arguments are being called.
def test_getDBconn():
    with patch("sqlite3.connect") as mock_connect:
        dbfile = "test_db.sqlite"
        getDBconn(dbfile)

        # Check the right arguments were called.
        mock_connect.assert_called_once_with(dbfile, timeout=1800)

        # Grab the mock connection object that would be created
        mock_conn = mock_connect.return_value

        # Check if row_factory was set correctly
        assert mock_conn.row_factory == sqlite3.Row


# Test to check getDBconn fetches data by associative array
def test_getDBconn_row_factory():
    dbfile = ":memory:"  # use an in-memory database for testing
    conn = getDBconn(dbfile)

    cursor = conn.cursor()
    # Create a test table
    cursor.execute("CREATE TABLE test (id INTEGER, name TEXT)")

    # Insert a test row
    cursor.execute("INSERT INTO test VALUES (?, ?)", (1, "PhenoDigm2"))

    # Query the database
    cursor.execute("SELECT * FROM test WHERE id = 1")
    row = cursor.fetchone()

    # Check that the row is an sqlite3.Row instance
    assert isinstance(row, sqlite3.Row)

    # Check that we can access elements by column name
    assert row["id"] == 1
    assert row["name"] == "PhenoDigm2"

    conn.close()


# Test for runProcess
def test_runProcess():
    # Create a mock object with a "run" method
    mock_p = MagicMock()

    # Call the function to be tested
    runProcess(mock_p)

    # Assert that the run method was called
    mock_p.run.assert_called_once()


# This seems like test code, does not actually test anything from the tools module

# from os.path import join, dirname
# import unittest
# from pd2 import tools as pd2tools
# from pd2 import loaddiseases as pd2diseases


# s2c = pd2diseases.omimTitle

# class ToolsTests(unittest.TestCase):
#     """Test for pd2 tools."""

#     def setUp(self):
#         """Setup for individual tests."""
#         pass

#     def test_captilize_simple(self):
#         input1 = "title in lowercase"
#         expected1 = "Title In Lowercase"
#         input2 = "LONG TITLE IN UPPER"
#         expected2 = "Long Title In Upper"
#         self.assertEqual(expected1, s2c(input1), "lower to caps")

#     def test_capitalize_hyphens(self):
#         input1 = "abc-xyz; GENE"
#         expected1 = "Abc-Xyz"
#         input2 = "ABC-XYZ; GENEabc"
#         expected2 = "Abc-Xyz"
#         self.assertEqual(expected1, s2c(input1), "lower to caps")
#         self.assertEqual(expected2, s2c(input2), "upper to caps")
