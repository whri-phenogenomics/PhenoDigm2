""" Test suite for the dbbuild module
    Tests the behaviour of phenodigm2.py build and index.

    @author: Diego Pava
"""

import pytest
import sqlite3
import json
import os
from pd2.dbbuild import PhenoDigmDB, runDBBuild, runDBIndexing


# Test PhenoDigmDB's maketable and makeindexes methods.
@pytest.fixture
def mock_db_no_schema():
    """Creates a mock SQLite database in memory and creates an indexed table.

    The database is yielded to the test functions, once done, the connection
    to the database is closed and the database is deleted automatically.

    Returns: An instance of PhenoDigmDB representing a mock database.
    """

    # Create a connection in memory and a PDdb instance
    conn = sqlite3.connect(":memory:")
    db = PhenoDigmDB(conn)

    # Create a table from tablename and fields
    tabname = "test_table"
    fields = ["test_id TEXT", "test_score REAL"]
    db.maketable(tabname, fields)

    # Create the indexes here
    indexes = {"t_i": "test_id"}
    db.makeindexes(tabname, indexes)

    # Yield db and close connection
    yield db
    conn.close()


# Test for maketable
def test_maketable(mock_db_no_schema):
    # Obtain general information from testable
    table_info = mock_db_no_schema.cursor.execute(
        "PRAGMA table_info('test_table')"
    ).fetchall()

    # Get the columns and types
    columns = [x[1] for x in table_info]
    types = [x[2] for x in table_info]

    # Assert the column name and type in table are correct
    assert "test_id" and "TEXT" in (columns[0] and types[0])
    assert "table_score" and "REAL" in (columns[1] and types[1])


# Test for makeindexes
def test_makeindexes(mock_db_no_schema):
    # Obtain index information from test_table
    index_info = mock_db_no_schema.cursor.execute(
        "PRAGMA index_list('test_table')"
    ).fetchall()

    # Assert the index name matches the expecteation
    assert "test_table__t_i" in index_info[0][1]


# Test main runDBBuild and runDBIndexing functions


# Mock config class to customise action and dbfile arguments
class MockConfig:
    def __init__(self, action, dbfile):
        self.action = action
        self.dbfile = dbfile


# Define a fixture that sets up and returns dirs with a schema file.
@pytest.fixture
def mock_dirs_and_schema(mock_pd2dirs):
    """Creates mock resource directory from mock_pd2dirs
    and a test db_schema.json file.

    Writes a test schema file to the mock resource directory for other
    fixtures to build databases from it instead of using
    the actual codebase.

    Returns: mock_pd2directories with mock resourcedir and schema file.
    """
    # Create resource directory
    resourcesdir = mock_pd2dirs["resourcesdir"]
    resourcesdir.mkdir()

    # Create a mock schema.json
    schema = {
        "test_test_association": {
            "fields": [
                "query TEXT",
                "match TEXT",
                "score_avg_norm REAL",
                "score_avg_raw REAL",
            ],
            "indexes": {"tta_q": "query"},
        }
    }

    # Path and write mock schema json file
    schema_file = resourcesdir / "db_schema.json"
    with open(schema_file, "w") as f:
        json.dump(schema, f)

    return mock_pd2dirs


# Define a fixture that generates a db from a schema, builds tables and indexes.
@pytest.fixture
def mock_db_with_schema(mock_dirs_and_schema, tmp_path, request):
    """Builds a SQLite database from a schema.json file
    and provides a connection to perform tests.

    Uses custom mock config object. The `action` argument is parametrized
    so multiple values can be passed at each test. The `dbfile` argument is
    a mock `.db` file under a temporary directory. The db is created using
    runDBBuild and indexed with runDBIindexing. Provides a cursor for testing.

    Returns: A cursor connected to a database built from a schema in `.json`
    """
    # Define a tmp_path for the db
    tmp_db = str(tmp_path / "test.db")

    # Define action parameters
    action = request.param

    # Define mock config
    config = MockConfig(action, tmp_db)

    # Build database and indexes
    runDBBuild(config)

    # Index the database
    runDBIndexing(config)

    # Connect to the created database to check contents
    conn = sqlite3.connect(tmp_db)
    cursor = conn.cursor()

    # Provides the test a cursor to perform tests
    yield cursor

    # Closes connection and removes database
    conn.close()
    os.remove(tmp_db)


# Fixture to test the runDBBuild edge cases without yielding cursor or indexing.
@pytest.fixture
def config_fixture(mock_dirs_and_schema, tmp_path, request):
    """Creates a fixture containing a customisable mock config object.

    The action param in the config object can be customised to pass
    multiple options to a test. The dbfile is set to a static test.db
    file.

    Returns: A config object for testing
    """
    # Define action parameter and tmp_db file path
    action = request.param
    tmp_db = str(tmp_path / "test.db")

    # Define mock config
    config = MockConfig(action, tmp_db)

    return config


# Tests the os.remove in runDBBuild when `build` or `new` are passed.
@pytest.mark.parametrize("config_fixture", ["build", "new"], indirect=True)
def test_run_DBBuild_param_build_new(config_fixture, monkeypatch):
    # Declare the config from the fixture
    config = config_fixture

    # Set a mock os.remove for runDBBuild
    # If os.remove is called on the db, it raises an exception instead of removing a file.
    def mock_remove(path):
        assert path == config.dbfile
        raise Exception("The file was mock deleted")

    # Monkeypatch the functionality of os.remove
    monkeypatch.setattr(os, "remove", mock_remove)

    # Run DBBuild, expect the specific exception if os.remove to be called.
    with pytest.raises(Exception, match="The file was mock deleted"):
        runDBBuild(config)


# Tests runDBBuild when neither `build` or `new` are passed.
@pytest.mark.parametrize(
    "config_fixture", ["not_build", "not_new", "other"], indirect=True
)
def test_run_DBBuild_param_not_build_not_new(config_fixture):
    # Declare the config from the fixture
    config = config_fixture

    # Run the function
    runDBBuild(config)

    # Check database exists:
    assert os.path.exists(config.dbfile)


# Tests the OSError is triggered on a non-existent database.
@pytest.mark.parametrize("action", ["build", "new"])
def test_run_DBBuild_raises_OSError(tmp_path, action, mock_pd2dirs):
    # Mock config with build and non-existent database
    config = MockConfig(action, str(tmp_path / "non-existent.db"))

    # OSError is raised if the file does not exist
    with pytest.raises(OSError):
        runDBBuild(config)


# Tests the content of a table created from runDBBuild
@pytest.mark.parametrize("mock_db_with_schema", ["build"], indirect=True)
def test_run_DBBuild(mock_db_with_schema):
    """  Takes a cursor and checks the table information created after runDBBuild.

        Tests the runDBBuild expected behaviour by checking the table
        was created correctly from a `.json` file. Checks the the table exists,
        the number of columns, the column names, and column types.

    """
    # Get the cursor from the fixture
    cursor = mock_db_with_schema

    # Check the table exists
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
    assert cursor.fetchone() is not None

    # Check the table is empty
    cursor.execute("SELECT * FROM 'test_test_association'")
    assert cursor.fetchone() is None

    # Obtain information about the created table
    table_info = cursor.execute("PRAGMA table_info('test_test_association')").fetchall()

    # Get the columns and types
    columns = [x[1] for x in table_info]
    types = [x[2] for x in table_info]

    # Assert column names in table
    assert "query" and "match" and "score_avg_norm" and "score_avg_raw" in columns
    # Assert number of columns
    assert len(table_info) == 4
    # Assert the type of columns
    assert "TEXT" and "REAL" in types


# Tests indexes were created properly from the runDBIndexing function
@pytest.mark.parametrize("mock_db_with_schema", ["index"], indirect=True)
def test_run_DBNIndexing(mock_db_with_schema):
    """  Takes a cursor and checks the index information created after runDBIndexing.

        Tests the runDBIndexing expected behaviour by checking the table
        was created correctly from a `.json` file. Checks the index
        in the json file exists within the database.

    """

    # Get the cursor from the fixture
    cursor = mock_db_with_schema

    # Obtain index information from test_table
    index_info = cursor.execute("PRAGMA index_list('test_test_association')").fetchall()

    # Assert the index name matches the expecteation
    assert "test_test_association__tta_q" in index_info[0][1]
