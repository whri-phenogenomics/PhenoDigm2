import pytest
import os
from pathlib import Path
import sqlite3
from pd2.dbmodels import (
    PhenodigmTable,
    PhenodigmSimpleGenerator,
    PhenodigmJoinGenerator,
    ModelGene,
    ModelGeneGeneMapping,
    ModelOntology,
    ModelOntologySynonym,
    ModelDisease,
    ModelDiseaseGeneMapping,
    ModelDiseaseUpdate,
    ModelModel,
    ModelModelGenotype,
    ModelIdPhenotype,
    ModelModelPhenotype,
    ModelDiseasePhenotype,
    ModelOntologyOntologyMapping,
    ModelAssociation,
)
from pd2.dbbuild import runDBBuild


# Create a fixutre to serve as a database for some PhenodigmTable objects
# This will run before each test of the class
@pytest.fixture
def setup_phenodigmtable(tmp_path):
    # Set a PhenodigmTable instance
    dbfile = tmp_path / "test.db"
    tabname = "query_table"
    phenod_table = PhenodigmTable(dbfile, tabname)
    # Define fieldnames as you would in a child classs
    phenod_table.fieldnames = ["id", "query"]

    # Create a connection and a cursor
    with sqlite3.connect(dbfile) as conn:
        cursor = conn.cursor()

        # Create a new table
        cursor.execute(
            f"""
            CREATE TABLE {tabname}(
                id INTEGER PRIMARY KEY,
                query TEXT
            )
        """
        )
        # Add 3 rows to the table (integer and text to test both types)
        cursor.executemany(
            f"INSERT INTO {tabname} VALUES (?, ?)",
            [(1, "query1"), (2, "query2"), (3, "query3")],
        )
    return phenod_table


class TestPhenodimgtable:
    # Test getConn
    def test_getConn(self, setup_phenodigmtable):
        phenod_table = setup_phenodigmtable
        # Create a connection to db through the PhenodigmTable instace.
        conn = phenod_table.getConn()

        # Verify a connection is returned
        assert conn is not None

        # Verify the returned connection is an instance of sqlite3.Connection
        assert isinstance(conn, sqlite3.Connection)

        # Verify row_factory is set to sqlite3.Row
        assert conn.row_factory == sqlite3.Row

        conn.close()

    # Test for countrows method
    def test_countrows(self, setup_phenodigmtable):
        phenod_table = setup_phenodigmtable
        # Call the method and test the expected number of rows
        result = phenod_table.countrows()
        assert result == 3

    # Test for save method
    def test_save(self, setup_phenodigmtable):
        phenod_table = setup_phenodigmtable
        # add values for id and query
        new_data = [4, "query4"]
        phenod_table.data.append(new_data)
        # Run save
        phenod_table.save()
        # Query the db to see if the new values were added.
        with phenod_table.getConn() as conn:
            cursor = conn.cursor()
            cursor.execute(
                f"""SELECT * FROM {phenod_table.tabname}
                        WHERE id = 4"""
            )
            result = cursor.fetchone()

        # Check the added row is not empty and it includes what is expceted
        assert result is not None
        assert result["id"] == 4
        assert result["query"] == "query4"

        # check self data is empty after save
        assert phenod_table.data == []

    # Test for clear
    def test_clear(self, setup_phenodigmtable):
        phenod_table = setup_phenodigmtable
        # Add multiple items into self.data and check they are there
        example_data = [[4, "query4"], [5, "query5"], [6, "query6"]]
        phenod_table.data = example_data
        assert phenod_table.data == example_data

        # Check self.data is empty after running clear.
        phenod_table.clear()
        assert phenod_table.data == []

    # Test for emptyTable
    def test_empty_table(self, setup_phenodigmtable):
        # call the table
        phenod_table = setup_phenodigmtable
        # run command
        phenod_table.emptyTable()

        # check it is empty
        with phenod_table.getConn() as conn:
            cursor = conn.cursor()
            cursor.execute(f"SELECT * FROM {phenod_table.tabname}")
            result = cursor.fetchall()

        assert result == []

    # TODO: Test for getall --> not used in PhenoDigm?
    # The function also seems incomplete


# Test for generator
class TestPhenodigmSimpleGenerator:
    # Test test_next
    def test_next(self, setup_phenodigmtable):
        # Define a PSG object, this is an iteratable
        test_generator = PhenodigmSimpleGenerator(setup_phenodigmtable)

        # Iterate over generator and test the contents.
        # The function returns row by row, hence why we iterate.
        assert [row["id"] for row in test_generator.next()] == [1, 2, 3]
        assert [row["query"] for row in test_generator.next()] == [
            "query1",
            "query2",
            "query3",
        ]


# Test for join generator
# The class needs two PhenodigmTable objects, a second fixture is created:
# TODO: this needs to be less repetitive
@pytest.fixture
def second_phenodigmtable(tmp_path):
    # Set a PhenodigmTable instance
    dbfile = tmp_path / "test.db"
    tabname = "match_table"
    phenod_table = PhenodigmTable(dbfile, tabname)
    # Define fieldnames as you would in a child classs
    phenod_table.fieldnames = ["id", "match"]

    # Create a connection and a cursor
    with sqlite3.connect(dbfile) as conn:
        cursor = conn.cursor()

        # Create a new table
        cursor.execute(
            f"""
            CREATE TABLE {tabname}(
                id INTEGER PRIMARY KEY,
                match TEXT
            )
        """
        )
        # Add 3 rows to the table (integer and text to test both types)
        cursor.executemany(
            f"INSERT INTO {tabname} VALUES (?, ?)",
            [(1, "match1"), (2, "match2"), (3, "match3")],
        )
    return phenod_table


class TestPhenodigmJoinGenerator:
    # Creates a join generator instance from two phenodigmtables before each test.
    @pytest.fixture(autouse=True)
    def setup_tables(self, setup_phenodigmtable, second_phenodigmtable):
        self.pdtable_1 = setup_phenodigmtable
        self.pdtable_2 = second_phenodigmtable
        self.join_generator = PhenodigmJoinGenerator(
            self.pdtable_1, self.pdtable_2, ["id", "id"]
        )

    # Test fieldsql
    def test_fieldsql(self):
        # Generate field names calling the tested function
        fields = self.join_generator.fieldsql(
            self.pdtable_1.fieldnames, self.pdtable_1.tabname
        )

        # Assert field names are as expected SQL aliases.
        assert fields == [
            "query_table.id AS query_table_id",
            "query_table.query AS query_table_query",
        ]

    # Test next method by iterating over the yielded rows
    def test_next(self):
        # Capture the output from next() and store it as a list
        rows = list(self.join_generator.next())

        # Iterate over the elements of next() and assert expected values
        assert [x["query_table_id"] for x in rows] == [1, 2, 3]
        assert [x["query_table_query"] for x in rows] == ["query1", "query2", "query3"]
        assert [x["match_table_id"] for x in rows] == [1, 2, 3]
        assert [x["match_table_match"] for x in rows] == ["match1", "match2", "match3"]


# Tests for the remaining classes that insert that into the database
# All classes require an existing database with corresponding tables
# Below, fixtures are created to set up the database for each test class.


# Prearing the database
# Return a path to a mock database for the whole module.
@pytest.fixture(scope="function")
def db_path(tmp_path_factory):
    """
    Fixture to return a path to a test database for the whole module.

    The fixutre uses `mp_path_factory` to create a tmp_dir that
    can be shared accross tests, as opposed to `tmp_path`, which
    is created and deleted for individual tests.

    Returns: a path for a test database.
    """
    # Define a temporary database
    return tmp_path_factory.mktemp("database") / "test.db"


# Create a mock db and tables to test all dbmodels.
@pytest.fixture(scope="function")
def mock_db_and_tables(db_path):
    """
    Fixture to build a mock database with pd2tools.DBBuild

    The fixture creates a single database that persists for
    all tests within the module instead of creating a database
    when each test runs. It is created using the actual db_schema.json.

    Does not return anything. Only creates necessary database to
    build instances of the data structures.

    """

    # MockConfig object for building necessary databases
    class MockConfig:
        def __init__(self, db, datadir, resource, dbfile, action):
            self.db = db
            self.datadir = datadir
            self.resource = resource
            self.dbfile = dbfile
            self.action = action

    # Define the resource dir to fetch db_schema
    resource_dir = Path.cwd() / "resources"

    # Define a mock configuration to satisfy the requirements of DBBbuild
    mock_config = MockConfig("", "", resource_dir, db_path, "build")

    # Run the database build and tables
    runDBBuild(mock_config)


# The tests require building an instance of the class and saving data to the db
# To avoid repetition, creating and instace and saving are passed as fixtures.
@pytest.fixture
def setup_model(request, db_path, mock_db_and_tables, test_data):
    # Take Model class from passed parameter
    model_class = request.param
    # Select database and create model instance
    model_instance = model_class(db_path)
    # Return the instance
    return model_instance


# Tests for ModelGene class
class TestModelGene:
    # Test data to be added to ModelGene instance
    @pytest.fixture
    def test_data(self):
        return {
            "gene_id": "test_id",
            "organism": "test_organism",
            "symbol": "test_symbol",
            "name": "test_name",
            "altname": "test_altname",
            "type": "test_type",
            "locus": "test_locus",
            "withdrawn": 1,
        }

    @pytest.mark.parametrize("setup_model", [ModelGene], indirect=True)
    def test_model_gene_add_data(self, test_data, setup_model):
        # Add the data unpacking kwargs
        setup_model.addData(**test_data)
        # Assert the content in the data attribute matches the order and values.
        assert setup_model.data == [list(test_data.values())]

    # Test for empty gene_id or symbol
    @pytest.mark.parametrize("setup_model", [ModelGene], indirect=True)
    @pytest.mark.parametrize("empty_param", [("gene_id"), ("symbol")])
    def test_model_gene_empty_gene_or_symbol(self, test_data, setup_model, empty_param):
        # Replace id or symbol with empty value
        test_data[empty_param] = ""

        # Checks adding empty params returns None
        assert setup_model.addData(**test_data) is None

        # Checks the data attribute is empty after calling addData
        assert setup_model.data == []


# Test ModelGeneGeneMapping
class TestModelGeneGeneMapping:
    # Test data to be added
    @pytest.fixture
    def test_data(self):
        return {"query": "test_query", "match": "test_match"}

    @pytest.mark.parametrize("setup_model", [ModelGeneGeneMapping], indirect=True)
    def test_model_gene_gene_mapping_add_data(self, test_data, setup_model):
        # Add the data unpacking kwargs
        setup_model.addData(**test_data)
        # Assert the content in the data attribute matches the order and values.
        assert setup_model.data == [list(test_data.values())]

    # Test for empty query or match
    @pytest.mark.parametrize("setup_model", [ModelGeneGeneMapping], indirect=True)
    @pytest.mark.parametrize("empty_param", [("query"), ("match")])
    def test_model_ggm_empty_query_or_match(self, test_data, setup_model, empty_param):
        # Replace id or symbol with empty value
        test_data[empty_param] = ""

        # Checks adding empty params returns None
        assert setup_model.addData(**test_data) is None

        # Checks the data attribute is empty after calling addData
        assert setup_model.data == []


class TestModelOntology:
    # Test data to be added
    @pytest.fixture
    def test_data(self):
        return {"id": "test_id", "term": 'test_term"'}

    @pytest.mark.parametrize("setup_model", [ModelOntology], indirect=True)
    def test_model_ontology(self, test_data, setup_model):
        # Add the data unpacking kwargs
        setup_model.addData(**test_data)

        # Assert the input data is the same as stored the data attribute
        # Note the change at test_term' is as expected.
        assert setup_model.data == [["test_id", "test_term'"]]


# TODO: Comment this better
class TestModelOntologySynonym:
    # Test data to be added
    # Synonym takes sets, lists or strings
    @pytest.fixture
    def test_data(self, request):
        return {"id": "test_id", "synonym": request.param}

    # Different values of 'synonym' to be passed. Note the quotation marks "" used.
    @pytest.mark.parametrize(
        "test_data",
        [
            ['"synonym1"', '"synonym 2"', '"syn onym3"'],
            {'"synonym1"', '"synonym 2"', '"syn onym3"'},
            '"synonym1"',
        ],
        indirect=True,
    )
    @pytest.mark.parametrize("setup_model", [ModelOntologySynonym], indirect=True)
    def test_model_ontology_synonym(self, test_data, setup_model):
        # Add the data unpacking kwargs
        setup_model.addData(**test_data)

        # Conditional assertion: one case for strings, another for the rest.
        # Assertion to check elements are added correctly,
        # and correct quotation marks are replaced ''.
        if isinstance(test_data["synonym"], str):
            assert setup_model.data == [["test_id", "'synonym1'"]]
        else:
            assert sorted(setup_model.data) == sorted(
                [
                    ["test_id", "'synonym1'"],
                    ["test_id", "'synonym 2'"],
                    ["test_id", "'syn onym3'"],
                ]
            )


# Test for ModelDisease
class TestModelDisease:
    # Test data to add, note keys are parameters of addData function,
    # not the class fieldnames.
    @pytest.fixture
    def test_data(self):
        return {
            "id": "test_id",
            "term": "test_term",
            "alts": "test_alts",
            "cname": "test_class",
        }

    # TODO: group of repeated tests.
    @pytest.mark.parametrize("setup_model", [ModelDisease], indirect=True)
    def test_model_disease(self, test_data, setup_model):
        # Add the data unpacking kwargs
        setup_model.addData(**test_data)

        # Assert the input data is the same as stored the data attribute
        assert setup_model.data == [list(test_data.values())]


# Tests for ModelDiseaseGeneMapping
class TestModelDiseaseGeneMapping:
    # Test data to add
    @pytest.fixture
    def test_data(self):
        return {
            "query": "test_query",
            "match": "test_match",
            "locus": "test_locus",
            "source": "test_source",
        }

    # TODO: group of repeated tests.
    @pytest.mark.parametrize("setup_model", [ModelDiseaseGeneMapping], indirect=True)
    def test_model_disease_gene_mapping(self, test_data, setup_model):
        # Add the data unpacking kwargs
        setup_model.addData(**test_data)

        # Assert the input data is the same as stored the data attribute
        assert setup_model.data == [list(test_data.values())]

    # Test for empty query or match
    @pytest.mark.parametrize("setup_model", [ModelDiseaseGeneMapping], indirect=True)
    @pytest.mark.parametrize("empty_param", [("query"), ("match")])
    def test_model_dgm_empty_query_or_match(self, test_data, setup_model, empty_param):
        # Replace id or symbol with empty value
        test_data[empty_param] = ""

        # Checks adding empty params returns None
        assert setup_model.addData(**test_data) is None

        # Checks the data attribute is empty after calling addData
        assert setup_model.data == []


# Need to readdress the database for this test.
class TestModelDiseaseUpdate:
    # insert mock data into `disease` table
    @pytest.fixture
    def disease_update_prep(self, db_path):
        model_disease = ModelDisease(db_path)
        data = {
            "id": "1",
            "term": "term_1",
            "alts": "alts_1",
            "cname": "class_1",
        }
        model_disease.addData(**data)
        model_disease.save()

        return model_disease

    # Define updated data to be inserted or not.
    # Different `id` values can be passed and a static `term` value.
    @pytest.fixture
    def test_data(self, request):
        return {"id": request.param, "term": "term_2"}

    # Test MDU class. pass different `id` values.
    @pytest.mark.parametrize("setup_model", [ModelDiseaseUpdate], indirect=True)
    @pytest.mark.parametrize("test_data", ["1", "2"], indirect=True)
    def test_model_disease_update(
        self, test_data, setup_model, disease_update_prep, db_path
    ):
        # Helper function to query the disease table by id.
        def fetch_by_id(id):
            with sqlite3.connect(db_path) as conn:
                cursor = conn.cursor()
                cursor.execute("SELECT * FROM disease where id=?", (id,))
                data = cursor.fetchone()
            return data

        # The disease_update_prep runs automatically, disease table has content.
        # The ModelDiseaseUpdate attempts to add new data
        setup_model.addData(**test_data)

        # The database is queried
        output = fetch_by_id(test_data["id"])

        # Case for 1
        # Data should not be added, expect original test_data
        if test_data["id"] == "1":
            assert output == ("1", "term_1", "alts_1", "class_1")
        # Case for 2
        # Data should be added, new id and new term added.
        else:
            assert output == ("2", "term_2", None, None)


# Test for ModelModel class
class TestModelModel:
    # Define test data
    @pytest.fixture
    def test_data(self):
        return {
            "id": "test_id",
            "source": "test_source",
            "species": "test_species",
            "genetic_background": "test_backgroud",
            "life_stage": "test_life_stage",
            "description": "test_description",
        }

    # TODO: this is also from the repetitive series
    # Test for ModelModel
    @pytest.mark.parametrize("setup_model", [ModelModel], indirect=True)
    def test_model_model(self, test_data, setup_model):
        # Add the data unpacking kwargs
        setup_model.addData(**test_data)

        # Assert the input data is the same as stored the data attribute
        assert setup_model.data == [list(test_data.values())]



class TestModelModelGenotype:
    # Define test data
    # Gene can be a str or a list
    @pytest.fixture
    def test_data(self, request):
        return {
            "id": "test_id",
            "gene": request.param,
            "description": "test_description",
        }

    # Test for ModelModelGenotype
    @pytest.mark.parametrize("setup_model", [ModelModelGenotype], indirect=True)
    @pytest.mark.parametrize(
        "test_data",
        [["gene_list1", "gene_list2", "gene_list3"], "gene_str1"],
        indirect=True,
    )
    def test_model_model_genotype(self, test_data, setup_model):
        # Add the data unpacking kwargs
        setup_model.addData(**test_data)

        # Conditional assertion: one case for strings, another for the lists.
        # Assertion to check elements are added correctly

        # If string, expect a single list with 3 terms
        # For lists, expect different 3 lists with each gene.
        if isinstance(test_data["gene"], str):
            assert setup_model.data == [["test_id", "gene_str1", "test_description"]]
        else:
            assert setup_model.data == [
                ["test_id", "gene_list1", "test_description"],
                ["test_id", "gene_list2", "test_description"],
                ["test_id", "gene_list3", "test_description"],
            ]


class TestModelIdPhenotype:
    # Define test data
    # phenotype can be a str or a list
    @pytest.fixture
    def test_data(self, request):
        return {"id": "test_id", "phenotype": request.param}

    
    @pytest.mark.parametrize("setup_model", [ModelIdPhenotype], indirect=True)
    @pytest.mark.parametrize(
        "test_data",
        [["phenotype_list1", "phenotype_list2", "phenotype_list3"], "phenotype_str1"],
        indirect=True,
    )
    # Test for ModelIdPhenotype
    def test_model_id_phenotype(self, setup_model, test_data):
        # Add the data unpacking kwargs
        setup_model.addData(**test_data)

        # Conditional assertion: one case for strings, another for the lists.
        # Assertion to check elements are added correctly

        # If string, expect a single list with 2 terms
        # For lists, expect different 3 lists with each phenotype
        if isinstance(test_data["phenotype"], str):
            assert setup_model.data == [["test_id", "phenotype_str1"]]
        else:
            assert setup_model.data == [
                ["test_id", "phenotype_list1"],
                ["test_id", "phenotype_list2"],
                ["test_id", "phenotype_list3"],
            ]


# Test ModelModelPhenotype
class TestModelModelPhenotype:
    def test_model_model_phenotype(self, db_path):
        # Setup: instantiate the class
        mmp = ModelModelPhenotype(db_path)

        # Assert the tab name is correct
        assert mmp.tabname == "model_phenotype"

    # Check inheritability?


# Test ModelDiseasePhenotype
class TestModelDiseasePhenotype:
    def test_model_disease_phenotype(self, db_path):
        # Setup: instantiate the class
        mdp = ModelDiseasePhenotype(db_path)

        # Assert the tab name is correct
        assert mdp.tabname == "disease_phenotype"

    # Check inheritability?


# Test ModelOntologyOntologyMapping
class TestModelOntologyOntologyMapping:
    # Define test data
    @pytest.fixture
    def test_data(self, request):
        return {
            "query": request.param[0],
            "match": request.param[1],
            "simJ": "test_simJ",
            "ic": "test_simj",
            "lcs": "test_lcs",
        }

    @pytest.mark.parametrize(
        "setup_model", [ModelOntologyOntologyMapping], indirect=True
    )
    @pytest.mark.parametrize(
        "test_data",
        [
            # Case for query < match and query != match
            ("TS_0000001", "TS_0005000")
        ],
        indirect=True,
    )
    def test_model_onto_onto_map_match_smaller_than_query_and_unequal(
        self, test_data, setup_model
    ):
        setup_model.addData(**test_data)

        assert setup_model.data == [
            ["TS_0000001", "TS_0005000", "test_simJ", "test_simj", "test_lcs"],
            ["TS_0005000", "TS_0000001", "test_simJ", "test_simj", "test_lcs"],
        ]
        print(setup_model.data)

    @pytest.mark.parametrize(
        "setup_model", [ModelOntologyOntologyMapping], indirect=True
    )
    @pytest.mark.parametrize(
        "test_data",
        [
            # Case for match < query
            ("TS_0005000", "TS_0000001"),
        ],
        indirect=True,
    )
    def test_model_onto_onto_mapping_match_smaller_than_query(
        self, test_data, setup_model
    ):
        setup_model.addData(**test_data)

        assert setup_model.data == []

    @pytest.mark.parametrize(
        "setup_model", [ModelOntologyOntologyMapping], indirect=True
    )
    @pytest.mark.parametrize(
        "test_data",
        [
            # Case for match == query?
            ("TS_0005000", "TS_0005000")
        ],
        indirect=True,
    )
    def test_model_onto_onto_map_equal_match_query(self, test_data, setup_model):
        setup_model.addData(**test_data)

        assert setup_model.data == [
            ["TS_0005000", "TS_0005000", "test_simJ", "test_simj", "test_lcs"]
        ]

    @pytest.mark.parametrize(
        "setup_model", [ModelOntologyOntologyMapping], indirect=True
    )
    @pytest.mark.parametrize(
        "test_data",
        [
            # Case for match !=query
            ("TS_0005000", "TS_0004900"),
        ],
        indirect=True,
    )
    def test_model_onto_onto_map_unequal_match_query(self, test_data, setup_model):
        setup_model.addData(**test_data)

        assert setup_model.data == []


# Test ModelAssociation
class TestModelAssociation:
    @pytest.fixture
    def test_data(self):
        return {
            "q": "test_q",
            "m": "test_m",
            "s_avg_norm": 0.1,
            "s_avg_raw": 1.5,
            "s_max_norm": 0.9,
            "s_max_raw": 2.7,
            "p1": "test_p1",
            "p2": "test_p2",
        }

    # Test for addData. This is not used in PhenoDigm to date.
    @pytest.mark.parametrize("setup_model", [ModelAssociation], indirect=True)
    def test_model_association_add_data(self, test_data, setup_model):
        setup_model.addData(**test_data)

        assert setup_model.data == [list(test_data.values())]

    @pytest.mark.parametrize("setup_model", [ModelAssociation], indirect=True)
    def test_model_association_add_data_array(self, setup_model):
        test_arr = ["test_data", "should", "be", "appended"]
        setup_model.addDataArray(test_arr)

        assert setup_model.data == [test_arr]
