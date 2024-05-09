"""Build PhenoDigm2 database.

@author: Tomasz Konopka
"""

import json
import os
from os.path import join as join
from . import tools as pd2tools 


# ############################################################################
#

class PhenoDigmDB:
    """For creating tables for the PhenoDigm sqlite DB."""

    def __init__(self, conn):
        """Initialize object to create tables"""
                        
        self.conn = conn            
        self.cursor = self.conn.cursor()
        
    def maketable(self, tabname, fields):
        """Create a new table and indexes in the DB.
        
        Arguments:
        tabname -- name of the table
        fields -- list with the required fields
        indexes -- dict defining indexes to create on table
        """                
                
        # Create table                
        sqli = ", ".join(fields)            
        sql = "CREATE TABLE " + tabname + " (" + sqli + ")"
        self.cursor.execute(sql)                
        self.conn.commit()

    def makeindexes(self, tabname, indexes):
        """Create indexes on an existing table."""
        
        # use in-memory indexing
        self.cursor.execute("PRAGMA temp_store=2")
        
        # (re-)create indexes
        for key, val in indexes.items():
            indexname = tabname+"__"+key
            drop = "DROP INDEX IF EXISTS "+indexname
            self.cursor.execute(drop)
            sql = "CREATE INDEX " + indexname + " on "
            sql += tabname + " (" + val + ")"
            self.cursor.execute(sql)
            self.conn.commit()


# ############################################################################
# Create PhenoDigm DB and its tables

def runDBBuild(config):
    """Set up a SQLite DB, return a db connection."""
    
    pd2tools.log("Creating database")
    
    # define the location of the database
    rootdir, datadir, resdir, dbdir = pd2tools.getPD2dirs(config)        
    schema_file = join(resdir, "db_schema.json")
        
    # remove the db (i.e. drop all tables if exists)
    try:
        if config.action in ("build", "new"):
            pd2tools.log("Removing existing database file", 2)
            os.remove(config.dbfile)
    except OSError:
        pass
    
    # obtain a connection to interact with the db
    conn = pd2tools.getDBconn(config.dbfile)
            
    # get definitions for the db
    with open(schema_file, "r") as f:
        db_schema = json.load(f)
        
    pd2tools.log("Creating tables", 2)
    # Create tables.
    db = PhenoDigmDB(conn)
    for table_name in db_schema.keys():
        if table_name.endswith("deprecated"):
            continue        
        table_schema = db_schema[table_name]
        db.maketable(table_name, table_schema["fields"])

    conn.close()


def runDBIndexing(config):
    """Create indexes for DB tables."""
    
    pd2tools.log("Indexing tables")
    
    # define the location of the database
    rootdir, datadir, resdir, dbdir = pd2tools.getPD2dirs(config)        
    schemafile = join(resdir, "db_schema.json")
            
    conn = pd2tools.getDBconn(config.dbfile)
                
    with open(schemafile, "r") as f:
        dbschema = json.load(f)
        
    P = PhenoDigmDB(conn)        
    for tablename in dbschema.keys():
        if tablename.endswith("deprecated"):
            continue        
        tabschema = dbschema[tablename]
        P.makeindexes(tablename, tabschema["indexes"])    

    conn.close()

