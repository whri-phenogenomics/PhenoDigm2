"""Utilities to export phenodigm tables into tab-separated format. 

@author: Tomasz Konopka
"""


import json
from os.path import join
from . import tools as pd2tools


# ############################################################################
# Helper functions

def printOneTable(conn, tabname, colnames, where=""):
    """Print contents of a table to stdout.
    
    conn -- a db connection
    tabname -- name of table to query
    colnames -- list with colnames to extract
    where -- optional string for sql where clause
    """
    
    # start with the header line
    print("\t".join(colnames))
    
    # select from db and print
    sql = "SELECT " + ", ".join(colnames) + " FROM " + tabname
    if where != "": 
        sql += " WHERE "+where
                    
    c = conn.cursor()
    c.execute(sql)
    for row in c:
        result = []
        for k in colnames:                        
            result.append(str(row[k]))
        rowstr = "\t".join(result)
        print(rowstr)        


# ############################################################################
# call this from external file

def exportTables(config):
    """Checks input, then prints contents of a table."""

    if config.table == "":
        print("Missing required option --table")
        return
    
    # read db schema and check that requested table exists
    rootdir, datadir, resourcesdir, dbdir = pd2tools.getPD2dirs(config)              
    schemafile = join(resourcesdir, "db_schema.json")
    with open(schemafile, "r") as f:
        dbschema = json.load(f)    
    if config.table not in dbschema.keys():
        print("Unknown table name: "+config.table)
        return
        
    # extract columns from the db schema definition
    tableconfig = dbschema[config.table]
    colnames = [f.split(" ")[0] for f in tableconfig["fields"]]
    
    # perform db operations in separate function
    conn = pd2tools.getDBconn(config.dbfile)    
    printOneTable(conn, config.table, colnames, config.where)    
    conn.close()

