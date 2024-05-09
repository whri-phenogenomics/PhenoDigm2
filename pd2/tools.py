"""Collection of helper functions. 

These are used throughout the other modules.

@author: Tomasz Konopka
"""


import os.path
import sqlite3
from datetime import datetime
from os.path import abspath as abspath
from os.path import join as join
import re


def readHeader(filename, sep=","):
    """Read header definitions from a file."""
    
    with open(filename, "r") as f:
        line = f.readline()
        
        fields = line.strip().split(sep)
        return fields


def time():
    """Return a time string"""
    return datetime.now().strftime('%Y-%m-%d %H:%M:%S')


def matches(pattern, s):
    """Checks whether a pattern fits string."""
    
    return len(re.findall(pattern, s)) > 0


def log(msg, indent=0):
    """Print a log message with a date and indentation."""    
    print("["+time()+"] "+(" "*indent) + msg, flush=True)
    
        
def getPD2dirs(config):
    """Get paths to resources and downloads directories.
    
    Returns paths to four directories (various parts of output locations)
    """
    
    rootdir = abspath(config.db)
    datadir = join(rootdir, "data_raw")
    resourcesdir = join(rootdir, "resources")
    procdir = join(rootdir, "data_processed")
    
    return rootdir, datadir, resourcesdir, procdir
 

def getPD2coredir(config):
    """Get path to a directory for a solr core."""
        
    dbname = abspath(config.db)
    if dbname.endswith("/"):
        dbname = dbname[:-1]
        
    coredir = "phenodigm2_"+os.path.basename(dbname)
    
    allcoresdir = abspath(config.solr_cores_dir)
    coredir = join(allcoresdir, coredir)
    confdir = join(coredir, "conf")
    datadir = join(coredir, "data")    
    
    return coredir, confdir, datadir                
        

def getDBfile(config):
    """Get full path to db file."""

    rootdir = abspath(config.db)
    dbname = os.path.basename(rootdir)
    return join(rootdir, "phenodigm2-"+dbname+".sqlite")        


def getDBconn(dbfile):
    """Get a connection to the databse"""
    
    conn = sqlite3.connect(dbfile, timeout=1800)    
    # This next line enables fetching data by associative array
    conn.row_factory = sqlite3.Row            
    return conn
    

def runProcess(p):
    """Generic function that runs a runnable object."""
    p.run()

