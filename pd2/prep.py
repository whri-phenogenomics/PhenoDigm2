"""Download data for Phenodigm.

@author: Tomasz Konopka
"""

import gzip
import json
import os
import os.path
from os.path import join as join
import requests
import time
import shutil
import xml.etree.ElementTree as XML
from . import tools as pd2tools 


# ############################################################################
# Helper functions for retrieving data

def fetchFromURL(urlpath, filename=None, outdir=None, 
                 skip_if_exists=True,
                 compress=["tab", "txt", "csv", "tsv", "rpt"]):
    """Download from URL into a file.
    
    Arguments:
    urlpath -- the data source URL or query
    filename -- target file name. If None, uses end of urlpath.
        Note: the function uses global resdir as root directory.
    outdir -- directory for download 
    skip_if_exists -- set True to avoid re-download; 
        set False to force re-download even if target already present.
    compress -- list of extensions, function will apply 
        compression when saving files of those types    
    """
    
    if filename is None:
        filename = os.path.basename(urlpath)    
    
    fullpath = os.path.join(outdir, filename)
    logpath = os.path.join(outdir, filename + ".log")
    
    # perhaps create dir if necessary    
    filedir = os.path.dirname(fullpath)
    if not os.path.exists(filedir):                    
        os.makedirs(filedir)
        
    # check whether to compress the output file
    tocompress = False    
    for _ in compress:
        if filename.endswith(_):
            tocompress = True 
    if tocompress:
        filename += ".gz"
        fullpath += ".gz"
    
    # perhaps avoid work                 
    if skip_if_exists and os.path.exists(fullpath):
        pd2tools.log("Skipping: "+filename, 2)
        return
        
    def writeChunks(req, filehandle):
        """Helper to transfer data from a request to file."""
        for chunk in req.iter_content(chunk_size=4096):
            if chunk:
                filehandle.write(chunk)
    
    # download the file
    pd2tools.log("Downloading: "+filename, 2)            
    r = requests.get(url=urlpath, stream=True)
    if tocompress:
        with gzip.GzipFile(fullpath, "wb") as f:
            writeChunks(r, f)                        
    else:
        with open(fullpath, "wb") as f:
            writeChunks(r, f)                    
    r.close()
    
    # write a log file with the download time
    with open(logpath, "w") as f:
        f.write("File\t"+filename+"\n")
        f.write("Source\t"+urlpath+"\n")        
        f.write("Downloaded_on\t"+time.strftime("%c")+"\n")


def fetchUsingXmlQuery(urlpath, querybase, filename=None,
                       querydir=None, outdir=None, skip_if_exists=True):
    """Build a query URL for a REST service, then execute and save.
    
    Arguments:
    urlpath -- base url for REST service
    querybase -- either empty string or nonempty string
            when empty - the urlpath is executed as-is
            when nonempty - function looks for the query content
            in a file querydir/querybase.xml             
    filename -- output file, if blank, uses a default based on querybase
    outdir -- directory for output data 
    skip_if_exists -- passed on to fetchFromURL
    """
        
    querydata = ""
    if querybase != "":    
        # get the query contents from a file on disk
        querypath = os.path.join(querydir, querybase+".xml")
        querydata = []
        with open(querypath, "r") as f:
            for fline in f:
                querydata.append(fline)
        querydata = "".join(querydata)
        querydata = querydata.replace("\n", "")    
        
    if filename is None or filename == "":
        filename = querybase+".txt"
                
    # Then carry out the data fetch    
    fetchFromURL(urlpath+querydata, filename=filename, 
                outdir=outdir, skip_if_exists=skip_if_exists)


# ############################################################################
# Carry out files downloads

def runDirPrep(config):
    """Prepare directory structure."""

    pd2tools.log("Preparing directory structure")

    rootdir, datadir, resourcesdir, dbdir = pd2tools.getPD2dirs(config)
    
    if not os.path.exists(resourcesdir):
        pd2tools.log("Missing directory resources")
        raise Exception("Missing resources directory")
    
    if not os.path.exists(datadir):
        os.makedirs(datadir)
    
    if not os.path.exists(dbdir):
        os.makedirs(dbdir)


def runDownloads(config):
    """Download all the preparation downloads.
    
    This function looks for files declared in 
    dependencies.json and catalog.xml."""
    
    # define the directories for outputs as ../downloads/
    rootdir, datadir, resourcesdir, dbdir = pd2tools.getPD2dirs(config)
        
    pd2tools.log("Downloading data inputs")

    # Parse dependency file for annotations
    dependencies_file = join(resourcesdir, "dependencies.json")
    with open(dependencies_file, "r") as f:
        dependencies = json.load(f)
    
    for nowid in dependencies:
        nowset = dependencies[nowid]
        nowurl = nowset["url"]
        nowdir = os.path.abspath(join(datadir, nowset["targetdir"]))                
        
        # Try to download queries        
        if "query" in nowset:
            nowresdir = join(resourcesdir, nowset["targetdir"])
            for q, f in zip(nowset["query"], nowset["filename"]):
                fetchUsingXmlQuery(nowurl, q, filename=f,
                                   querydir=nowresdir,
                                   outdir=nowdir)
        else:
            # If "query" is not present, download individual files            
            if "filename" in nowset:
                # Fetch the OMIM API from secrets or as a local global variable
                if nowid == 'OMIM':
                    secret = os.getenv('OMIM_API_KEY')
                    nowurl += f'{secret}/' 
                for f in nowset["filename"]:
                     fetchFromURL(nowurl + f, outdir=nowdir)                                    
        
    # Parse dependency file for ontology data
    catalog_file = join(resourcesdir, "catalog.xml")
    catalog_data = XML.parse(catalog_file)    
            
    for xmlnode in list(catalog_data.getroot()):
        if xmlnode.tag.endswith("uri"):
            nowurl = xmlnode.attrib["name"]
            nowfile = xmlnode.attrib["uri"]
            fetchFromURL(nowurl, filename=nowfile, outdir=datadir)

    # Copy some resource files into downloads directory
    for f in ["hp-importer.owl", "mp-importer.owl", "catalog.xml"]:
        shutil.copy(join(resourcesdir, f), join(datadir, f))

