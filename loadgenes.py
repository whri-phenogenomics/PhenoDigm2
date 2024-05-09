""" Load gene data from disk files into the sqlite db.

@author: Tomasz Konopka
"""

import csv
import gzip
from os.path import join as join

from . import tools as pd2tools
from . import dbmodels as pd2models
from . import parsers as pd2parsers
from . import dss as pd2dss

readHeader = pd2tools.readHeader


# ############################################################################
# Intermediate functions for reading/parsing/merging data

def loadGeneData(dbfile, annodir, headersdir):    
    """Load plain gene definitions from MGI into db.
    
    Returns a set of ids that have been recorded into the db.
    """    
        
    file_mgi = join(annodir, "MGI_EntrezGene.rpt.gz")     
    header_mgi = join(headersdir, "MGI_EntrezGene.header")
    file_hgnc = join(annodir, "hgnc_complete_set.txt.gz")
    header_hgnc = join(headersdir, "hgnc_complete_set.header")
    file_hgnc2 = join(annodir, "withdrawn.txt.gz")
    header_hgnc2 = join(headersdir, "withdrawn.header")
   
    # use a custom parser to read from file into a memory object
    parser = pd2parsers.GeneParser()
    parser.parseMGI(file_mgi, header_mgi)
    parser.parseHGNCWithdrawn(file_hgnc2, header_hgnc2)
    parser.parseHGNC(file_hgnc, header_hgnc)
    
    # transfer information into db
    gene_data = pd2models.ModelGene(dbfile)
    for geneid in parser.genes.keys():
        genelike = parser.genes[geneid]
        # add the official and secondary gene symbols
        official = genelike.official_symbol
        if official != "" and official is not None:
            gene_data.addData(geneid, genelike.organism, 
                              symbol=official,
                              name=genelike.name, 
                              altname=genelike.altname, 
                              type=genelike.type,
                              locus=genelike.locus, withdrawn=False)
        # add unofficial symbols
        for unofficial in genelike.symbols:
            if unofficial == official:
                continue
            gene_data.addData(geneid, genelike.organism, 
                              symbol=unofficial, name="-", 
                              altname=genelike.altname, type="-", 
                              locus="-", withdrawn=True)
    
    # actually send to db
    gene_data.save()
        
    return set(parser.genes.keys())    
        
 
def readEnsemblOrthologs(annodir, headersdir, genes):
    """Read a file from Ensemble and extract orthologs.
    
    returns: dictionary with orthologs.
    """
                
    file_mh = join(annodir, "human_mouse_mapping.txt.gz")     
    header_mh = join(headersdir, "human_mouse_mapping.header")
    
    # figure out what data is in what position in the table
    mhheader = readHeader(header_mh, "\t")
    HGNCindex = mhheader.index("hgnc_id")
    HGNCsymbol = mhheader.index("hgnc_symbol")
    MGIindex = mhheader.index("mgi_id")
    MGIsymbol = mhheader.index("mgi_symbol")
    
    # transfer data from data file into memory objects
    orthologs = dict()    
    with gzip.open(file_mh, "rt") as f:
        for line in f:        
            fields = line[:-1].split("\t")            
            # save definition of gene ids and symbols
            mgi = fields[MGIindex]            
            hgnc = fields[HGNCindex]
            if mgi == "" or hgnc == "":
                continue
            
            # save association between mouse-human 
            if mgi not in orthologs:                
                orthologs[mgi] = set()       
            orthologs[mgi].add(hgnc)
            
            if mgi not in genes or hgnc not in genes:
                mh = "(" + mgi + ", " + hgnc + ")"
                print("Reading Ensembl orthologs, did not understand "+mh)                
            
    return orthologs


def readMGIOrthologs(annodir, headersdir, genes):
    """Read a file from MGI and extract orthologs.
    
    This file presents orthologs in multiple lines, 
    e.g. an ortholog relation between geneA and geneX would be
    x mouse geneA
    x human geneX
    
    annodir, headersdir - paths to directories on disk
    genes - dict linking gene ids to pd2dss.Gene objects    
    """
    
    # load orhtologs from MGI
    file_hom = join(annodir, "HOM_MouseHumanSequence.rpt.gz")     
    header_hom = join(headersdir, "HOM_MouseHumanSequence.header")
    
    homheader = readHeader(header_hom, ",")
    orthog_index = homheader.index("HomoloGene_ID")
    organism_index = homheader.index("Common_Organism_Name")
    mouseid_index = homheader.index("Mouse_MGI_ID")
    humanid_index = homheader.index("HGNC_ID")    
     
    orthologset = dict()
    with gzip.open(file_hom, "rt") as f:
        reader = csv.reader(f, delimiter="\t", quotechar="\"")
        for fields in reader:                    
            organism = fields[organism_index].split(",")[0]            
            orthog = fields[orthog_index]
            
            if orthog not in orthologset:
                orthologset[orthog] = pd2dss.OrthologMapping()            
            if organism == "mouse":
                geneid = fields[mouseid_index]                
            elif organism == "human":
                geneid = fields[humanid_index]
            else:
                continue
            
            if geneid == "":
                continue
            
            if geneid not in genes:
                print("MGI MouseHumanSequence contains new geneid: "+geneid)
                print("in pair: "+geneid+" "+orthog) 
                                                                   
            orthologset[orthog].addPair(organism, geneid)
        
    return orthologset
    
    
def loadOrthologData(dbfile, annodir, headersdir, genes):
    """Transfer mouse-human gene associations.
    
    dbfile - file with target db
    annodir, headersdir - paths to directories on disk
    genes - set with ids that have been defined in genes table.
    """
            
    # load orthologs from Ensembl and from MGI
    orthologs = readEnsemblOrthologs(annodir, headersdir, genes)        
    orthologset = readMGIOrthologs(annodir, headersdir, genes)
        
    # transfer ortholog definitions into db   
    orthologs_data = pd2models.ModelGeneGeneMapping(dbfile)    
    for oo in orthologset.keys():
        oodata = orthologset[oo].getPairs("mouse", "human")
        for x in oodata:
            i, j = x[0], x[1]
            if i not in orthologs:
                orthologs[i] = set()
            if j not in orthologs[i]:                
                orthologs[i].add(j)
    
    for mouseid in orthologs.keys():
        for humanid in orthologs[mouseid]:            
            orthologs_data.addData(mouseid, humanid) 
            orthologs_data.addData(humanid, mouseid)
                       
    orthologs_data.save()


# ############################################################################
# Run this function from outside module to load all data

def runLoadGenes(config):
    """Load gene-centered data into db."""
    
    # identify file system directories
    t1, downdir, resdir, dbdir = pd2tools.getPD2dirs(config)

    # load orthologs from the MGI file
    annodir = join(downdir, "annotations")
    headersdir = join(resdir, "annotations")                
    dbfile = config.dbfile
    
    # get gene definitions, use gene defs during ortholog loading
    pd2tools.log("Loading: gene definitions", 2)
    genes = loadGeneData(dbfile, annodir, headersdir)
    pd2tools.log("Loading: gene orthologs", 2)  
    loadOrthologData(dbfile, annodir, headersdir, genes)

