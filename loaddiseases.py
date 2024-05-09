""" Load disease data from disk files into the sqlite db.

@author: Tomasz Konopka
"""

import csv
import gzip
import json
from os.path import join as join

from . import tools as pd2tools
from . import dbmodels as pd2models
from . import parsers as pd2parsers
from . import dss as pd2dss

# shortcuts
readHeader = pd2tools.readHeader
PhenoSet = pd2dss.PhenotypeSet 


# ############################################################################
# Helper functions that parse/process disease data files

def omimTitle(s):
    """Change an Omim disease name to avoid all-caps"""
      
    sparts = s.split(";")
    return sparts[0].title()
    
            
def loadDiseaseData(dbfile, annodir, headersdir):
    """Transfer from OMIM downloads into table disease."""
                    
    # load disease classifications from orphanet
    file_orphaclasses = join(headersdir, "orphanet_classes.json")
    with open(file_orphaclasses,"r") as f:
        class2file = json.loads(f.read())
    odcp = pd2parsers.OrphanetClassificationParser()
    for oneclass in class2file.keys():
        onefile = join(annodir, class2file[oneclass])
        odcp.parseFile(onefile, oneclass)    

    # load the disease defintions from orphanet
    file_orphanet = join(annodir, "en_product1.xml")
    orphanet = pd2parsers.OrphanetDiseaseParser(file_orphanet)
    file_orphagenes = join(annodir, "en_product6.xml")
    orphagenes = pd2parsers.OrphanetGeneParser(file_orphagenes)
        
    # load disease definitions from omim
    file_mimtitles = join(annodir, "mimTitles.txt.gz")     
    header_mimtitles = join(headersdir, "mimTitles.header")
    file_morbidmap = join(annodir, "morbidmap.txt.gz")     
    header_morbidmap = join(headersdir, "morbidmap.header")
    file_mim2gene = join(annodir, "mim2gene.txt.gz")
    header_mim2gene = join(headersdir, "mim2gene.header")                
    omim = pd2parsers.OmimParser()
    omim.parseMim2gene(file_mim2gene, header_mim2gene)
    omim.parseTitles(file_mimtitles, header_mimtitles)
    omim.parseMorbid(file_morbidmap, header_morbidmap)
    
    # transfer from OMIM into the db
    ddata = pd2models.ModelDisease(dbfile)
    for nowmim in omim.mimdata.keys():
        nowdata = omim.mimdata[nowmim]
        nowtitle, nowalt = nowdata[0], nowdata[1]                
        # perhaps use the orphanet classifications, synonyms
        nowclass = ""
        if nowmim in orphanet.synonyms:
            noworpha = orphanet.synonyms[nowmim]
            if noworpha in odcp.classes:
                nowclass = ",".join(list(odcp.classes[noworpha]))
            
        nowtitle = omimTitle(nowtitle)
        ddata.addData(id=nowmim, term=nowtitle, alts=nowalt, cname=nowclass)

    # transfer from ORPHANET into db
    for noworpha in orphanet.titles.keys():        
        nowtitle = orphanet.titles[noworpha]
        nowclass = ""
        if noworpha in odcp.classes:
            nowclass = ",".join(list(odcp.classes[noworpha]))
        ddata.addData(id=noworpha, term=nowtitle, alts="", cname=nowclass)
         
    ddata.save()
    
    # transfer disease-gene associations
    human2id = getSymbol2Id(dbfile, "human")    
    dgmap_data = pd2models.ModelDiseaseGeneMapping(dbfile)
    for nowmim in omim.symbolsloci.keys():
        # keep track of gene ids recorded (avoid duplicate entries)
        nowids = set()        
        for symbollocus in omim.symbolsloci[nowmim]:
            nowsymbol = symbollocus.query
            nowlocus = symbollocus.locus
            if nowsymbol in human2id:                
                nowid = human2id[nowsymbol]
                if nowid not in nowids:
                    dgmap_data.addData(nowmim, nowid, nowlocus, "OMIM, HGNC")
                    nowids.add(nowid)
            
    for noworpha in orphagenes.genelocus.keys():
        for gl in orphagenes.genelocus.get(noworpha):
            dgmap_data.addData(noworpha, gl.query, gl.locus, "ORPHANET")
                                                    
    dgmap_data.save()


def getSymbol2Id(dbfile, organism):
    """Query DB and obtain a mapping from active symbols to ids.
    
    This function can link withdrawn symbols to current ids.    
    """
    
    result = dict()
    sql = "SELECT id, symbol, withdrawn FROM gene WHERE organism=?"
    
    conn = pd2tools.getDBconn(dbfile)
    c = conn.cursor()                
    c.execute(sql, (organism, ))        
    for row in c.fetchall():
        rowsymbol = row["symbol"]        
        if rowsymbol in result:
            # the symbol is already defined
            # overwrite only if this entry is valid            
            if row["withdrawn"] == 0:
                result[rowsymbol] = row["id"]
        else:            
            result[rowsymbol] = row["id"]   
                
    return result


def loadDiseaseGenes(dbfile, annodir, headersdir):
    """Load associations between disease ids and gene ids."""
            
    # read mapping from gene symbols to gene ids
    human2id = getSymbol2Id(dbfile, "human")
        
    file_do = join(annodir, "MGI_DO.rpt.gz")
    dgmap_data = pd2models.ModelDiseaseGeneMapping(dbfile) 
    with gzip.open(file_do, "rt") as f:    
        reader = csv.DictReader(f, delimiter="\t", quotechar="\"")
        for fields in reader:            
            id = fields["DO Disease ID"]
            omim = fields["OMIM IDs"].split("|")
            taxon = fields["NCBI Taxon ID"]
            symbol = fields["Symbol"]
            mgi = fields["Mouse MGI ID"]
            
            geneid = ""
            if taxon == "9606":  # human
                # must convert symbol to HGNC id
                if symbol in human2id:
                    geneid = human2id[symbol]                    
            elif taxon == "10090":  # mouse
                # the mgi id is the gene id
                geneid = mgi
                
            if geneid != "":
                dgmap_data.addData(id, geneid, "", "MGI")
                for oneomim in omim:
                    dgmap_data.addData(oneomim, geneid, "", "MGI")
            
    dgmap_data.save()
    

def checkDiseaseIndex(did):
    """Check that a disease_id, e.g. OMIM:1234 is well formed."""

    # disease ids cannot be empty or with spaces
    if did == "" or (' ' in did):
        return False

    # don't want commas
    if ',' in did:
        return False

    # ids must start with a prefix and have meaningful
    okprefix = ["OMIM:", "DECIPHER:", "ORPHA:"]
    for p in okprefix:
        if did.startswith(p) and len(did) > len(p):
            return True
    return False
    
          
def loadDiseasePhenotypeData(dbfile, annodir, headersdir, dbdir):
    """Fill table with disease-hp associations."""

    file_pheno = join(annodir, "phenotype.hpoa")
    header_pheno = join(headersdir, "phenotype_annotation.header")
    
    header = readHeader(header_pheno, "\t")
    diseaseindex = header.index("DatabaseID")
    titleindex = header.index("DiseaseName")
    qualifierindex = header.index("Qualifier")
    phenoindex = header.index("HPO_ID")

    # Scan file, transfer all annotations into a memory dict
    disease_phenotypes = dict()
    diseases = pd2models.ModelDiseaseUpdate(dbfile) 
    with open(file_pheno, "rt") as f:
        reader = csv.reader(f, delimiter="\t", quotechar="\"")
        for fields in reader:
            if fields[0][0] == "#":
                continue
            if fields[qualifierindex] == "NOT":
                continue
            nowindex = fields[diseaseindex].strip()
            if not checkDiseaseIndex(nowindex):                
                continue
            if nowindex not in disease_phenotypes:
                nowtitle = fields[titleindex].strip()
                disease_phenotypes[nowindex] = PhenoSet(nowtitle)
                diseases.addData(nowindex, nowtitle)
            disease_phenotypes[nowindex].add(fields[phenoindex])
        
    # Transfer info from memory dict into a db table
    disease_phenotype_data = pd2models.ModelDiseasePhenotype(dbfile)
    for nowd in sorted(disease_phenotypes.keys()):        
        nowhps = disease_phenotypes[nowd].phenotypes
        disease_phenotype_data.addData(nowd, nowhps)            
    disease_phenotype_data.save()
    
    # This block (until end of function) produces text files
    # Those text files are eventually used by owltools
    out_diseases = join(dbdir, "Hs-disease-labels.txt")
    with open(out_diseases, "w") as f:
        for nowd in sorted(disease_phenotypes.keys()):
            nowtitle = disease_phenotypes[nowd].title
            f.write(str(nowd)+"\t"+nowtitle+"\n")
    
    out_phenotypes = join(dbdir, "Hs-disease-to-phenotype-O.txt")        
    with open(out_phenotypes, "w") as f:
        for nowd in sorted(disease_phenotypes.keys()):
            for hh in disease_phenotypes[nowd].phenotypes:                        
                f.write(str(nowd)+"\t"+hh+"\n")


# ############################################################################
# Run this function from outside module to load disease data

def runLoadDiseases(config):
    """Load disease data from disk into db."""
    
    # identify file system directories
    t1, downdir, resdir, dbdir = pd2tools.getPD2dirs(config)

    # load orhtologs from the MGI file
    annodir = join(downdir, "annotations")
    headersdir = join(resdir, "annotations")                
    dbfile = config.dbfile
    
    pd2tools.log("Loading: disease definitions", 2)
    loadDiseaseData(dbfile, annodir, headersdir)
    pd2tools.log("Loading: disease-gene mappings", 2)
    loadDiseaseGenes(dbfile, annodir, headersdir)
    pd2tools.log("Loading: disease-hp mappings", 2)   
    loadDiseasePhenotypeData(dbfile, annodir, headersdir, dbdir)    

