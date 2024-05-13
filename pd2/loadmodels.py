""" Load mouse data from disk files into the sqlite db.

@author: Tomasz Konopka
"""

import csv
import gzip
import os.path
from os.path import join as join

from . import tools as pd2tools
from . import dbmodels as pd2models
from . import dss as pd2dss
from .dbextractors import *

# shortcuts
readHeader = pd2tools.readHeader
PhenoSet = pd2dss.PhenotypeSet 


# ############################################################################
# Helper functions that process mouse model data files

class MouseModel:
    """Container class for a mouse model."""
    
    # a particular mp term: states no abnormal phenotypes 
    nophenotype = "MP:0002169"
    
    def __init__(self, id, source, strain, description, life_stage="NA"):
        """Create a new mouse model with a description.
        
        By default, the model will have one associated
        phenotype that indicates 'no phenotype'
        """ 
        
        self.id = str(id)        
        self.strain = strain        
        self.description = description
        self.life_stage = life_stage
        self._mp = set()
        self._mp.add(self.nophenotype)
        self._marker = set()
        self._source = set()
        self._source.add(source)
    
    def addMP(self, mpid):
        """Add a new phenotype into current set of phenotypes.
        
        This tries to remove the nophenotype default.
        """
        
        if mpid != "" and mpid != self.nophenotype:
            if self.nophenotype in self._mp:
                self._mp.remove(self.nophenotype)        
            self._mp.add(mpid)

    def addMarker(self, marker):
        self._marker.add(marker)

    def addSource(self, source):
        self._source.add(source)
    
    def nonempty(self, iterable):
        """Avoid empty strings, e.g for MP, markers, etc."""
        
        return [_ for _ in iterable if _ != ""]        
        
    @property
    def source(self):
        return ",".join(self.nonempty(self._source))    
        
    @property
    def mp(self):
        return list(self.nonempty(self._mp))
    
    @property
    def marker(self):
        return list(self.nonempty(self._marker))


def load_life_stages(life_stage_file):
    """create a map from a life-stage accession to a group name/code"""
    
    acc_index, group_index = None, None
    result = dict()    
    with open(life_stage_file, "rt") as f:                            
        for fields in csv.reader(f, delimiter="\t", quotechar="\""):                    
            if acc_index is None:                
                acc_index = fields.index("life_stage_acc")
                group_index = fields.index("life_stage_group")                                                
                continue
            result[fields[acc_index]] = fields[group_index]            
    return result
        
        
def pval(pvalstring):
    """Convert a string into a p value."""
        
    try:
        pval = float(pvalstring)
        if pval > 1:
            pval = 1
        elif pval < 0:
            pval = 0
        return pval
    except ValueError:        
        return 1
        

def saveMouseModels(dbfile, models):
    """Transfer data from an in-memory set of models into db.
    This is used to transfer data from IMPC and MGI loaders."""
    
    model_data = pd2models.ModelModel(dbfile)
    genotype_data = pd2models.ModelModelGenotype(dbfile)
    model_mp_data = pd2models.ModelModelPhenotype(dbfile)
    for nowuid in models:
        nowmodel = models[nowuid]
        model_data.addData(id=nowmodel.id,
                           source=nowmodel.source,
                           species="Mus musculus",
                           genetic_background=nowmodel.strain,
                           life_stage=nowmodel.life_stage,
                           description=nowmodel.description)
        genotype_data.addData(nowmodel.id, nowmodel.marker, nowmodel.description)
        model_mp_data.addData(nowmodel.id, nowmodel.mp)
    model_data.save()
    model_mp_data.save()
    genotype_data.save()


def readIMPCModels(infile, models, addphenotypes=False, pthreshold=1e-4,
                   life_stages=None):
    """Collect information about models from a file
    
    infile -- file to open (must be gziped with certain columns)
    models -- a dict with models
    addphenotypes -- boolean determines whether phenotypes from the file
                     should be recorded
    pthreshold -- p value threshold
    life_stages -- dict with translation between an accession code and life
                   stage code name
    
    returns -- a modified object models                 
    """
            
    with gzip.open(infile, "rt") as f:    
        impc_reader = csv.reader(f, delimiter=",", quotechar="\"")
        impcheader = None
        for fields in impc_reader:
            # learn the table structure from the first line in the file
            if impcheader is None:
                impcheader = fields
                sourceindex = impcheader.index("resource_name")
                strainindex = impcheader.index("strain_name")
                geneindex = impcheader.index("marker_accession_id")
                alleleindex = impcheader.index("allele_symbol")
                modelindex = impcheader.index("allele_accession_id")
                zygindex = impcheader.index("zygosity") 
                mpindex = impcheader.index("mp_term_id")
                pvalindex = impcheader.index("p_value")
                life_index = impcheader.index("life_stage_acc")                                    
                continue                        
            
            # processing for data rows
            zygosity = fields[zygindex][0:3]
            if zygosity == "":
                continue
            if fields[modelindex] == "" or fields[life_index] == "":
                continue

            # create a unique identifier for the mouse model
            life_stage_group = life_stages[fields[life_index]]
            modelid = fields[modelindex] + "#" + zygosity + "#" + life_stage_group
            
            # record a new mouse model               
            if modelid not in models:
                description = fields[alleleindex] + " " + zygosity + " " + life_stage_group
                nowmodel = MouseModel(modelid,
                                      fields[sourceindex], fields[strainindex],
                                      description)
                nowmodel.life_stage = life_stage_group                                   
                models[modelid] = nowmodel   
            
            models[modelid].addMarker(fields[geneindex]) 
            models[modelid].addSource(fields[sourceindex])
            
            # perhaps record a phenotype
            if addphenotypes:
                if pval(fields[pvalindex]) <= pthreshold:
                    models[modelid].addMP(fields[mpindex])
    
    return models


def loadIMPCMouseModelData(db_file, anno_dir, p_threshold, headers_dir):
    """Parse IMPC download and fill tables model and model_phenotype.""" 

    # pick a file to get model definitions from
    sr = "IMPC_ALL_statistical_results"
    file_sr = join(anno_dir, sr + ".csv.gz")
    file_sr_dev = join(anno_dir, sr + "_dev.csv.gz")
    if os.path.exists(file_sr_dev):
        file_sr = file_sr_dev    
    
    # pick a file to get phenotypes from
    gp = "IMPC_ALL_genotype_phenotype"
    file_gp = join(anno_dir, gp + ".csv.gz")
    file_gp_dev = join(anno_dir, gp + "_dev.csv.gz")        
    if os.path.exists(file_gp_dev):
        file_gp = file_gp_dev    
    
    # file with definitions of life stages
    life_stage_file = join(headers_dir, "impc_life_stages.csv")
    life_stages = load_life_stages(life_stage_file)
    
    # scan the statistical_results file to obtain a full list of models
    models = dict()
    models = readIMPCModels(file_sr, models, False, life_stages=life_stages)
    
    # scan the genotype-phenotype association file to get phenotypes
    # set p-value threshold to 1 to catch p_value=NA in this file
    models = readIMPCModels(file_gp, models, True, 1, life_stages=life_stages)
                    
    # transfer into db
    saveMouseModels(db_file, models)                


def isIMPC(s):
    """helper gives true if input points to IMPC source."""
    
    if "EUCOMM" in s or "KOMP" in s:
        return True
    if "IMPC" in s or "NCOM" in s:
        return True
    return False

  
def loadMGIMouseModelData(db_file, anno_dir, headers_dir, db_dir):
    """Fill table with mouse model-phenotype associations."""
    
    # get information about IMPC models (genotypes and phenotypes)
    modelgenes = getDbMapSets(pd2models.ModelModelGenotype(db_file)) 
    modelphen = pd2models.ModelIdPhenotype(db_file)    
    modelphenotypes = getDbMapSets(pd2models.ModelModelPhenotype(db_file))
    
    file_pheno = join(anno_dir, "MGI_GenePheno.rpt.gz")
    header_pheno = join(headers_dir, "MGI_GenePheno.header")

    header = readHeader(header_pheno, ",")        
    geneindex = header.index("MGI_Marker_Accession_ID")    
    phenoindex = header.index("Mammalian_Phenotype_ID")
    bgindex = header.index("Genetic_Background")
    pubmedindex = header.index("PubMed_ID")
    alleleindex = header.index("Allelic_Composition")
    modelindex = header.index("MGI_Genotype_Accession_ID")
    maxindex = max([geneindex, phenoindex, bgindex, 
                    pubmedindex, alleleindex, modelindex])

    # perform two scans of the file
    # in the first pass, identify models that are supported by pubmed links
    model_pmid = set()    
    with gzip.open(file_pheno, "rt") as f:    
        reader = csv.reader(f, delimiter="\t", quotechar="\"")
        for fields in reader:
            # this len check is necessary because MGI omits last columns
            # sometimes and the csv reader provides short lists.
            # accessing elements in the last columns then causes exceptions
            if len(fields) <= maxindex:
                print("omitting raw line: " + str(fields))
                continue                            
            modelid = fields[modelindex].strip()
            if fields[pubmedindex] != "":
                model_pmid.add(modelid)

    # in a second pass, transfer model phenotypes into memory
    # this pass skips over some models that do not have any pubmed references
    models = dict()
    gene_phenotypes = dict()         
    with gzip.open(file_pheno, "rt") as f:    
        reader = csv.reader(f, delimiter="\t", quotechar="\"")
        for fields in reader:            
            nowindex = fields[geneindex].strip()
            # skip entries with strange gene ids
            if nowindex == "" or (' ' in nowindex):
                print("skipping a model with index " + str(nowindex))
                continue
            
            modelid = fields[modelindex].strip()
            nowallele = fields[alleleindex]
            # skip entries that describe IMPC models without additional refs
            if modelid not in model_pmid and isIMPC(nowallele):
                continue

            # perhaps record a new mouse model
            if modelid not in models:
                nowmodel = MouseModel(modelid, "MGI", fields[bgindex], nowallele)
                for g in nowindex.split("|"):
                    nowmodel.addMarker(g)
                models[modelid] = nowmodel

            # record the phenotype
            models[modelid].addMP(fields[phenoindex])
            
            # JAX files can have multiple genes associated
            # with the phenotypes (mouse with more than one aberration)
            # ignore those items
            if len(nowindex.split(",")) > 1:
                continue
            if len(nowindex.split("|")) > 1:
                continue            
                                    
            if nowindex not in gene_phenotypes:
                gene_phenotypes[nowindex] = PhenoSet("")
            gene_phenotypes[nowindex].add(fields[phenoindex])

    # transfer into db
    saveMouseModels(db_file, models)                

    # output gene phenotypes
    out_phenotypes = join(db_dir, "Mm-gene-to-phenotype-O.txt")        
    with open(out_phenotypes, "w") as f:
        for nowg in sorted(gene_phenotypes.keys()):
            for hh in gene_phenotypes[nowg].phenotypes:                        
                f.write(str(nowg) + "\t" + hh + "\n")
    
    # Scan for gene titles
    file_genes = join(anno_dir, "MGI_EntrezGene.rpt.gz")
    header_genes = join(headers_dir, "MGI_EntrezGene.header")
  
    header = readHeader(header_genes, ",")
    geneindex = header.index("MGI_Marker_Accession_ID")
    titleindex = header.index("Marker_Symbol")
    statusindex = header.index("Status")
    
    with gzip.open(file_genes, "rt", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t", quotechar="\"")
        for fields in reader:
            nowindex = fields[geneindex].strip()
            nowstatus = fields[statusindex].strip()
            if nowstatus == "W" or nowindex == "":
                continue
            
            if nowindex in gene_phenotypes:
                gene_phenotypes[nowindex].title = fields[titleindex]
    
    out_labels = join(db_dir, "Mm-gene-labels.txt")
    with open(out_labels, "w") as f:
        for nowg in sorted(gene_phenotypes.keys()):
            nowt = gene_phenotypes[nowg].title
            f.write(str(nowg) + "\t" + nowt + "\n")
    

# ############################################################################
# Run this function from outside module to load all data

def runLoadModels(config):
    """Load model data from disk into db."""
    
    # identify file system directories
    t1, downdir, resdir, dbdir = pd2tools.getPD2dirs(config)

    # load orthologs from the MGI file
    annodir = join(downdir, "annotations")
    headersdir = join(resdir, "annotations")                
    dbfile = config.dbfile
    
    pd2tools.log("Loading: mouse models from IMPC", 2)
    loadIMPCMouseModelData(dbfile, annodir, config.impc_pval, headersdir)
    pd2tools.log("Loading: mouse models from MGI", 2)
    loadMGIMouseModelData(dbfile, annodir, headersdir, dbdir)

