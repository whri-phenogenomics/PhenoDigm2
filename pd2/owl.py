"""Run owltools components.

@author: Tomasz Konopka
(owltools commands from Damian Smedley)
"""

import csv
import gzip
import os
import multiprocessing as mp
from os.path import join as join
from os.path import basename as basename
import subprocess

from . import tools as pd2tools
from . import dbmodels as pd2models

readHeader = pd2tools.readHeader
runProcess = pd2tools.runProcess


# ############################################################################
# Functions to load data from phenodigm-cache into db

def loadOneOwlCacheFile(config, owlpair):
    """Transfer data from one owltools-cache file into db."""
        
    _, _, resourcesdir, dbdir = pd2tools.getPD2dirs(config)
    cachefile = "ontology-embeddings-cache-" + owlpair + ".txt.gz"
    pd2tools.log("Loading: " + cachefile, 2)
    cachefilepath = os.path.join(dbdir, cachefile)
    annodir = os.path.join(resourcesdir, "annotations")

    header_file = os.path.join(annodir, "ontology_embeddings_cache.header")
    header = readHeader(header_file, "\t")
    idindex = header.index("id")
    hitindex = header.index("id_hit")
    csindex = header.index("cs")
    # simindex = header.index("simJ")
    # icindex = header.index("ic")
    # lcsindex = header.index("lcs")
    
    # transfer from disk into memory
    mapdata = pd2models.ModelOntologyOntologyMapping(config.dbfile)      
    with gzip.open(cachefilepath, "rt") as f:
        reader = csv.reader(f, delimiter="\t", quotechar="\"")        
        for fields in reader:
            # idtext = fields[idindex].strip().replace("_", ":")  
            # hittext = fields[hitindex].strip().replace("_", ":")
            # mapdata.addData(idtext, hittext
            mapdata.addData(
                            fields[idindex].strip(),
                            fields[hitindex].strip(),
                            fields[csindex].strip())
                            # fields[simindex].strip(),
                            # fields[icindex].strip(),
                            # fields[lcsindex].strip())
    # transfer into db
    mapdata.save()


def loadOwlCacheFiles(config, owlpairs=("hp-hp", "hp-mp")):
    """Transfer data from owltools-cache-XX-YY.txt files into db."""
    
    for onepair in owlpairs:    
        loadOneOwlCacheFile(config, onepair)


# ############################################################################
# Functions for preparing owltools-cache files

def runMergingOwltools(config, importer_f, instances_f, labels_f,
                       out_f):
    """Run owltools to build merged ontology owl files.
    
    importer_f -- base filename for the importer.owl file
    instances_f -- base filename for --load-instances
    labels_f -- base filename for --load-labels
    out_f -- full path for output file
    """
    
    # perhaps skip work
    if os.path.exists(out_f):
        pd2tools.log("Skipping: "+basename(out_f), 2)
        return
    
    pd2tools.log("Generating " + basename(out_f), 2)
    
    # get file system directories
    t1, datadir, resdir, dbdir = pd2tools.getPD2dirs(config)

    # turn filenames into full paths
    importer_f = join(datadir, importer_f)
    instances_f = join(dbdir, instances_f)
    labels_f = join(dbdir, labels_f)            
    catalog_f = join(datadir, "catalog.xml")
        
    env = os.environ.copy()
    env["OWLTOOLS_MEMORY"] = config.owltools_mem
    
    # generate the owltools command
    ocmd = [config.owltools, "--catalog-xml", catalog_f, importer_f,
            "--merge-imports-closure", "--load-instances", instances_f,
            "--load-labels", labels_f, "--merge-support-ontologies",
            "-o", out_f]
    
    # run owltools and record output into the logfile
    print(" ".join(ocmd))
    with open(out_f + ".log", "w") as f:
        subprocess.call(ocmd, env=env, stdout=f)


class OwlCompProcess:
    """A runnable class for an Owltools onto-comparison job."""
    
    def __init__(self, config, owl_fs, xopt, out_f):
        """Constructor remembers settings.
         
        config -- configuration object, will be used to get dirs
        owl_fs -- list of input owl files (full paths)
        xopt -- options for -x, e.g. HP,HP or MP,HP or HP,MP
        out_f -- output filename
        """
            
        self.config = config
        self.owl_fs = owl_fs
        self.xopt = xopt
        self.out_f = out_f
        
    def run(self):        
        """Run owltools to create 'owltools-cache-X.txt' files."""
        
        # perhaps skip work
        if os.path.exists(self.out_f + ".gz"):
            pd2tools.log("Skipping: "+basename(self.out_f), 2)
            return
    
        pd2tools.log("Generating " + basename(self.out_f), 2)
        
        # get file system directories
        config = self.config
        t1, datadir, resdir, dbdir = pd2tools.getPD2dirs(config)
        
        # will need an updated environment to set java heap size
        env = os.environ.copy()
        env["OWLTOOLS_MEMORY"] = config.owltools_mem
        
        # start building the owltools command 
        ocmd = [config.owltools]
        for f in self.owl_fs:
            ocmd.append(f)
        if len(self.owl_fs) > 1:
            ocmd.append("--merge-support-ontologies")
        
        rest = ["--sim-save-phenodigm-class-scores",
                "-m", config.owltools_min_ic,
                "-x", self.xopt, "-a", self.out_f]
        for r in rest:
            ocmd.append(str(r))            
        
        print(" ".join(ocmd))
        with open(self.out_f+".log", "w") as f:            
            subprocess.call(ocmd, env=env, stdout=f)

        # compress
        if not self.out_f.endswith(".gz"):
            subprocess.call(["gzip", self.out_f])

    
# ############################################################################
# run this from outside the module

def runOwltools(config):
    """Run owltools commands to produce ontology-ontology scores."""
    
    # identify file system directories
    t1, datadir, resdir, dbdir = pd2tools.getPD2dirs(config)

    pd2tools.log("Running owltools")    

    # Hs_f = join(dbdir, "Hs-all.owl")
    # Mm_f = join(dbdir, "Mm-all.owl")
    # cross_f = join(datadir, "obo/upheno/hp-mp/mp_hp-align-equiv.owl")

    # # create merged ontology files
    # runMergingOwltools(config, "hp-importer.owl",
    #                    "Hs-disease-to-phenotype-O.txt",
    #                    "Hs-disease-labels.txt", Hs_f)
    # runMergingOwltools(config, "mp-importer.owl",
    #                    "Mm-gene-to-phenotype-O.txt",
    #                    "Mm-gene-labels.txt", Mm_f)
        
    # # avoid parallel comp because owltools needs a lot of memory
    # numcores = 1

    # # create jobs as runnable class objects    
    # oca = "owltools-cache-"
    # j1 = OwlCompProcess(config, [Hs_f], "HP,HP",
    #                     join(dbdir, oca+"hp-hp.txt"))
    # j2 = OwlCompProcess(config, [Mm_f], "MP,MP",
    #                     join(dbdir, oca+"mp-mp.txt"))
    # j3 = OwlCompProcess(config, [Mm_f, Hs_f, cross_f], "HP,MP",
    #                     join(dbdir, oca+"hp-mp.txt"))    
    # j4 = OwlCompProcess(config, [Mm_f, Hs_f, cross_f], "MP,HP",
    #                     join(dbdir, oca+"mp-hp.txt"))    

    # # create phenodigm-cache files from ontology comparisons
    # with mp.Pool(numcores) as pool:
    #     pool.map(runProcess, [j3, j4, j1, j2])

    # transfer data from owltools cache into db
    pd2tools.log("Transferring owltools data into db")
    loadOwlCacheFiles(config)

