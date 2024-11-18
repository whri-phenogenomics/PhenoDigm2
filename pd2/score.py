"""Final phase of PhenoDigm2 db build. Generates scores from phenotypes.

Uses owltools-cache files to produce phenodigm scores.

@author: Tomasz Konopka
"""

import time
import gzip
import math
import multiprocessing as mp
import os
import os.path
from . import tools as pd2tools
from . import dbmodels as pd2models
from . import dss as pd2dss
from . import scoring as pd2scoring
from . import dbextractors as pd2extractors

# shorthand
runProcess = pd2tools.runProcess


# ############################################################################
# produce summary tables
    
def fillOOMap(dbfile, source_ontology, target_ontology):
    """Create an OntologyOntologyMap object from db."""
    
    ootab = "ontology_ontology_mapping"
    sql = "SELECT query, match, simJ, ic, lcs FROM "+ootab
    sql += " WHERE query LIKE '"+source_ontology+"%'"
    sql += " AND match LIKE '"+target_ontology+"%'"

    oomap = pd2dss.OntologyOntologyMap()

    conn = pd2tools.getDBconn(dbfile)        
    cur = conn.cursor()
    cur.execute(sql)       
    for row in cur:                
        oomap.add(row)        
    conn.close()
    
    return oomap


def fillIdPhenotypeMap(dbfile, tablename):
    """Create an IdPhenotypeMap object from db."""
    
    map = pd2dss.IdPhenotypeMap()    
    sql = "SELECT id, phenotype FROM "+tablename
    conn = pd2tools.getDBconn(dbfile)    
    cur = conn.cursor()
    cur.execute(sql)    
    for row in cur:        
        map.add(row["id"], row["phenotype"])
    conn.close()
        
    return map        
    

def showOccasionalUpdate(counter, interval=10, total=100, prefix="", time=None):
    """Helper function for writing status messages."""
    if counter % interval == 0 or counter == total:
        counterperc = round(100*counter/total)        
        if counterperc == 100 and time is not None:
            minutes = str(round(time/60))            
            pd2tools.log(prefix + "... done (" + minutes + "min)", 4)
        elif counterperc < 100:
            pd2tools.log(prefix + "... " + str(counterperc) + "%", 4)
    

class ScorePacket:
    """Container, collects data into one object for multiprocessing."""
    
    def __init__(self, packetid, config, typeA, typeB, keysA, keysB):
        """Create an object that defines a scoring calculation set.
                
        packetid -- a code identifying this packet
        config -- configuration for phenodigm2        
        typeA - disease/model, will convert to "HP" or "MP"
        typeB - disease/model, will convert to "HP" or "MP"
        keysA - (a subset of) the keys in tableA to process in packet.
        keysB - (a subset of) the keys in tableB to process in packet.
        """
                
        # plain copy of the input parameters                
        self.packetid = packetid
        self.config = config        
        type2onto = {"model": "MP", "disease": "HP"}
        self.ontologyA = type2onto[typeA]
        self.ontologyB = type2onto[typeB]
        self.tableA = typeA+"_phenotype"
        self.tableB = typeB+"_phenotype"            
        self.keysA = keysA
        self.keysB = keysB
        # determine cache file for this packet
        rootdir, datadir, resourcesdir, dbdir = pd2tools.getPD2dirs(config)
        self.cachebase = "scores-"+typeA+"-"+typeB+"-"+packetid+".txt.gz"
        self.cachefile = os.path.join(dbdir, self.cachebase)
        self.musthaves = pd2dss.MapSets()                

    def setMustHaves(self, musthaves):
        """set an object musthaves of class MapSets."""
        
        self.musthaves = musthaves

    def scorePass(self, scores):
        """Determine if a set of scores quality for writing to output.
        
        score - an array with four elements"""
        
        # check based on normalized scores
        temp0 = scores[0] > self.config.phenodigm_min_perc
        temp2 = scores[2] > self.config.phenodigm_min_perc
        if temp0 or temp2:
            return True
        
        # check based on raw scores
        temp1 = scores[1] > self.config.phenodigm_min_raw
        temp3 = scores[3] > self.config.phenodigm_min_raw
        if temp1 or temp3:
            return True
        
        return False
        
    def inMustHaves(self, idA, idB):
        """check if combination of ids must be recorded."""
        return self.musthaves.haspair(idA, idB)

    def run(self):
        """Carry out work for an AssociationTable."""                        
                        
        # shorthand for this function
        conf = self.config
        
        # perhaps avoid work
        if os.path.exists(self.cachefile):
            if conf.verbose:
                pd2tools.log("Skipping: "+self.cachebase, 4)
            return
                
        pmsg = "Computing: "+self.cachebase
        if conf.verbose:
            pd2tools.log(pmsg+"... starting", 4)    
        
        start_time = time.time()
        
        # fetch background data from database
        dbfile = conf.dbfile
        oomap = fillOOMap(dbfile, self.ontologyB, self.ontologyA)              
        pAmap = fillIdPhenotypeMap(dbfile, self.tableA)    
        pBmap = fillIdPhenotypeMap(dbfile, self.tableB)    
        ps = pd2scoring.PhenodigmScoring(oomap)
        # thresholds for writing to output 
        minperc = conf.phenodigm_min_perc 
        minraw = conf.phenodigm_min_raw
        
        # get objects ids (for keysB, use the packet)
        keysA, keysB = self.keysA, self.keysB        
        keysAlen = len(keysA)
            
        # show progress only for single-core jobs
        counter, interval = 0, keysAlen*2
        if conf.cores == 1:
            interval = math.ceil(keysAlen/20)
                    
        with gzip.open(self.cachefile, 'wb') as f:        
            # calculate scores for each pair    
            for idA in keysA:                
                phenosA = pAmap.getPhenotypes(idA)
                phenosA = oomap.subsetKeys(phenosA, 2)                             
                ideal_phenos = oomap.ideal_match_phenotypes(phenosA)
                for idB in keysB:                                                     
                    phenosB = pBmap.getPhenotypes(idB)                                        
                    scores = ps.score(phenosA, phenosB, ideal_phenos)
                    if self.scorePass(scores) or self.inMustHaves(idA, idB):                          
                        phenos = ps.all_scoring_phenotypes(phenosA, phenosB)                                  
                        outs = idA+"\t"+idB+"\t" + scores2str(scores)+"\t"                    
                        outs += phenos[0]+"\t"+phenos[1]+"\n"                                    
                        f.write(outs.encode('utf-8'))                       
                                
                counter += 1
                if conf.verbose:
                    showOccasionalUpdate(counter, interval, keysAlen, pmsg)
            
            # send a final entry (signals the packet completed)
            f.write("Done\n".encode('utf-8'))
        
        end_time = time.time()
        if conf.verbose:
            showOccasionalUpdate(counter, interval, keysAlen, 
                                 pmsg, end_time-start_time)                
                

def scores2str(scores):
    """Formats scores into a long tab-separated string."""
    
    result = [str(round(_, 2)) for _ in scores]
    return "\t".join(result)
                    

def getPartSize(size, cores):
    """Implements a 'magic formula' to split a whole of given size
    into parts as a function of processing cores. 
    
    Returns the number of items (<= size) that should be 
    assigned to each part."""
    
    cores = math.ceil(cores)
    if cores <= 1:    
        return size        
    result = math.ceil(size/(6*cores))    
    return min(size, result)


class AssociationScoring:
    """Object for computing phenodigm scores for models, diseases.
    
    Uses disk caching.
    
    Can be used to populate a db table, or produce single output.    
    """

    # Set keysA, keysB to specific model/disease arrays
    # to limit calculation, or leave None to compute for all possible
    keysA = None
    keysB = None
    # set dbsave to true to store computations in db
    # set to false to output calculations to screen
    dbsave = True    
    # determine whether output p2tools.log messages
    verbose = True
    # set a nontrivial packetid to dump all output into one cache file 
    packetid = None

    def __init__(self, config):
        """Init by remembering the config dbtable."""
        self.config = config        
        
        # collect information about genes in diseases and models
        mgmodel = pd2models.ModelModelGenotype(config.dbfile)
        modelgenes = pd2extractors.getDbMapSets(mgmodel)         
        diseasegenes = pd2extractors.getAllDiseaseGenes(config.dbfile)
        # look for disease/models that have genes in common
        # self.musthaves = pd2dss.MapSets()        
        # for onedisease in diseasegenes.keys():
        #     dgenes = diseasegenes.getset(onedisease)
                        
        #     for onemodel in modelgenes.keys():
        #         mgenes = modelgenes.getset(onemodel)
        #         if not mgenes.isdisjoint(dgenes):
        #             # the model and the disease have genes in common
        #             self.musthaves.add(onedisease, onemodel)
        #             self.musthaves.add(onemodel, onedisease)
        self.musthaves = pd2dss.OptimisedMapSets()
        for onedisease in diseasegenes.keys():
            for onemodel in modelgenes.keys():
                if self.musthaves.quick_intersect(onedisease, onemodel):
                    # the model and the disease have genes in common
                    self.musthaves.add(onedisease, onemodel)
                    self.musthaves.add(onemodel, onedisease)                 

    def scoreAssociations(self, typeA="model", typeB="model"):
        """Compute association scores
        
        The arguments determine what ontologies and phenotype collections
        are used during the comparison. This single function is suitable
        for HP-HP, HP-MP, and MP-MP comparisons.
        
        The variable self.dbsave determines whether output goes to db or
        is printed.
        
        Variables self.keysA and self.keysB can focus the calcualtions
        only on a subset of the associations. 
        """
        
        self.config.verbose = self.verbose
                    
        tableA, tableB = typeA + "_phenotype", typeB + "_phenotype"
        tabnameAB = typeA + "_" + typeB + "_association"
        if self.verbose:   
            pd2tools.log(tabnameAB, 2)
                      
        # read phenotypes for AB collections
        dbfile = self.config.dbfile    
        pAmap = fillIdPhenotypeMap(dbfile, tableA)    
        pBmap = fillIdPhenotypeMap(dbfile, tableB)    
        
        # lists for looping over collections
        keysA, keysB = self.keysA, self.keysB
        if keysA is None:
            keysA = sorted(list(pAmap.map.keys()))
        if keysB is None:
            keysB = sorted(list(pBmap.map.keys()))
        
        keysAlen, keysBlen = len(keysA), len(keysB)  
        if self.verbose:              
            pd2tools.log("Size: " + str(keysAlen) + " x " + str(keysBlen), 4)
            
        # for transferring results into database
        mda_data = pd2models.ModelAssociation(dbfile, tabnameAB)  
        mda_data.emptyTable()    
                        
        # split keys into parts for multiprocessing         
        partsize = getPartSize(keysAlen, self.config.cores)
        numparts = math.ceil(keysAlen/partsize)    
            
        # split the calculation into packets        
        packets = []
        for x in range(0, keysAlen, partsize):                            
            xkeysA = keysA[x:x+partsize]
            xpart = math.ceil(x/partsize)+1
            xid = self.packetid   
            if self.packetid is None:      
                xid = str(xpart) + "-" + str(numparts)
            xpack = ScorePacket(xid, self.config, typeA, typeB, 
                                xkeysA, keysB)
            xpack.setMustHaves(self.musthaves)
            packets.append(xpack)   
                    
        # calculate scores for each packet
        with mp.Pool(self.config.cores) as pool:    
            pool.map(runProcess, packets)
                    
        # transfer data into db
        if self.verbose:        
            pd2tools.log("Inserting into db", 4)
        
        # write out a header on screen
        if not self.dbsave:            
            print("\t".join(mda_data.fieldnames))
        
        # write out content of packets
        for xpack in packets:  
            with gzip.open(xpack.cachefile, 'rb') as f:
                for line in f:
                    fields = line.decode("utf-8")[:-1].split("\t") 
                    if fields[0] == "Done":
                        continue                                            
                    if not self.dbsave:
                        print(line.decode("utf-8")[:-1])
                        continue
                    if fields[0] == fields[1]:
                        continue
                    mda_data.addDataArray(fields)
            if self.dbsave:                                     
                mda_data.save()
            
        # print a db summary    
        if self.verbose and self.dbsave:
            numrows = mda_data.countrows()    
            pd2tools.log("Inserted "+str(numrows)+" rows", 4)
        

# ############################################################################
# run this from outside the module

def runScoring(config):
    """Compute phenodigm scores; insert into *_association tables."""
                      
    pd2tools.log("Computing phenodigm scores")            
    pd2tools.log("Using cores: "+str(config.cores))

    # set the types of comparisons to compute         
    ABpairs = [("disease", "model"), ("disease", "disease"), 
               ("model", "disease"), ("model", "model")]
    if config.fast:
        ABpairs = [("disease", "model")]
    
    mdas = AssociationScoring(config)                   
    for A, B in ABpairs:
        mdas.scoreAssociations(A, B)

