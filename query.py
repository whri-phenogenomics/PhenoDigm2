"""Manual extraction of data from the db. 

Collection of functions to help manual exploration 
of data, e.g. when a user would like to extract one
phenodigm score or understand a comparison.

@author: Tomasz Konopka
"""

import hashlib
from . import tools as pd2tools
from . import score as pd2score


# ############################################################################
# Individual query functionality

def getArgArray(arrstring):
    """Convert a single string value into an array,
    e.g. xyz -> ["xyz"] but "abc xyz" -> ["abc", "xyz"]
    """
            
    result = [_.strip() for _ in arrstring.split(" ")]
    while "" in result:
        result.remove("")        
    return result


class PhenodigmQuery:
    """Holds a db conn; and interprets query from config object."""
    
    # column names for the association question
    scorecolumns = ["query", "match", "score_avg_norm", "score_avg_raw",
                    "score_max_norm", "score_max_raw",
                    "query_phenotype", "match_phenotype"]         

    def __init__(self, config):
        self.verbose = config.verbose        
        self.conn = pd2tools.getDBconn(config.dbfile)
        self.c = self.conn.cursor()
        self.models = getArgArray(config.model) 
        self.diseases = getArgArray(config.disease)
        self.genes = getArgArray(config.gene)
        self.terms = getArgArray(config.term)
        
    def print_rows(self):
        """Internal function, prints all rows got by the cursor."""
        for row in self.c.fetchall():
            result = []
            for k in row.keys():
                result.append(str(row[k]))                
            print("\t".join(result))
            
    def queryUsingSelect(self, sql, values, label, header):
        """Generic select and print function."""
        
        if self.verbose:
            pd2tools.log("Phenodigm query ["+label+"]")
        
        print(header)
        for v in values:
            self.c.execute(sql, (v,))
            self.print_rows()        
           
    def queryTerms(self):
        if len(self.terms) == 0:
            return
        sql = "SELECT id, term FROM ontology WHERE id=?"
        self.queryUsingSelect(sql, self.terms, "term", "id\tterm")                    
    
    def queryModels(self):
        if len(self.models) == 0:
            return
        cols = ["id", "source", "species", 
                "genetic_background", "description"]
        sql = "SELECT "+",".join(cols)+" FROM model WHERE id=?"
        self.queryUsingSelect(sql, self.models, "model", "\t".join(cols))
        pass                    
        
    def queryGenes(self):
        if len(self.genes) == 0:
            return
        cols = ["id", "source", "species",
                     "genetic_background", "description"]
        geneslike = ["%"+g+"%" for g in self.genes]        
        sql = "SELECT "+",".join(cols) + " FROM model "
        sqlwhere = " WHERE description LIKE ?"
        self.queryUsingSelect(sql+sqlwhere, geneslike, "gene", "\t".join(cols))
        pass                    
            
    def queryDiseases(self):
        if len(self.diseases) == 0:
            return
        sql = "SELECT id, term FROM disease WHERE id=?"
        self.queryUsingSelect(sql, self.diseases, "disease", "id\tterm")
                                    
    def queryOwltoolsScores(self):
        """Display owltools scores for pairs of terms."""
        
        if self.verbose:
            pd2tools.log("Phenodigm query [owltools]")
                
        columns = ["query", "match", "simJ", "ic", "lcs"]
        sql = "SELECT " + ", ".join(columns) + " FROM "
        sql += "ontology_ontology_mapping WHERE query=? and match=?"
        
        print("\t".join(columns))
        for t1 in self.terms:
            for t2 in self.terms:
                self.c.execute(sql, (t1, t2))
                self.print_rows()
    
    def queryDiseasePhenotypes(self):
        if len(self.diseases) == 0:
            return
        sql = "SELECT id, phenotype FROM disease_phenotype WHERE id=?"
        self.queryUsingSelect(sql, self.diseases, "phenotype", "id\tphenotype")

    def queryModelPhenotypes(self):
        if len(self.models) == 0:
            return
        sql = "SELECT id, phenotype FROM model_phenotype WHERE id=?"
        self.queryUsingSelect(sql, self.models, "phenotype", "id\tphenotype")

    def queryOnePhenodigm(self, x1, x2, tabprefix):
        tabname = tabprefix+"_association"
        sql = "SELECT " + (",".join(self.scorecolumns)) + " FROM " + tabname
        sql += " WHERE query=? and match=?"
        self.c.execute(sql, (x1, x2))
        self.print_rows()
        
    def queryScore(self):
        """Delegate query phenodigm calcs for models, diseases."""
        
        if self.verbose:
            pd2tools.log("Phenodigm query [association]")
         
        models = self.models
        diseases = self.diseases
        
        print("\t".join(self.scorecolumns))         
        # query pairwise (within model-model)
        for m1 in models:
            for m2 in models:
                if m1 != m2:
                    self.queryOnePhenodigm(m1, m2, "model_model")
                
        # query pairwise (within disease-disease)
        for d1 in diseases:
            for d2 in diseases:
                if d1 != d2:
                    self.queryOnePhenodigm(d1, d2, "disease_disease")
        
        # query pairwise (within model-disease)        
        for m1 in models:
            for d1 in diseases:
                self.queryOnePhenodigm(m1, d1, "model_disease")          
        for d1 in diseases:
            for m1 in models:
                self.queryOnePhenodigm(d1, m1, "disease_model")      


class PhenodigmReCompute:
    """Holds a db conn; interprets compute request from config."""
    
    def __init__(self, config):          
        self.config = config
        self.config.cores = 1
        self.config.phenodigm_min_perc = -1
        self.config.phenodigm_min_raw = -1
        self.models = getArgArray(config.model) 
        self.diseases = getArgArray(config.disease)        
        
    def packetHashId(self, a, b):
        ss = ".".join(a) + "." + ".".join(b)
        sshash = hashlib.md5(ss.encode())
        return "query-"+sshash.hexdigest()[0:12]        
        
    def recomputeScores(self):
        """compute or re-compute the requested phenodigm scores."""
        
        # interpret the request, is it model-model, disease-disease, model-disease
        is0m, is0d = len(self.models) == 0, len(self.diseases) == 0

        if is0m and is0d:
            print("Aborting: must specify at least one model or disease")
            return

        typeA, typeB = "", ""   
        keysA, keysB = "", ""     
        if is0d:
            typeA, typeB = "model", "model"
            keysA, keysB = self.models, self.models
        elif is0m:
            typeA, typeB = "disease", "disease"
            keysA, keysB = self.diseases, self.diseases
        else: 
            typeA, typeB = "model", "disease"
            keysA, keysB = self.models, self.diseases
                
        # make a calculator object
        mdas = pd2score.AssociationScoring(self.config)
        mdas.verbose = False
        mdas.dbsave = False
        
        # calculate the forward relation
        mdas.keysA, mdas.keysB = keysA, keysB        
        mdas.packetid = self.packetHashId(keysA, keysB)        
        mdas.scoreAssociations(typeA, typeB)        
        # calculate the backward relation
        if typeA == typeB:
            return 
        print("")
        mdas.keysA, mdas.keysB = keysB, keysA        
        mdas.packetid = self.packetHashId(keysB, keysA)
        mdas.scoreAssociations(typeB, typeA)
        

# ############################################################################
# Gateway

def queryPhenodigm(config):
    """Gateway to the query function."""
    
    pe = PhenodigmQuery(config)
    
    if config.sim:
        pe.queryOwltoolsScores()
    elif config.phenotype:
        pe.queryDiseasePhenotypes()
        pe.queryModelPhenotypes()
    elif config.association or config.score:
        pe.queryScore()
    else:        
        pe.queryTerms()    
        pe.queryDiseases()    
        pe.queryModels()    
        pe.queryGenes()


def computeScores(config):
    """Gateway to re-computation of phenodigm scores."""
    
    prc = PhenodigmReCompute(config)
    prc.recomputeScores()

