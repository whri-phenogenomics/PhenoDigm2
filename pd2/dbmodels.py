"""Classes for inserting data into database tables.

All classes start with prefix 'Model' and extend PhenodigmTable.
They handle insertion of data into different types of db tables.

@author: Tomasz Konopka
"""

import sqlite3


# ############################################################################

class PhenodigmTable:
    """A base class to hold data for a Phenodigm db table.
    
    This class opens/closes a connection quickly at each operation.
    This should make it possible to use this class in a 
    multiprocessing environment.
    """
    
    # insertN determines the number of rows that are sent at once
    # too high number can cause trouble with some sqlite3 installs
    insertN = 480

    # name of table to insert into
    tabname = "set this in the child class"
        
    # field names for the table in an array
    fieldnames = ["set", "this", "in", "child"]
        
    def __init__(self, dbfile, tabname=None):
        """Set up a connection and cursor for db operations.
        
        dbfile -- filename for the database
        tabname -- name of table (overrides name set by class)
        """
        
        self.data = []
        self.dbfile = dbfile
        
        if tabname is not None:
            self.tabname = tabname

    def getConn(self):
        """Create a connection that returns rows with associative arrays."""
        
        conn = sqlite3.connect(self.dbfile, timeout=1800)            
        conn.row_factory = sqlite3.Row            
        return conn

    def countrows(self):
        """Return the number of rows in a table.
        
        Note this is the number of rows saved in the db,
        not the rows in the self.data object.
        """
        
        f0 = self.fieldnames[0]
        sql = "SELECT COUNT(" + f0 + ") FROM " + self.tabname
        conn = self.getConn()
        with conn:
            result = conn.cursor().execute(sql).fetchone()[0]        
        conn.close()
        return result

    def save(self):
        """Send contents of self.data to the database table."""
        
        # prepare an insert statement ready for binding
        sql1 = "INSERT INTO " + self.tabname + " "        
        sql2 = "( " + ", ".join(self.fieldnames)+" )"
        questionmarks = ["?" for _ in self.fieldnames]
        sql3 = "("+ ", ".join(questionmarks)+" )"
        sql = sql1 + sql2 + " VALUES " + sql3                

        conn = self.getConn()
        with conn:
            c = conn.cursor()
            # execute the insert in batches of insertN
            for x in range(0, len(self.data), self.insertN):                            
                xdata = self.data[x:x+self.insertN]
                c.executemany(sql, xdata)
        conn.close()                                                                                
        self.clear()

    def clear(self):
        """Remove everything from the current data store."""
        self.data = []

    def emptyTable(self):
        """Remove everything from the DB table."""
        
        sql = "DELETE FROM "+self.tabname
        conn = self.getConn()
        with conn:
            conn.cursor().execute(sql)
        conn.close()

    def getall(self):
        """Retrieve all the data from this table."""
        
        herefields = ", ".join(self.fieldnames)
        sql = "SELECT "+herefields+" FROM "+self.tabname

    def __repr__(self):
        return f"{self.data}"
        

class PhenodigmSimpleGenerator:
    """An class that retrieves rows from a phenodigm table."""
    
    def __init__(self, phenodigmtable):
        self.table = phenodigmtable

    def next(self):
        """Retrieve all the data from the table, one at a time."""
                
        # create select statement for the table
        fields = ", ".join(self.table.fieldnames)
        sql = "SELECT "+fields+" FROM "+self.table.tabname
        
        # execute query and yield one row at a time
        with self.table.getConn() as conn:
            cur = conn.cursor()
            cur.execute(sql)                
            for row in cur:                
                yield row  
                                  
        
class PhenodigmJoinGenerator():
    """Class that performs a select join on two tables."""
    
    def __init__(self, table1, table2, joinon):
        """Initialize the generator with two table models to join.
        
        table1, table2 - two PhenodigmTable objects
        joinon - tuble with two elements, fields to join on
        """        
        
        self.table1 = table1
        self.table2 = table2
        self.joinon = joinon

    def fieldsql(self, fieldnames, prefix):
        """create sql-like snippets to select fields from a table
        
        e.g. ["a", "b"] and prefix "x" -> ["x.a as x_a", "x.b as x_b"]
        """
                    
        r = [prefix+"."+_ + " AS "+prefix+"_"+_ for _ in fieldnames]            
        return r

    def next(self):
        """Perform a select and retrieve all rows."""
                
        # get table names                
        tab1 = self.table1.tabname
        f1 = self.fieldsql(self.table1.fieldnames, self.table1.tabname)                
        f1all = ", ".join(f1)
        
        tab2 = self.table2.tabname
        f2 = self.fieldsql(self.table2.fieldnames, self.table2.tabname)
        f2all = ", ".join(f2)
                        
        sqlselect = "SELECT "+f1all+", "+f2all +" FROM "
        sqlon = tab1+"."+self.joinon[0]+"="+tab2+"."+self.joinon[1]
        sql = sqlselect + tab1 + " JOIN " + tab2 + " ON "+ sqlon
                
        # execute the query and yield one row at a time
        with self.table1.getConn() as conn:
            cur = conn.cursor()
            cur.execute(sql)                
            for row in cur:                
                yield row


# ############################################################################

class ModelGene(PhenodigmTable):
    """Model for inserting data into the gene table."""
    
    tabname = "gene"
    fieldnames = ["id", "organism", "symbol", "name", "altname", 
                  "type", "locus", "withdrawn"]
        
    def addData(self, gene_id=None, organism=None, symbol=None,
                name=None, altname=None, 
                type=None, locus=None, withdrawn=1):
        """Record one link between an id and a symbol"""
        
        if gene_id == "" or symbol == "":
            return
        
        self.data.append([gene_id,  organism, symbol, name, altname, 
                          type, locus, withdrawn])
    

class ModelGeneGeneMapping(PhenodigmTable):
    """Model for inserting data into gene_ortholog."""
    
    tabname = "gene_gene_mapping"
    fieldnames = ["query", "match"]
    
    def addData(self, query=None, match=None):
        """Record one association between non-human and human genes."""
        
        if query == "" or match == "":
            return
        
        self.data.append([query, match])


class ModelOntology(PhenodigmTable):
    """Model for inserting data into an ontology definition table."""
    
    tabname = "ontology"
    fieldnames = ["id", "term"]
    
    def addData(self, id=None, term=None):
        self.data.append([id, term.replace("\"", "'")])        
        
        
class ModelOntologySynonym(PhenodigmTable):
    """Model for inserting data into an ontology synonyms table."""
    
    tabname = "ontology_synonym"
    fieldnames = ["id", "synonym"]
    
    def addData(self, id=None, synonym=None):
        if isinstance(synonym, list) or isinstance(synonym, set):
            for s in synonym:
                self.data.append([id, s.replace("\"", "'")])
        else:    
            self.data.append([id, synonym.replace("\"", "'")])        

    
class ModelDisease(PhenodigmTable):
    """Model for inserting data into the disease table."""
    
    tabname = "disease"
    fieldnames = ["id", "term", "alts", "class"]
    
    def addData(self, id=None, term=None, alts=None, cname=None):
        self.data.append([id, term, alts, cname])
        

class ModelDiseaseGeneMapping(PhenodigmTable):
    """Model for storing associations between a disease and a gene."""
    
    tabname = "disease_gene_mapping"
    fieldnames = ["query", "match", "locus", "source"]

    def addData(self, query=None, match=None, locus=None, source=None):
        
        # avoid inserting empty mappings
        if query=="" or match=="":
            return
        
        self.data.append([query, match, locus, source])


class ModelDiseaseUpdate(PhenodigmTable):
    """An unusual model in this series. Checks if a given 
    id is present in the disease table, and adds if not present."""

    # Is this a problem? Does the behaviour change unfavourably if fixed?
    fielnames = ["id", "term"]
    checksql = "SELECT COUNT(id) FROM disease where id=?"
    insertsql = "INSERT INTO disease (id, term) VALUES (?, ?)"
        
    def addData(self, id=None, term=None):
        """Check if id exists in table, and perhaps insert."""
                        
        with self.getConn() as conn:
            c = conn.cursor()
            result = c.execute(self.checksql, [id]).fetchone()[0]
            if result == 0:
                c.execute(self.insertsql, [id, term])


class ModelModel(PhenodigmTable):
    """Model for inserting into table with model details, 
    i.e. table with allele, zygosity, source, genetic_background"""

    tabname = "model"
    fieldnames = ["id", "source", "species",
                  "genetic_background", "life_stage", "description"]
        
    def addData(self, id=None, source=None, species=None,
                genetic_background=None, life_stage=None,
                description=None):
        self.data.append([id, source, species,
                          genetic_background, life_stage, description])


class ModelModelGenotype(PhenodigmTable):
    """For insert into table with model genotype,
    i.e. description of modified genes within a model with a model id"""
    
    tabname = "model_genotype"
    fieldnames = ["id", "gene_id", "description"]
    
    def addData(self, id=None, gene=None, description=None):
        if isinstance(gene, list):
            for g in gene:
                self.data.append([id, g, description])
        else:
            self.data.append([id, gene, description])


class ModelIdPhenotype(PhenodigmTable):
    """For insert into disease_hp, associations of disease to phenotype."""
    
    tabname = "disease_phenotype"
    fieldnames = ["id", "phenotype"]
    
    def addData(self, id=None, phenotype=None):
        if isinstance(phenotype, list):
            for p in phenotype:
                self.data.append([id, p])
        else:
            self.data.append([id, phenotype])
                    

class ModelModelPhenotype(ModelIdPhenotype):
    """A special case of ModelIdPhenotype for mouse models."""
    
    tabname = "model_phenotype"
    
         
class ModelDiseasePhenotype(ModelIdPhenotype):
    """A special case of ModelIdPhenotype for diseases."""
    
    tabname = "disease_phenotype"         
         

class ModelOntologyOntologyMapping(PhenodigmTable):
    """For insert into ontology-ontology term comparison tables."""    
    
    tabname = "ontology_ontology_mapping"
    fieldnames = ["query", "match", "simJ", "ic", "lcs"]
    
    def addData(self, query, match, simJ, ic, lcs):
        """For this model, accept insertions of particular type only.
        
        Accept insertions of type A, B, data.
        But do not acccept of type B, A, data.
        """
        
        if match < query:
            return            
        self.data.append([query, match, simJ, ic, lcs])
        if match != query:
            self.data.append([match, query, simJ, ic, lcs])
                        
        
class ModelAssociation(PhenodigmTable):
    """For insert into model_disease_association."""

    tabname = "model_disease_association"
    fieldnames = ["query", "match", 
                  "score_avg_norm", "score_avg_raw", 
                  "score_max_norm", "score_max_raw",                 
                  "query_phenotype", "match_phenotype"]

    def addData(self, q, m, 
                s_avg_norm, s_avg_raw, s_max_norm, s_max_raw, 
                p1, p2):
        self.data.append([q, m, 
                          s_avg_norm, s_avg_raw, s_max_norm, s_max_raw, 
                          p1, p2])

    def addDataArray(self, arr):
        self.data.append(arr)

