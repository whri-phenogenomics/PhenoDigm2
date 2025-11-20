"""Containers and data structures for PhenoDigm2.

These data structures are used during the pheno-scoring calcualtions.

See also parsers.py for a collection of classes/data structure 
relevant for processing input data files.

@author: Tomasz Konopka
"""

import math

    
class PhenotypeSet:
    """Container keeps track of a set of phenotypes.
    
    A simple set only holds items. This container holds
    a set plus a string for the set title.
    """
    
    def __init__(self, title):
        self._phenotypes = set()
        self.title = title

    def add(self, hp_id):
        self._phenotypes.add(hp_id)

    @property
    def phenotypes(self):
        return sorted(list(self._phenotypes))


class IdPhenotypeMap:
    """Container holding id-phenotype associations (2-way).
    
    Example uses are for disease-phenotype and mousemodel-phenotype
    associations."""
    
    def __init__(self):
        self.map = dict()
        self.invmap = dict()            

    def add(self, id, phenotype):
        """Add one association to this set."""
        
        if id not in self.map:
            self.map[id] = set()
        self.map[id].add(phenotype)
        
        if phenotype not in self.invmap:
            self.invmap[phenotype] = set()
        self.invmap[phenotype].add(id)                    

    def getPhenotypes(self, id):
        """Get a set of phenotypes associated with an id."""
        
        if id in self.map:
            return self.map[id]
        else:
            return set()        

    def getIds(self, phenotype):
        return self.invmap[phenotype]

    def __str__(self):
        """Create a representation of dict of sets."""
        
        result = []
        keys = sorted(list(self.map.keys()))       
        for id in keys:            
            temp = [_ for _ in self.map[id]]
            result.append(id+": {"+ ", ".join(sorted(temp))+"}")
        return "\n".join(result)

# NOTE: This might not be needed since cs becomes the score.
def ooscore(simJ, ic):
    """get a score for an ontology/ontology mapping."""
    return math.sqrt(simJ*ic)

        
class OOMapData:
    """A container holding scores associated with a
    ontology-ontology mapping."""
    
    def __init__(self, id, cs):
        self.id = id
        self.cs = cs
        self.score = cs

    def checkUpdate(self, obj):
        """Update the values of simJ, ic, lcs if new
        set of values increases the score."""
                        
        if obj.score > self.score:
            self.score = obj.score
            self.id = obj.id

    def __str__(self):
        # result = self.id+ ": ["+str(self.simJ)+", "+str(self.ic)
        # result += ", "+str(self.lcs)+", "+str(self.score)+"]"

        result = self.id + ": ["+str(self.score)+"]"
        return result
            

class OOInvData:
    """Simple container with score and label;
    for holding inverse mapping."""

    def __init__(self, score, id):
        self.score = score
        self.id = id

    def __str__(self):
        return "(" + str(self.score) + ", " + str(self.id) + ")"
        

class OntologyOntologyMap:
    """Container holding many ontology-ontology mappings.
    
    Onto-onto mappings are stored in a depth-2 dict with 
    ontology terms as keys.
    
    Best reverse mappings are stored in a depth-1 dict.     
    """
    
    def __init__(self):
        """Initialize an empty mapping."""
        
        self.map = dict()
        self.bestmapped = dict()
        self.keys1 = set()
        self.keys2 = set()
        
    def add(self, row):
        """Store a set of mappings into a 2D dict.
        
        Arguments:
            rows -- an object with keys query, hit, simJ, ic, lcs.
        """
        
        id1 = row["query"]
        id2 = row["match"]
        
        self.keys1.add(id1)
        self.keys2.add(id2)
        
        if id1 not in self.map:
            self.map[id1] = dict()
            
        # rowentry = OOMapData(id1, row["simJ"], row["ic"], row["lcs"])
        rowentry = OOMapData(id1, row["cs"])
        self.map[id1][id2] = rowentry

        if id2 not in self.bestmapped:
            self.bestmapped[id2] = OOInvData(rowentry.score, id1)              
        else:
            if rowentry.score > self.bestmapped[id2].score:
                self.bestmapped[id2] = OOInvData(rowentry.score, id1)                

    def get(self, id1, id2):
        """Get mapping data from id1 to id2."""
        
        if id1 not in self.map:            
            return None
        temp = self.map[id1]                        
        if id2 in temp:
            return temp[id2]
        else:
            return None

    def subsetKeys(self, phenolist, keylevel=1):
        """Look through a phenolist and return a subset that occurs in this oomap."""
                            
        if keylevel == 1:
            keys = self.keys1
        elif keylevel == 2:
            keys = self.keys2
        else:
            raise Exception("keylevel must be either 1 or 2")
        
        result = []
        for pheno in phenolist:
            if pheno in keys:
                result.append(pheno)
        return result
    
    def getBestmapped(self, id2):
        """Return annotation for best inverse mapping from id2."""
        
        if id2 not in self.bestmapped:
            return None
        
        return self.bestmapped[id2]
        
    def ideal_match_phenotypes(self, phenotypes):
        """Get list of disease phenotypes that would 
        be an idealized match for given list. 
        
        What hypothetical disease would fit the model perfectly?       
        """        
                    
        best_matches = set()        
        for pheno in phenotypes:
            if pheno in self.bestmapped:
                best_matches.add(self.bestmapped[pheno].id)
            
        return list(best_matches)

    def __str__(self):
        result1 = "Map:\n"
        for k in self.map.keys():            
            for k2 in self.map[k].keys():
                result1 += k + "  " + k2 + ": " + str(self.map[k][k2]) + "\n"
                
        result2 = "\nBestmapped:\n"
        for k in self.bestmapped:
            result2 += k + ": " + str(self.bestmapped[k]) + "\n"
        return result1+result2


class OrthologMapping:
    """An object that can hold many-to-many gene pairing.
    It can be constructed without full info on pairing;
    new orthologs can be added as we go along."""
    
    def __init__(self):
        self.map = dict()

    def addPair(self, itemA, itemB):
        """Add one gene to the ortholog object"""
        
        if itemA not in self.map:
            self.map[itemA] = set()
        self.map[itemA].add(itemB)
        
    def getPairs(self, itemA, itemB): 
        """Get an array of pairwise links from A to B."""
        
        result = []
        
        # avoid work if one of the organisms is not defined 
        mapkeys = self.map.keys()
        if itemA not in mapkeys or itemB not in mapkeys:
            return result
        
        # generate all pairs
        for a in self.map[itemA]:
            for b in self.map[itemB]:
                if b is not "":
                    result.append([a, b])
        return result


class MapSets:
    """a map of sets.
    
    (Perhaps this should extend dict)."""
    
    def __init__(self):
        self.map = dict()

    def keys(self):
        return self.map.keys()

    def add(self, a, b):
        """Add a connection only from a to b."""
        
        if a not in self.map:                    
            self.map[a] = set()            
        self.map[a].add(b)
        
    def get(self, a):
        """retrieve a set, or None if a is not in map."""
        
        if a not in self.map:
            return None        
        return self.map[a]

    def getset(self, a):
        """retrieve a set for key a. 
        
        In contrast to get(), this always returns a set perhaps empty."""
        
        if a not in self.map:
            return set()
        return self.map[a]

    def has(self, a):
        """retrieve a boolean for whether a is in the map."""
        
        return a in self.map

    def haspair(self, a, b):
        """retrieve a boolean for whether map has [a][b]"""
        
        if a in self.map:
            return b in self.map[a]
        return False

    def reverse(self):
        """create a new MapSets. 
        
        if self has x -> (y1, y2), etc.
        the new MapSets instance will have y1 -> x, y2 ->x, etc."""
        
        result = MapSets()        
        for a in self.map.keys():
            for b in self.map[a]:
                result.add(b, a)

        return result    


class GeneInfo:
    """A container for building up info about a gene.
    
    The class can also describe a genetic feature."""
    
    def __init__(self, id, organism, symbol=None):
        self.id = id
        self.organism = organism
        self.official_symbol = symbol
        self.type = "Gene"
        self.symbols = set()  # these are withdrawn/alt symbols
        self.name = None
        self.locus = None     
        self.altname = "-"   

    def addSymbol(self, symbol, withdrawn=False):
        """Set a symbol for this gene."""
        
        if symbol is None: 
            return
        symbol = symbol.strip()
        if symbol == "":
            return
                
        if withdrawn:
            # just add a symbol to the withdrawn set
            self.symbols.add(symbol)
        else:
            # if the new symbol is official, replace existing official
            if self.official_symbol is not None:
                self.symbols.add(self.official_symbol)
            self.official_symbol = symbol

    def isValid(self):
        """Checks if an official symbol has been assigned."""
        return self.official_symbol is not None


class WithLocus:
    """A container for a pair or strings. 
    The strings can be used for a gene symbol and a locus, 
    or for a gene id and a locus."""
    
    def __init__(self, query, locus):
        self.query = query.strip()
        self.locus = locus
        
    def __repr__(self):
        return "[" + self.query + ", " + self.locus + "]"

