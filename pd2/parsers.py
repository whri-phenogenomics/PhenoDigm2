""" Classes for parsing particular files, e.g. OBO, Omim, Orphanet.

These classes read data from disk and store it in objects in memory.
(Any db interaction with phenodigmDB should be carried out
elsewhere)

@author: Tomasz Konopka
"""

import csv
import gzip
from collections import OrderedDict
import xml.etree.ElementTree as XML

from . import tools as pd2tools
from . import dss as pd2dss

readHeader = pd2tools.readHeader


# ############################################################################
# Parsing of gene definitions from MGI and HGNC 

class GeneParser:
    """Parser for a gene definition file."""
       
    def __init__(self):        
        # genes will be a dict linking an id to a Gene object
        self.genes = dict()
                
    def parseMGI(self, datafile, headerfile):
        """Parse gene defs from MGI_EntrezGene.rpt and its header."""
   
        # get indexes in table that correspond to columns of interest
        header = readHeader(headerfile, ",")       
        MGIid = header.index("MGI_Marker_Accession_ID")
        MGIsymbol = header.index("Marker_Symbol")
        MGIstatus = header.index("Status")
        MGIname = header.index("Marker_Name")
        MGIsynonyms = header.index("Synonyms")
        MGItype = header.index("Type")
        MGIchr = header.index("Chromosome")
        MGIstart = header.index("Genome_Coordinate_Start")
        MGIend = header.index("Genome_Coordinate_End")
                
        # read from files into a dictionary        
        with gzip.open(datafile, "rt") as f:
            reader = csv.reader(f, delimiter="\t", quotechar="\"")
            for fields in reader:                    
                mgi = fields[MGIid]
                if mgi not in self.genes:
                    self.genes[mgi] = pd2dss.GeneInfo(mgi, "mouse")
                if fields[MGIstatus] == "W":
                    self.genes[mgi].addSymbol(fields[MGIsymbol], True)
                else:
                    self.genes[mgi].addSymbol(fields[MGIsymbol], False)
                    locus = fields[MGIstart]+"-"+fields[MGIend]
                    self.genes[mgi].locus = fields[MGIchr]+":"+locus
                    self.genes[mgi].type = fields[MGItype]
                    self.genes[mgi].name = fields[MGIname]
                    self.genes[mgi].altname = fields[MGIsynonyms]

    def parseHGNCWithdrawn(self, datafile, headerfile):
        """Parse withdrawn gene defs from hgnc."""
                    
        header = readHeader(headerfile, ",")
        Hid = header.index("HGNC_ID")
        Hsymbol = header.index("WITHDRAWN_SYMBOL")
        Hstatus = header.index("STATUS")
        Hname = header.index("MERGED_INTO_REPORTS")
        
        # read from files into a dictionary        
        with gzip.open(datafile, "rt") as f:
            reader = csv.reader(f, delimiter="\t", quotechar="\"")
            for fields in reader:                    
                hgnc = fields[Hid]        
                # skip the header
                if hgnc == "HGNC_ID":
                    continue
                
                if hgnc not in self.genes:                    
                    self.genes[hgnc] = pd2dss.GeneInfo(hgnc, "human")
       
                self.genes[hgnc].addSymbol(fields[Hsymbol], True) 
                if fields[Hname] != "":
                    self.genes[hgnc].altname = fields[Hname]

    def parseHGNC(self, datafile, headerfile):                        
        """Parse gene defs from hgnc_complete_set.txt and header."""
        
        header = readHeader(headerfile, ",")
        Hid = header.index("hgnc_id")
        Hsymbol = header.index("symbol")
        Hname = header.index("name")
        Haltname = header.index("prev_name")
        Hlocus = header.index("location")
        Htype = header.index("locus_type")
        Halias = header.index("alias_symbol")
        Hstatus = header.index("status")
        Hprev = header.index("prev_symbol")
        
        # read from files into a dictionary        
        with gzip.open(datafile, "rt") as f:
            reader = csv.reader(f, delimiter="\t", quotechar="\"")
            for fields in reader:                    
                hgnc = fields[Hid]
                # skip the header
                if hgnc == "hgnc_id":
                    continue
                
                if hgnc not in self.genes:
                    self.genes[hgnc] = pd2dss.GeneInfo(hgnc, "human")
                if fields[Hstatus] != "Approved":
                    self.genes[hgnc].addSymbol(fields[Hsymbol], True)
                else:
                    self.genes[hgnc].addSymbol(fields[Hsymbol], False)                    
                    self.genes[hgnc].locus = fields[Hlocus]
                    self.genes[hgnc].type = fields[Htype]
                    self.genes[hgnc].name = fields[Hname]
                    self.genes[hgnc].altname = fields[Haltname]                    
                    for onealias in fields[Halias].split("|"):
                        self.genes[hgnc].addSymbol(onealias, True)
                    for oneprev in fields[Hprev].split("|"):
                        self.genes[hgnc].addSymbol(oneprev, True)


# ############################################################################
# Parsing of obo (ontology files)

class OboSynonymsParser:
    """Bare-bones parser to extract id-term-synonym mappings
    from an obo file."""
    
    def __init__(self, filepath):        
        self.terms = OrderedDict()
        self.synonyms = OrderedDict()
        self.parse(filepath)            
    
    def parse(self, filename):
        """Scan an obo file and looks for 'id': and 'name:' rows"""
        
        nowid = ""
        nowname = ""                
        with open(filename, "r") as f:
            for line in f:
                line2 = line[:-1].split(": ", 1)                
                # if line starts with id:, this is a start of [Term]                
                if line2[0] == "id":
                    nowid = line2[1]                
                # ignore most lines, but record name: and synonym:     
                if line2[0] == "name":
                    nowname = line2[1]
                    self.terms[nowid] = nowname                    
                if line2[0] == "synonym":
                    nowsyn = line2[1]
                    if nowid not in self.synonyms:
                        self.synonyms[nowid] = set()                    
                    # extract the synonym text from between quotes
                    secondq = nowsyn.find('"', 2)
                    if secondq > 0:
                        nowsyn = nowsyn[1:secondq]
                    if nowsyn != nowname:
                        self.synonyms[nowid].add(nowsyn)


# ############################################################################
# Disease classifications from orphanet xmls

class OrphanetClassificationParser:
    """Parse information from orphanet linking disease ids to classes."""
    
    def __init__(self):
        """Create a dict that will connect OrphaNumbers to disease classes."""
        self.classes = dict()

    def parseFile(self, filepath, classname):
        """Read one file, transfer orphanumbers in the file into 
        this instance's dict."""        
        
        # load entire file contents
        filedata = XML.parse(filepath)        
        
        # scan all the tags, identify orphanumbers
        for xmlnode in filedata.iter():
            if xmlnode.tag == "OrphaCode":
                xmlnumber = "ORPHA:" + xmlnode.text
                if xmlnumber not in self.classes:
                    self.classes[xmlnumber] = set()                    
                self.classes[xmlnumber].add(classname)                
            

# ############################################################################
# Disease definitions from orphanet (ids, titles, synonyms)

class OrphanetDiseaseParser:
    """Parse Orphanet xml with disease ids, titles, synonyms"""
    
    # the orphanumber that means a match is exact
    orpha_exact = 377807
    # some disease relations are non-exact
    # (narrower to broader, or broader to narrower)
    orpha_ntbt = 377808
    orpha_btnt = 377809
    
    def __init__(self, filepath=None):
        """Create mappings from OrphaNumbers to disease definitions."""
        
        self.synonyms = dict()
        self.titles = dict()
        self.exact = dict()
        
        if filepath is not None:
            self.parse(filepath)

    def parseOneExternal(self, xmlnode):
        """Parse content of a <ExternalReference> tag.        
        Returns a 2-tuple.
        First item is synonym id.
        Second item is a logical value, True if match is exact, False if not
        """
                
        source = xmlnode.find("Source").text
        ref = xmlnode.find("Reference").text            
        relation = xmlnode.find("DisorderMappingRelation")
        name = relation.find("Name").text
        result = source + ":" + ref
        if "Exact" in name or "equivalent" in name:
            return result, True
        if "BTNT" in name or "NTBT" in name:
            return result, False
        return None, True
                    
    def parseOneDisorder(self, xmlnode):
        """Parse content of one <Disorder> tag."""

        orphanumber = xmlnode.find("OrphaCode").text
        orphaid = "ORPHA:" + orphanumber
        title = xmlnode.find("Name").text
        self.titles[orphaid] = title.title()

        external = xmlnode.find("ExternalReferenceList")
        for oneexternal in external:
            onematch = self.parseOneExternal(oneexternal)
            oneid = onematch[0]
            isexact = onematch[1]
            if oneid is not None:
                # perhaps skip if an exact match already found
                if oneid in self.exact and self.exact[oneid]:
                    continue
                # otherwise record this synonym
                self.synonyms[oneid] = orphaid
                self.exact[oneid] = isexact

    def parse(self, filepath):
        """Read one file, transfer into synonyms dict."""        
        
        # load entire file contents
        filedata = XML.parse(filepath)        
        fileroot = filedata.getroot()
        
        # expect a certain structure
        # JDBOR -> DisorderList -> many <Disorder>
        for dl in fileroot:
            for disorder in dl:
                try:
                    self.parseOneDisorder(disorder)
                except Exception:
                    pass
                
        
# ############################################################################
# Loading of disease gene associations from Orphanet

def parse_loci(xmlnode):
    """Parse content of a <Gene> tag. returns a string summarizing locus."""
    
    result = ""
    for locus in xmlnode.find("LocusList"):
        if result != "":
            result += ", "
        result += locus.find("GeneLocus").text            
    return result


def parse_genelist_loci(xmlnode):
    """create a dictionary linking internal Orphanet gene ids to loci
    
    Note: The GeneList node is now removed from the orphanet download file.
    Hence this function is no longer used.
    """
    
    geneloci = dict()
    genelist = xmlnode.find("GeneList")        
    for onegene in genelist:
        tempid = onegene.attrib["id"]
        # in some version of the file, parsing the locus doesn't work
        onelocus = parse_loci(onegene)
        geneloci[tempid] = pd2dss.WithLocus("", onelocus)     
    return geneloci


class OrphanetGeneParser:
    """Parse Orphanet xml with mappings from disease to loci, genes"""

    def __init__(self, filepath=None):
        """Create mappings from OrphaNumbers to disease loci and genes."""
                
        self.genelocus = pd2dss.MapSets()                
        if filepath is not None:
            self.parse(filepath)

    def parseOneDisorderAssociation(self, xmlnode):
        """Parse content of a <DisorderGeneAssociation> tag.
        Returns an xml id and string with HGNC id."""
        
        result = ""        
        gene = xmlnode.find("Gene")                
        for oneexternal in gene.find("ExternalReferenceList"):
            onesource = oneexternal.find("Source").text            
            if onesource == "HGNC":                
                result = "HGNC:"+oneexternal.find("Reference").text
        return gene.attrib["id"], result

    def parseOneDisorder(self, xmlnode):
        """Parse the orphanet xml for one <Disorder> tag."""

        orphanumber = xmlnode.find("OrphaCode").text
        orphaid = "ORPHA:"+orphanumber                                
            
        associationlist = xmlnode.find("DisorderGeneAssociationList")
        for oneassoc in associationlist:
            tempid, geneid = self.parseOneDisorderAssociation(oneassoc)                                    
            genelocus = pd2dss.WithLocus(geneid, "")            
            self.genelocus.add(orphaid, genelocus)

    def parse(self, filepath):
        """Read one file, transfer into synonyms dict."""        
        
        # load entire file contents
        filedata = XML.parse(filepath)        
        fileroot = filedata.getroot()
        
        # expect a certain structure
        # JDBOR -> DisorderList -> many <Disorder>
        for dl in fileroot:
            for disorder in dl:
                try:
                    self.parseOneDisorder(disorder)
                except Exception as e:
                    print("OrphanetGene parser: skipping " + str(disorder))
                    pass        


# ##################################################################
# Loading of disease data (OMIM)

class MorbidEntryParser:
    """A Parser for one entry in the morbidMap"""
    
    def __init__(self, header):
        # definitions for format of one morbid entry        
        self.mimindex = header.index("MimNumber")        
        self.cytoindex = header.index("CytoLocation")
        self.phenoindex = header.index("Phenotype")
        self.geneindex = header.index("GeneSymbols")
        # result of parsing
        # mimnumber is sometimes the disease mimumber from first column
        # but when that is unavailable, it is the locus mimnumber
        self.mimnumber = ""
        # lcsnumber is a mimnumber that comes from the MimNumber column 
        self.lcsnumber = ""
        # genes are the symbols associated with the locus
        self.genes = []
        # locus is a string with genomic location
        self.locus = ""
        # category and title are components that make up the phenotype desc
        self.category = 0
        self.title = ""       
       
    def parse(self, line):
        """Parse a line (string).
        
        This updates values for mimnumber, genes, locus."""
        
        fields = line[:-1].split("\t")
                
        # get the naive field values
        self.mimnumber = self.lcsnumber = "OMIM:"+fields[self.mimindex]
        self.genes = fields[self.geneindex].split(", ")
        self.locus = fields[self.cytoindex]        
        self.category = 0
        
        # parse the phenotype field, i.e. first column
        phendata = fields[self.phenoindex].split(" ")
        phenlen = len(phendata)
        # the category is always a number at the end, i.e. (X)
        self.category = phendata[phenlen-1][1:2]
        # the disease number of sometimes just before the category
        penultimate = phendata[phenlen-2]
        if penultimate.isdigit() and len(penultimate) > 3:
            self.mimnumber = "OMIM:"+penultimate
            self.title = " ".join(phendata[:-2])
        else:
            self.title = " ".join(phendata[:-1])
                

class OmimParser:
    """Parse disease data from omim mimTitles and morbidmap."""
    
    def __init__(self):        
        # mimdata will hold a map between disease id and array(title, alt)
        self.mimdata = dict()
        # mimdefined holds what omim ids are defined 
        # (this is relevant because parser skips some lines in mimTitles)
        self.mimdefined = set()
        # mimgenes will hold a map from mim numbers to HGNC symbols
        self.mimgenes = dict()
        # symbolsloci is a dict mapping a disease id to gene symbols+locus        
        self.symbolsloci = dict()

    def parseMim2gene(self, datafile, headerfile):
        """Parse file mim2gene.txt and its header, record mim hgnc links."""
        
        mheader = readHeader(headerfile, "\t")
        mimindex = mheader.index("MimNumber")
        typeindex = mheader.index("EntryType")
        hgncindex = mheader.index("HGNC")
        
        # transfer mimnumbers and hgnc symbols into dict                            
        with gzip.open(datafile, "rt") as f:
            for line in f:
                # ignore comment lines, gene definitions, deprecated
                if line.startswith("#"):
                    continue
                
                fields = line[:-1].split("\t")
                nowmim = "OMIM:"+fields[mimindex]
                nowtype = fields[typeindex]
                
                if nowtype == "gene" or nowtype == "gene/phenotype":                                
                    self.mimgenes[nowmim] = fields[hgncindex]

    def parseTitles(self, datafile, headerfile):
        """Parse file mimTitles.txt and its header."""
        
        mheader = readHeader(headerfile, "\t")
        prefixindex = mheader.index("Prefix")        
        mimindex = mheader.index("MimNumber")
        titleindex = mheader.index("Title")
        altindex = mheader.index("AltTitle")
        incindex = mheader.index("IncTitle")

        # transfer titles and alt titles from file into dict                            
        with gzip.open(datafile, "rt") as f:
            for line in f:
                # ignore comment lines, gene definitions, deprecated
                if line.startswith("#"):
                    continue
                if line.startswith("Asterisk") or line.startswith("Plus"):
                    continue
                if line.startswith("Caret"):
                    continue
                              
                fields = line[:-1].split("\t")
                nowmim = "OMIM:"+fields[mimindex]
                nowtitle, nowalt = "", ""
                if len(fields) > titleindex:
                    nowtitle = fields[titleindex]
                if len(fields) > altindex:
                    nowalt = fields[altindex]
                if len(fields) > incindex:
                    if len(fields[incindex]) > 0:
                        if len(nowalt) > 0:
                            nowalt += " | "
                        nowalt += fields[incindex]
                
                # record the mim-title association
                self.mimdata[nowmim] = [nowtitle, nowalt]                
                self.mimdefined.add(nowmim)

    def parseMorbid(self, datafile, headerfile):
        """Extract locus information from morbidmap.txt."""
        
        mheader = readHeader(headerfile, "\t")
        mparser = MorbidEntryParser(mheader)
        
        with gzip.open(datafile, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue                
                mparser.parse(line)                
                if mparser.mimnumber in self.mimdefined:
                    nowmim = mparser.mimnumber
                    nowlcs = mparser.lcsnumber
                    
                    if nowmim not in self.symbolsloci:
                        self.symbolsloci[nowmim] = set()
                                        
                    if nowlcs in self.mimgenes:
                        # if gene mimnumber is in mim2gene, use that
                        nowgene = self.mimgenes[nowlcs]
                        genelocus = pd2dss.WithLocus(nowgene, 
                                                     mparser.locus)
                        self.symbolsloci[nowmim].add(genelocus)
