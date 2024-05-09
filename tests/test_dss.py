'''
Tests for contents of pd2/scoring.py
(Classes to produce scores)
'''


import unittest
from os.path import join, dirname
from pd2 import dss


class IdPhenotypeMapTests(unittest.TestCase):
    """Test cases for class IdPhenotypeMap."""
        
    def setUp(self):
        """Setup for individual tests."""
        self.idp = dss.IdPhenotypeMap()
        self.idp.add("i0", "HP:1")
        self.idp.add("i1", "HP:4")
        self.idp.add("i0", "HP:4")
                
                                                                                                                                
    def test_avg(self): 
        idpstr = "i0: {HP:1, HP:4}\ni1: {HP:4}"             
        self.assertEqual(idpstr, str(self.idp), "string rep")
        
    def test_retrieving(self):          
        i0 = set()
        i0.add("HP:1")
        i0.add("HP:4")        
        self.assertEqual(i0, self.idp.getPhenotypes("i0"), 
                         "retrieving phenotypes")
        hp4 = set()
        hp4.add("i0")
        hp4.add("i1")
        self.assertEqual(hp4, self.idp.getIds("HP:4"), 
                         "retrieving ids")
        

class OOMapDataTests(unittest.TestCase):
    """Test cases for class OOMapData."""
        
    def setUp(self):
        """Setup for individual tests."""
        pass
                                                                                                                                                
    def test_score(self):
        ood = dss.OOMapData("abc", 1, 4, "xyz")                     
        self.assertEqual(2, ood.score, "basic score check")

    def test_str(self):
        ood = dss.OOMapData("abc", 1, 4, "xyz")
        oodstr = "abc: [1, 4, xyz, 2.0]"
        self.assertEqual(oodstr, str(ood), 
                         "simple str representation")

    def test_update(self):        
        ood = dss.OOMapData("abc", 1, 4, "xyz")
        oodstr = str(ood)
        ood2 = dss.OOMapData("abc2", 1, 1, "xyz2")
        ood3 = dss.OOMapData("abc3", 2, 8, "xyz3")        
        ood.checkUpdate(ood2)
        self.assertEqual(oodstr, str(ood), 
                         "update with lower score should be ignored")        
        ood.checkUpdate(ood3)
        self.assertEqual(str(ood3), str(ood), 
                         "update with higher score should replace all fields")


class OntologyOntologyMapTests(unittest.TestCase):
    """Test cases for class OntologyOntologyMap."""
                                
    def setUp(self):
        """Setup common map for individual tests."""
        
        # an example with maps across two different ontologies
        a1 = {"query": "hp0", "match": "mp1",
                "ic": 1, "simJ": 1, "lcs": "xyz"}
        a2 = {"query": "hp0", "match": "mp2",
                "ic": 4, "simJ": 4, "lcs": "xyz"}
        a3 = {"query": "hp1", "match": "mp1",
                "ic": 1, "simJ": 4, "lcs": "xyz"}
        self.oomap = dss.OntologyOntologyMap()
        self.oomap.add(a1)
        self.oomap.add(a2)
        self.oomap.add(a3)
        
        aaa1 = {"query": "hp0", "match": "hp0",
                "ic": 1, "simJ": 4, "lcs": "xyz"}
        aaa2 = {"query": "hp0", "match": "hp1",
                "ic": 1, "simJ": 2, "lcs": "xyz"}
        aaa3 = {"query": "hp1", "match": "hp1",
                "ic": 1, "simJ": 9, "lcs": "xyz"}                        
        self.oomap2 = dss.OntologyOntologyMap()
        self.oomap2.add(aaa1)
        self.oomap2.add(aaa2)        
        self.oomap2.add(aaa3)
                

    def test_map(self):                             
        """Check for basic number of associations to terms."""
        self.assertEqual(2, len(self.oomap.map["hp0"]), "basic length check")
        self.assertEqual(1, len(self.oomap.map["hp1"]), "basic length check")

    def test_get1(self):                             
        """Check for existence of forward associations."""                
        self.assertTrue(self.oomap.get("hp0", "mp1"), "basic has check")        
        self.assertFalse(self.oomap.get("hp1", "mp2"), "basic has check")    

    def test_inv1(self):
        """Check for bestmapped, i.e. inverse associations."""        
        mp2inv = self.oomap.bestmapped["mp2"]
        self.assertEqual("hp0", mp2inv.id, "simple inverse id")
        self.assertEqual(4, mp2inv.score, "simple inverse score")
        
    def test_inv2(self):
        """Check for bestmapped, i.e. inverse associations."""        
        mp1inv = self.oomap.bestmapped["mp1"]
        self.assertEqual("hp1", mp1inv.id, "replaced inverse id")
        self.assertEqual(2, mp1inv.score, "replaced inverse id")

    def test_ideal_match1(self):
        """Test identifying best inverse matches."""                       
        phenos = self.oomap.ideal_match_phenotypes(["mp1"])
        self.assertEqual(["hp1"], phenos, 
                         "ideal connects with high score hp terms (hp1)")        
        
    def test_ideal_match2(self):
        """Test identifying best inverse matches."""                  
        phenos = self.oomap.ideal_match_phenotypes(["mp2"])
        self.assertEqual(["hp0"], phenos, 
                         "ideal connects with only hp term (hp1)")        
        
    def test_ideal_match3(self):
        """Test identifying best inverse matches."""                      
        phenos = self.oomap.ideal_match_phenotypes(["mp99"])
        self.assertEqual([], phenos, 
                         "ideal is empty because model does not exist")        

    def test_ideal_match4(self):
        """Test for self-cross oomap."""
        
        ideal0 = self.oomap2.ideal_match_phenotypes(["hp0"])                              
        ideal1 = self.oomap2.ideal_match_phenotypes(["hp1"])
        ideal01 = self.oomap2.ideal_match_phenotypes(["hp0", "hp1"])
        
        self.assertEqual(["hp0"], ideal0, 
                         "ideal should return itself for hp0")        
        self.assertEqual(["hp1"], ideal1, 
                         "ideal should return itself for hp1")
        self.assertEqual(["hp0", "hp1"], sorted(ideal01), 
                         "ideal should return itself hp0, hp1")


class OrthologMappingTests(unittest.TestCase):
    """Test cases for class OntologyOntologyMap."""
                                    
    def test_mapEmpty(self):                             
        """Check for mapping with incomplete data."""
        mm = dss.OrthologMapping()
        mm.addPair("mouse", "a")
        expected = []
        output = mm.getPairs("mouse", "human")
        
        self.assertEqual(output, expected, "no orthology")        

    def test_one2one(self):
        """Check mapping for one-to-one orthology."""
        
        mm = dss.OrthologMapping()
        mm.addPair("mouse", "a")
        mm.addPair("human", "x")
        expected1 = [["a", "x"]]
        output1 = mm.getPairs("mouse", "human")
        expected2 = [["x", "a"]]
        output2 = mm.getPairs("human", "mouse")
        self.assertEqual(expected1, output1, "one-to-one orthologs, mouse-human")
        self.assertEqual(expected2, output2, "one-to-one orthologs, human-mouse")
        
        
    def test_many2many(self):
        
        mm= dss.OrthologMapping()
        mm.addPair("mouse", "a")
        mm.addPair("mouse", "b")
        mm.addPair("human", "x")
        mm.addPair("human", "y")
        expected = [["a", "x"], ["a", "y"], ["b", "x"], ["b", "y"]]
        output = sorted(mm.getPairs("mouse", "human"))
        self.assertEqual(expected, output, "many-to-many orthologs")
        
                
class GeneInfoTests(unittest.TestCase):
    """Test cases for class GeneInfo."""
    
    def test_official1(self):
        """Check tracking of official gene symbols."""
        
        gg = dss.GeneInfo("id0", "human")
        self.assertEqual(gg.official_symbol, None, 
                         "symbol should be empty at first")
        gg.addSymbol("NAME0")
        self.assertEqual(gg.official_symbol, "NAME0", 
                         "adding first symbol")
        self.assertEqual(gg.symbols, set(), 
                         "first symbol is official symbol")
        gg.addSymbol("NAME1")
        self.assertEqual(gg.official_symbol, "NAME1", 
                         "adding new symbol replaces previous official")
        self.assertEqual(gg.symbols, set(["NAME0"]), 
                         "the old official symbol should move to alt set")
        
        
    def test_official2(self):
        """Check adding of withdrawn symbols."""
        
        gg = dss.GeneInfo("id0", "human", "Official")
        gg.addSymbol(None) # should not do anything        
        gg.addSymbol("Withdrawn", True)
        self.assertEqual(gg.official_symbol, "Official", 
                         "adding withdrawn symbol does not affect official")
        self.assertEqual(gg.symbols, set(["Withdrawn"]), 
                         "withdrawn symbol should appear in .symbols")
        
        
        
class MapSetsTests(unittest.TestCase):
    """Test cases for class MapSets."""
    
    def setUp(self):
        """Create a MapSets with a few elements."""
        
        self.ms = dss.MapSets()
        self.ms.add("a", "x")
        self.ms.add("a", "y")
        self.ms.add("a", "x")
        self.ms.add("b", "y")
        
    
    def test_has(self):
        """Check creation of simple maps sets."""
                
        self.assertEqual(self.ms.has("a"), True, 
                         "instance has key that has been inserted")
        self.assertEqual(self.ms.has("x"), False, 
                         "instance does not have other keys ")
        self.assertEqual(self.ms.has("A"), False, 
                         "has is case sensitive")
    
        
    def test_get(self):
        """Check extraction of sets."""
            
        self.assertEqual(self.ms.get("a"), set(["x", "y"]), "a has two distinct elements")
        self.assertEqual(self.ms.get("c"), None, "return none when key is not in map")        
        
            
    def test_reverse(self):
        """Check reversal of sets."""
        
        output = self.ms.reverse()
            
        expected = dss.MapSets()
        expected.add("x", "a")
        expected.add("y", "a")
        expected.add("y", "b")
                
        self.assertEqual(expected.get("a"), output.get("a"), "a not in reverse")
        self.assertEqual(expected.get("x"), output.get("x"), "x should link to one element")
        self.assertEqual(expected.get("y"), output.get("y"), "y should link to two elements")
        
        
  
        