'''
Tests for contents of pd2/scoring.py
(Classes to produce scores)
'''


from os.path import join, dirname
import unittest
from pd2 import scoring
from pd2 import dss


class AvgMaxScoreTests(unittest.TestCase):
    """Test cases for class AvgMaxScore."""
        
    def setUp(self):
        """Setup for individual tests."""
        pass                                                
                                                                                        
    def test_avg(self): 
        ms = scoring.AvgMaxScore()
        ms.add(3)                       
        self.assertEqual(3, ms.avg_score, "one-time add (avg)")
        self.assertEqual(3, ms.max_score, "one-time add (max)")
        ms.add(5)                       
        self.assertEqual(4, ms.avg_score, "two adds (avg)")
        self.assertEqual(5, ms.max_score, "two adds (max)")

    def test_ignore(self):          
        ms = scoring.AvgMaxScore()
        ms.add(3)
        ms.add(0)        
        self.assertEqual(3, ms.avg_score, "adding zeros should not affect average")
    
    
    
class PhenodigmModelDiseaseTests(unittest.TestCase):
    """Test cases for class PhenodigmModelDisease."""
    
    def assertDelta(self, first, second, msg, delta=0.1):
        """Shortcut function for testing float similarity."""
        self.assertAlmostEqual(first, second, msg=msg, delta=delta)
    
    def setUp(self):
        """Setup common objects for individual tests."""
        
        ## set up a cross-ontology map        
        a1 = {"query": "hp0", "match": "mp1",
                "ic": 1, "simJ": 1, "lcs": "xyz"} ## score=1
        a2 = {"query": "hp0", "match": "mp2",
                "ic": 4, "simJ": 4, "lcs": "xyz"} ## score=4
        a3 = {"query": "hp1", "match": "mp1",
                "ic": 1, "simJ": 4, "lcs": "xyz"} ## score=2
        self.oomap = dss.OntologyOntologyMap()
        self.oomap.add(a1)
        self.oomap.add(a2)
        self.oomap.add(a3)
        
        # set up disease phenotypes
        self.phenos = dict()        
        self.phenos["Dw"] = ["hp0"]
        self.phenos["Dx"] = ["hp1"]
        self.phenos["Dy"] = ["hp0", "hp1"]
        self.phenos["D99"] = ["hp99"] ## a phenotype not in oomap        
        
        # set up a model phenotypes
        self.phenos["Mi"] = ["mp1"]
        self.phenos["Mj"] = ["mp2"]
        self.phenos["M99"] = ["mp99"] ## a phenotype not in oomap
        
        ## precompute ideal best matches
        self.ideal = dict()
        impfun = self.oomap.ideal_match_phenotypes
        for x in ["Mi", "Mj", "M99"]:
            self.ideal[x] = impfun(self.phenos[x])        
        
        ## set a scoring object
        self.pmd = scoring.PhenodigmScoring(self.oomap)
                

    def atest_score_nomatches(self):
        """Tests between non-matching diseases/models."""
        
        Mi = self.phenos["Mi"]
        M99 = self.phenos["M99"]        
        D99 = self.phenos["D99"]        
        Dx = self.phenos["Dx"]
        
        scores = self.pmd.score(Mi, D99)
        self.assertEqual(0, scores[0], "disease without valid hp; avg")
        self.assertEqual(0, scores[1], "disease without valid hp; max")
        self.assertEqual(0, scores[2], "disease without valid hp; mean")
        
        scores = self.pmd.score(M99, Dx)
        self.assertEqual(0, scores[0], "model with unmapped mp; avg")
        self.assertEqual(0, scores[1], "model with unmapped mp; max")
        self.assertEqual(0, scores[2], "model with unmapped mp; mean")
        
    
    def test_score_nonexistent(self):
        """Test with some keys that do not exist."""
        
        Mi = self.phenos["Mi"]
        Mi_ideal = self.ideal["Mi"]
        Dx = self.phenos["Dx"]
        empty = []
        
        scores = self.pmd.score(Mi, empty, Mi_ideal)
        self.assertEqual(0, scores[2], "non-existent disease")
        scores = self.pmd.score(empty, Dx, empty)
        self.assertEqual(0, scores[2], "non-existent model")
        scores = self.pmd.score(empty, empty, empty)
        self.assertEqual(0, scores[2], "non-existent model and disease")
                    
    def test_score1(self):
        """Test for 1-1 valid association."""
        
        Mj = self.phenos["Mj"]
        Dw = self.phenos["Dw"]
        Mj_ideal = self.ideal["Mj"]
        
        ## Mj -> mp2 -> hp0
        ## Dw -> hp0        
        scores = self.pmd.score(Mj, Dw, Mj_ideal)
        self.assertEqual(100, scores[2], "Mj and Dw are one-term easy matches")
        
    def test_score2(self):
        """Test for valid but non-matching associations."""
        
        Mj = self.phenos["Mj"]
        Mj_ideal = self.ideal["Mj"]
        Dx = self.phenos["Dx"]
        
        ## Mj -> mp2 -> hp0
        ## Dx -> hp1
        scores = self.pmd.score(Mj, Dx, Mj_ideal)
        self.assertEqual(0, scores[2], "Mj and Dx are non-matches")
                
    def test_score3(self):
        """Test for a simple model, complex disease."""
        
        Mj = self.phenos["Mj"]
        Mj_ideal = self.ideal["Mj"]
        Dy = self.phenos["Dy"]
        
        ## Mj -> mp2 -> hp0
        ## Dy -> hp0, hp1 -> partial match
        
        scores = self.pmd.score(Mj, Dy, Mj_ideal)        
        self.assertDelta(66.66, scores[0],                         
            msg="Mj and Dy are partial matches (avg)")
        self.assertDelta(100, scores[2],
            msg="Mj and Dy are partial matches (max)")
        
        
    def test_score4(self):
        """Test model with multiple inverse matches."""
        
        Mi = self.phenos["Mi"]
        Mi_ideal = self.ideal["Mi"]
        Dx = self.phenos["Dx"]
        Dw = self.phenos["Dw"]
        
        ## Mi -> mp1 -> hp1 is best, but also links to hp0
        ## Dx -> hp1 (perfect match with Mi)
        ## Dw -> hp0 (matches with Mi but not perfect)

        scores = self.pmd.score(Mi, Dx, Mi_ideal)   
        self.assertDelta(100, scores[2], msg="Match because disease has hp1")
        
        scores = self.pmd.score(Mi, Dw, Mi_ideal)        
        self.assertDelta(50, scores[2], 
                         msg="Lower score because disease has hp0 (max)")
        self.assertDelta(50, scores[0],  
                         msg="Lower score because because disease has hp0 (avg)")

    def test_common(self):
        """Test common phenotype results."""

        Mi = self.phenos["Mi"]
        Dx = self.phenos["Dx"]
        Dw = self.phenos["Dw"]                                                
        
        self.assertEqual(sorted(["hp1_mp1"]), self.pmd.common_phenos(Mi, Dx))
        self.assertEqual(sorted(["hp0_mp1"]), self.pmd.common_phenos(Mi, Dw))
        