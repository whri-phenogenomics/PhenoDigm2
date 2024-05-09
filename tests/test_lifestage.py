'''
Tests for contents of pd2/scoring.py
(Classes to produce scores)
'''

import unittest
from os.path import join, dirname
from pd2.loadmodels import load_life_stages

life_stage_file = join("resources", "annotations", "impc_life_stages.csv")


class LifeStageTests(unittest.TestCase):
    """Test cases for class LifeStage."""                    
                                                                                                                                
    def test_embryo(self):
        """several life stage accession ids should map to embryo"""
        
        translator = load_life_stages(life_stage_file)
        self.assertEqual(translator["IMPCLS:0001"], "embryo") 
        self.assertEqual(translator["IMPCLS:0002"], "embryo")
        self.assertEqual(translator["IMPCLS:0003"], "embryo")
        self.assertEqual(translator["IMPCLS:0004"], "embryo")
        
    def test_late(self):
        """one accession id should map to a 'late' group name."""
        
        translator = load_life_stages(life_stage_file)
        self.assertEqual(translator["IMPCLS:0007"], "late")
        
