'''
Tests for contents of pd2/tools.py
(functions with generic tools)
'''


from os.path import join, dirname
import unittest
from pd2 import tools as pd2tools
from pd2 import loaddiseases as pd2diseases


s2c = pd2diseases.omimTitle

class ToolsTests(unittest.TestCase):
    """Test for pd2 tools."""
        
    def setUp(self):
        """Setup for individual tests."""
        pass                                                
                                                                                        
    def test_captilize_simple(self): 
        input1 = "title in lowercase"
        expected1 = "Title In Lowercase"     
        input2 = "LONG TITLE IN UPPER"
        expected2 = "Long Title In Upper"                   
        self.assertEqual(expected1, s2c(input1), "lower to caps")
       
    def test_capitalize_hyphens(self):          
        input1 = "abc-xyz; GENE"
        expected1 = "Abc-Xyz"
        input2 = "ABC-XYZ; GENEabc"
        expected2 = "Abc-Xyz"
        self.assertEqual(expected1, s2c(input1), "lower to caps")
        self.assertEqual(expected2, s2c(input2), "upper to caps")


