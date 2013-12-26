import unittest

import example_3state
import example_random_graph
from test_graph_transformation import *
from test_kmc import *

class TestExamples(unittest.TestCase):
    def test1(self):
        example_3state.main()
    
    def test2(self):
        example_random_graph.main()
        
    def test3(self):
        example_random_graph.readme_example()

if __name__ == "__main__":
    unittest.main()