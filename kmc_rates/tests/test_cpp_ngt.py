import unittest
import numpy as np
import networkx as nx

from kmc_rates.ngt import NGT
from test_graph_transformation import _MakeRandomGraph, _three_state_rates



class TestNgtCpp3(unittest.TestCase):
    def setUp(self):
        self.rates = _three_state_rates()
        # all rates after graph renormalization should be 1.0
        self.final_rate = 1.0
        self.final_rate_SS = 1.5

    def _test_rate(self, i, j):
        reducer = NGT(self.rates, [i], [j], debug=True)
        reducer.compute_rates()
        rAB = reducer.get_rate_AB()
        rBA = reducer.get_rate_BA()
        self.assertAlmostEqual(rAB, self.final_rate, 7)
        self.assertAlmostEqual(rBA, self.final_rate, 7)

        rAB_SS = reducer.get_rate_AB_SS()
        rBA_SS = reducer.get_rate_BA_SS()
        self.assertAlmostEqual(rAB_SS, self.final_rate_SS, 7)
        self.assertAlmostEqual(rAB_SS, self.final_rate_SS, 7)

    def test01(self):
        self._test_rate(0,1)

    def test12(self):
        self._test_rate(1,2)

    def test02(self):
        self._test_rate(0,2)

class TestNgtCppRandom(unittest.TestCase):
    def do_test(self, A, B, nnodes=20, nedges=20):
        maker = _MakeRandomGraph(nnodes=20, nedges=20, node_set=A+B)
        maker.run()
        reducer = NGT(maker.rates, A, B, debug=False)  
        reducer.compute_rates()
        rAB = reducer.get_rate_AB()
        rBA = reducer.get_rate_BA()
             
    def test(self):
        A, B = [0], [1]
        self.do_test(A, B)
 
    def test_setA(self):
        A, B = [0, 1, 2], [3]
        self.do_test(A, B)
  
    def test_setAB(self):
        A, B = [0, 1, 2], [3, 4, 5, 6]
        self.do_test(A, B)


if __name__ == "__main__":
    unittest.main()
