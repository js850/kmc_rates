import unittest
from test_graph_transformation import _three_state_rates, _MakeRandomGraph

from kmc_rates.rates_linalg import CommittorLinalg, MfptLinalgSparse, TwoStateRates

class TestLinalg3(unittest.TestCase):
    def setUp(self):
        self.rates = _three_state_rates()
        # all rates after graph renormalization should be 1.0
        self.final_rate = 1.0

    def _test_rate(self, i, j):
        reducer = TwoStateRates(self.rates, [i], [j])
        reducer.compute_rates()
        rAB = reducer.get_rate_AB() 
        self.assertAlmostEqual(rAB, self.final_rate, 7)
        
        reducer.compute_committors()
        rAB_ss = reducer.get_rate_AB_SS()
        print "kSS", rAB_ss
        self.assertAlmostEqual(rAB_ss, 1.5, 7)
        
    def test01(self):
        self._test_rate(0,1)

    def test12(self):
        self._test_rate(1,2)

    def test02(self):
        self._test_rate(0,2)

class TestLinalgRandom(unittest.TestCase):
    def do_check(self, A, B, nnodes=20, nedges=20):
        maker = _MakeRandomGraph(nnodes=20, nedges=20, node_set=A+B)
        rates = maker.make_rates()
        reducer = MfptLinalgSparse(rates, B)
        reducer.compute_mfpt()

    def test(self):
        A, B = [0], [1]
        self.do_check(A, B)
 
    def test_setA(self):
        A, B = [0, 1, 2], [3]
        self.do_check(A, B)
  
    def test_setAB(self):
        A, B = [0, 1, 2], [3, 4, 5, 6]
        self.do_check(A, B)

    def test_weakly_connected(self):
        rates = {}
        nnodes = 5
        for n in range(nnodes-1):
            rates[(n,n+1)] = 1.#np.random.rand()           
            rates[(n+1,n)] = 1.#np.random.rand()
        # add a disconnected edge
#         rates[(nnodes, nnodes+1)] = 1.
#         rates[(nnodes+1, nnodes)] = 1.
        
        for a in xrange(nnodes-1):
            calc = MfptLinalgSparse(rates, [a])
            times = calc.compute_mfpt()
            self.assertGreater(min(times.itervalues()), 0)


if __name__ == "__main__":
    unittest.main()

