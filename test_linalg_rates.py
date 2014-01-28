import unittest
import numpy as np

from rates_linalg import MfptLinalgSparse

class TestMfptLinalgSparse(unittest.TestCase):
    
    def test1(self):
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
            print np.linalg.inv(calc.matrix.todense())
            self.assertGreater(times.min(), 0)
            print times[0]


        
if __name__ == "__main__":
    unittest.main()