import unittest
import numpy as np

from rates_linalg import MfptLinalgSparse

class TestMfptLinalgSparse(unittest.TestCase):
    
    def test1(self):
        rates = {}
        nnodes = 50
        for n in range(nnodes-1):
            rates[(n,n+1)] = np.random.rand()           
            rates[(n+1,n)] = np.random.rand()
        
        for a in xrange(nnodes):
            calc = MfptLinalgSparse(rates, [a])
            times = calc.compute_mfpt()
            self.assertGreater(times.min(), 0)


        
if __name__ == "__main__":
    unittest.main()