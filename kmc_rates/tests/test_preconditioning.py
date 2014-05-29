import unittest
import itertools

import numpy as np

from kmc_rates._preconditioning import MSTSpectralDecomposition

class TestSpectralDecomp3(unittest.TestCase):
    def setUp(self):
        Ei = dict()
        Ei[3] = 0.
        Ei[1] = 2.
        Ei[2] = 3.
        Eij = dict()
        Eij[(1,2)] = 4.
        Eij[(2,3)] = 5.
        self.Ei = Ei
        self.Eij = Eij
    
    def test1(self):
        T = 0.1
        spect = MSTSpectralDecomposition(self.Ei, self.Eij, T=T)
        eval = spect.eigenvalues
        
        self.assertEqual(len(eval), 3)
        self.assertAlmostEqual(eval[0], 0, 4)
        self.assertAlmostEqual(eval[1], 1.26633e-14, 4)
        self.assertAlmostEqual(eval[2], 1.67027379804344e-05, 4)
        
        from numpy import sqrt
        evecs = [[1./sqrt(3), 1./sqrt(2), 0],
                 [1./sqrt(3), 1./sqrt(2), 1],
                 [1./sqrt(3), 0./sqrt(2), 0]]
        evecs = np.array(evecs)
        
        for v1, v2 in itertools.izip(evecs.flatten(), spect.eigenvectors.flatten()):
            self.assertAlmostEqual(v1, v2, 3)


if __name__ == "__main__":
    unittest.main()