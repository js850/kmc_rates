import unittest
import itertools

import numpy as np

from kmc_rates._preconditioning import MSTSpectralDecomposition


def make_random_energies_complete(nnodes):
    Ei = {}
    Eij = {}
    for i in xrange(nnodes):
        Ei[i] = np.random.uniform(-1,1)
    for i in xrange(nnodes):
        for j in xrange(i):
            Eij[(j,i)] = max(Ei[i], Ei[j]) + np.random.uniform(.1, 1)
    return Ei, Eij 

def make_rate_matrix(Ei, Eij, T=.05):
    node_list = sorted(Ei.iterkeys())
#    node2i = dict([(node,i) for i, node in enumerate(node_list)])

    n = len(Ei)
    m = np.zeros([n,n])
    for i in xrange(n):
        ni = node_list[i]
        for j in xrange(n):
            if i == j: continue
            nj = node_list[j]
            try:
                Ets = Eij[(ni,nj)]
            except KeyError:
                try:
                    Ets = Eij[(nj,ni)]
                except KeyError:
                    # there is no edge i,j
                    continue
            m[i,j] = np.exp(-(Ets - Ei[ni])/T)
    for i in xrange(n):
        m[i,i] = - m[i,:].sum()
    return m

def get_eigs(Ei, Eij, T=0.05):
    from pele.utils.hessian import sort_eigs
    m = make_rate_matrix(Ei, Eij, T=T)
    lam, v = np.linalg.eig(m)
    lam, v = sort_eigs(lam, v, reverse=True)
    print "exact eigenvalues", -lam
    print "exact eigenvectors"
    print v
    return m, lam, v
#    print v


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

class TestSpectralDecompRandom(unittest.TestCase):
    def test(self):
        n = 20
        T = 0.02
        Ei, Eij = make_random_energies_complete(n)
        spect = MSTSpectralDecomposition(Ei, Eij, T=T)
        m, evals, evecs = get_eigs(Ei, Eij, T=T)
        for v1, v2 in itertools.izip(spect.eigenvalues, -evals):
            self.assertAlmostEqual(v1, v2, 2)

#         for v1, v2 in itertools.izip(spect.eigenvectors.flatten(), evecs.flatten()):
#             self.assertAlmostEqual(v1, v2, 1)


if __name__ == "__main__":
    unittest.main()