"""
methods for preconditioning the rate matrix to improve the performance 
of the linear algebra solvers
"""
from itertools import izip
import sys

import numpy as np
import scipy
import networkx as nx

def print_matrix(m, fmt="%9.3g"):
    print np.array_str(m, precision=3, max_line_width=300)
    return
    np.savetxt(sys.stdout, m, fmt=fmt)

class MSTSpectralDecomposition(object):
    """spectral decomposition of a graph using the minimum spanning tree in the limit of small T"""
    def __init__(self, Ei, Eij, T=.1):
        self.Ei = Ei
        self.Eij = Eij
        self.T = T
        
        self.sorted_indices = sorted(((E, i) for i, E in self.Ei.iteritems()))
        print self.sorted_indices
        self.sorted_indices = [i for E, i in self.sorted_indices]
        
        self.run()

    def compute_mst(self):
        self.graph = nx.Graph()
        for i, E in self.Ei.iteritems():
            self.graph.add_node(i, E=E)
        for (i, j), E in self.Eij.iteritems():
            self.graph.add_edge(i, j, E=E)
        self.mst = nx.minimum_spanning_tree(self.graph, weight='E')
        return self.mst
    
    def get_edge_energy(self, u, v):
        try:
            return self.Eij[(u,v)]
        except KeyError:
            return self.Eij[(v,u)]
    
    def compute_barrier_function(self, mst, s, barrier_function, escape_function):
        barrier_function[s] = 0.
        escape_function[s] = 0. 
        for parent, child in nx.bfs_edges(mst, s):
            E = self.get_edge_energy(parent, child)
            if parent == s or E > barrier_function[parent]:
                barrier_function[child] = E
            else:
                barrier_function[child] = barrier_function[parent]
            escape_function[child] = barrier_function[child] - self.Ei[child]
    
    def recompute_barrier_function(self, mst, s, barrier_function, escape_function):
        unew = dict()
        vnew = dict()
        self.compute_barrier_function(mst, s, unew, vnew)
        for i, unewi in unew.iteritems():
            if unewi < barrier_function[i]:
                barrier_function[i] = unewi
            escape_function[i] = min(escape_function[i], barrier_function[i] - self.Ei[i])
        barrier_function[s] = 0.
        escape_function[s] = 0.
            
    def get_connected_sink(self, mst, s1, sinks):
        s2 = None
        for parent, child in nx.bfs_edges(mst, s1):
            if child in sinks:
                s2 = child
                break
        if s2 is None:
            print nx.connected_components(mst)
        assert s2 is not None
        return s2
    
    def find_cutting_edge_old(self, mst, s1, barrier_function):
        u1 = barrier_function[s1]
        i = j = None
        for parent, child in nx.bfs_edges(mst, s1):
            print "u[parent], u[child]", s1, ":", parent, child, barrier_function[parent], barrier_function[child], u1
            if barrier_function[child] < u1:
                assert barrier_function[parent] == u1
                i, j = child, parent
                break
        if i is None:
            raise Exception("error in finding cutting edge")
        E = self.get_edge_energy(i, j)
        print "cutting edge", i, j, E
        return E, i, j
    
    def find_cutting_edge(self, mst, s1, barrier_function, sinks):
        s2 = self.get_connected_sink(mst, s1, sinks)
        path = nx.shortest_path(mst, s1, s2)
        print "path", path
        E, i, j = max( [(self.get_edge_energy(i,j), i, j) for i, j in izip(path, path[1:])] )
        assert barrier_function[i] == barrier_function[s1]
#         assert barrier_function[]
        return E, i, j
        
    
    def compute_kth_eigenvector_components(self, mst, s):
        """compute the non zero components of the kth eigenvector"""
        evec = dict()
        for i in nx.node_connected_component(mst, s):
            evec[i] = 1.
        return evec

    def compute_kth_eigenvector_components_committor(self, mst, s, sold):
        from kmc_rates.rates_linalg import CommittorLinalg
        rates = dict()
        for i, j in nx.bfs_edges(mst, s):
            Ets = self.get_edge_energy(i, j)
            rates[(i,j)] = np.exp(-(Ets - self.Ei[i]) / self.T)
            rates[(j,i)] = np.exp(-(Ets - self.Ei[j]) / self.T)
        com = CommittorLinalg(rates, [sold], [s])
        com.compute_committors()
        return com.committor_dict
    
    def improve_eigenvector(self, evec, eval):
        # this basically does one step of conjugate gradient minimization
        from tests.test_preconditioning import make_rate_matrix
        m = make_rate_matrix(self.Ei, self.Eij, T=self.T)
        # symmetrise m
        m = m * self.pi[:,np.newaxis]
        r = eval * evec - np.dot(m, evec)
        p = r
        alpha = np.dot(r,r) / np.dot(p, np.dot(m, p))
        evec += alpha * p
        if True:
            rnew = eval * evec - np.dot(m, evec)
#             assert np.linalg.norm(rnew) <= np.linalg.norm(r), "%g %g" % (np.linalg.norm(rnew), np.linalg.norm(r))
    
    def normalize_eigenvectors(self, eigenvectors):
        # normalize them
        for k in xrange(len(self.Ei)):
            v = eigenvectors[:,k]
            eigenvectors[:,k] /= np.sqrt(np.dot(v, v * self.pi))
    
    def matricize_eigenvectors(self, evecs, slist):
        nodes = sorted(self.Ei.iterkeys())
        node2index = dict(((n, i) for i, n in enumerate(nodes)))
        self.node_list = nodes
        self.node2index = node2index
        self._make_equilibrium_occupation_probabilities()
        
        eigenvectors = np.zeros([len(nodes), len(nodes)])
        eigenvectors[:,0] = 1.
        for k, evec in evecs.iteritems():
#             if True:
#                 # improve accuracy of the eigenvectors. (I'm not sure this is actually correct)
#                 assert min(evec.itervalues()) > 0
#                 pnorm = self.pi_dict[slist[k-1]]
#                 small = sum((self.pi_dict[node] for node in evec.iterkeys()))
#                 small = np.sqrt(small / pnorm)
#                 big = -1./small
#                 eigenvectors[node2index[slist[k-1]],k] = small
            for node, v in evec.iteritems():
                i = node2index[node]
#                 eigenvectors[i,k] = v * big
                eigenvectors[i,k] = v
#             self.improve_eigenvector(eigenvectors[:,k], k)
        
        print "eigenvectors before any fiddling"
        print_matrix(eigenvectors)
        self.normalize_eigenvectors(eigenvectors)

        return eigenvectors
    
    def _print_component(self, mst, s):
        for parent, child in nx.bfs_edges(mst, s):
            last = child
        # now do depth first search from child
        print last,
        for parent, child in nx.dfs_edges(mst, child):
            print "-", child,
        print ""
    
    def _make_equilibrium_occupation_probabilities(self):
        Ered = np.array([self.Ei[node] / self.T for node in self.node_list])
        Ered -= Ered.max()
        pi = np.exp(-Ered)
        pi /= pi.sum()
        self.pi = pi
        self.pi_dict = dict((self.node_list[i], p) for i, p in enumerate(pi))
        print "equilibrium occupation probability"
        print_matrix(self.pi)
    
    def _test_orthogonality(self): 
        # test to see if the eigenvectors are orthonormal
        print "\ntesting to see if the eigenvectors are normalized"
        for k in xrange(len(self.Ei)):
            v = self.eigenvectors[:,k]
            r = np.sum(v * v * self.pi)
            if abs(r-1) > .1:
                print "k", k, ":", r
        print "testing to see if the eigenvectors are orthogonal"
        for k1 in xrange(len(self.Ei)):
            for k2 in xrange(k1):
                v1 = self.eigenvectors[:,k1]
                v2 = self.eigenvectors[:,k2]
                r = np.sum(v1 * v2 * self.pi)
                if np.abs(r) > .01:
                    print "k1 k2", k1, k2, ":", r

        from tests.test_preconditioning import make_rate_matrix
        m = make_rate_matrix(self.Ei, self.Eij, T=self.T)
        print "testing to see if the the eigenvectors satisfy the eigenvalue equation"
        for k in xrange(len(self.eigenvalues)):
            v = self.eigenvectors[:,k]
            l = np.dot(v*self.pi, np.dot(m, v))
            print "k", k, ":", l, "=?=", self.eigenvalues[k], "normalized |diff|", np.abs((l-self.eigenvalues[k])/l)

    
    def run(self):
        mst = self.compute_mst()
        print self.mst.number_of_edges()
        
        # the barrier function gives the minimax barrier separating
        # a node from a sink node
        barrier_function = dict()
        # the escape function is usually u[i] - E[i]
        escape_function = dict() # escape function

        delta = np.zeros(len(self.Ei)) # used to compute eigenvalues
        evecs = dict() # eigenvectors
        
        k = 0
        s1 = self.sorted_indices[0]
        sinks = set([s1])
        slist = [s1]
        barrier_function[s1] = 0.
        escape_function[s1] = 0.
        
        self.compute_barrier_function(mst, s1, barrier_function, escape_function)
            
        for k in xrange(1, len(self.Ei)):
            print "\niteration", k
            # find the next sink
            s = max(escape_function.iterkeys(), key=lambda i:escape_function[i])
#            sold = slist[-1]
#            print slist
            print "current sink", s
            print "past sinks", slist
            print "total # edges", mst.number_of_edges()
            print "size of component connected to current sink", len(nx.node_connected_component(mst, s))
            self._print_component(mst, s)
            
            
            if True:
                s2 = self.get_connected_sink(mst, s, sinks)
            
            for i, ui in barrier_function.iteritems(): 
                print "  barrier_function[",i,"] = ", ui
            
            # find the new cutting edge
            Epq, p, q = self.find_cutting_edge(mst, s, barrier_function, sinks)
            print "p q Epq", p, q, Epq
            print "s Es   ", s, self.Ei[s], "sold Esold", s2, self.Ei[s2] 
            print "eval", "exp(-(", Epq, "-", self.Ei[s], ")/", self.T, ") = ", 
            print np.exp(-(Epq - self.Ei[s]) / self.T)
            
            # save the eigenvalue deltas
            delta[k] = Epq - self.Ei[s]
            
            # compute the eigenvector components in the more precise way
            evec_alternate = self.compute_kth_eigenvector_components_committor(mst, s, s2)
            
            # remove the cutting edge
            print "removing edge (",p,",",q,")"
            mst.remove_edge(p, q)
            
            # compute the eigenvector components in the simplified way
            evecs[k] = self.compute_kth_eigenvector_components(mst, s)
            
            if True:
                print "evec simple", evecs[k]
                print "  alternate", evec_alternate

            use_alternate_evec = False
            if use_alternate_evec:
                evecs[k] = evec_alternate
            
            
            # recompute barrier_function and escape_function
            self.recompute_barrier_function(mst, s, barrier_function, escape_function)
            
            sinks.add(s)
            slist.append(s)

        print ""
                
        # process eigenvalues
#        evals = np.cumsum(delta)
        self.eigenvalues = -np.exp(-delta / self.T)
        self.eigenvalues[0] = 0.
        
        # process eigenvectors
        self.eigenvectors = self.matricize_eigenvectors(evecs, slist)

        self._test_orthogonality()

        if True:  
            print "\ntrying to improve eigenvectors by CG minimization"      
            for k in xrange(len(self.eigenvalues)):
                self.improve_eigenvector(self.eigenvectors[:,k], self.eigenvalues[k])
            self.normalize_eigenvectors(self.eigenvectors)
            self._test_orthogonality()
        
    
        
        print "\neigenvalues"
        print_matrix(self.eigenvalues)

        print "eigenvectors (normalized to norm==1)"
        v = self.eigenvectors / (np.sqrt(np.sum(self.eigenvectors**2, axis=0)))[np.newaxis,:]
        print_matrix(v)
        print "eigenvectors (normalized to sum(v*v*pi) = 1)"
        print_matrix(self.eigenvectors)

class MSTPreconditioning(object):
    def __init__(self, Ei, Eij, T=0.1):
        self.spect = MSTSpectralDecomposition(Ei, Eij, T=T)
        self.run()
    
    def run(self):
        n = len(self.spect.Ei)
        evals = self.spect.eigenvalues
        evecs = self.spect.eigenvectors
        mat = np.zeros([n, n])
        mat_inv = np.zeros([n, n])
        pi = self.spect.pi
        
        for k in xrange(1,n):
            v = evecs[:,k]
            mat += evals[k] * np.outer(v, v * pi)
            mat_inv -= 1./evals[k] * np.outer(v, v * pi)
        print "approximate K matrix"
        print_matrix(mat)
        print "approximate inv K"
        print_matrix(mat_inv)

        if True:       
            m, eval, evec = get_eigs(self.spect.Ei, self.spect.Eij, T=self.spect.T)
            print "exact K matrix"
            print_matrix(m)
            print "condition number of m", np.linalg.cond(m[:-1, :-1])
        


        Kcond = np.dot(mat_inv, m)
        print "conditioned matrix"
        print_matrix(Kcond)
        print "condition number of conditioned matrix", np.linalg.cond(Kcond[:-1,:-1])
            
 
#
# only testing below here
#

#def test7():
#    E = dict()
#    E[1] = -13.
#    E[2] = -12.
#    E[3] = -10.5
#    E[5] = -9.2
#    E[6] = -9.
#    E[7] = -8.5
#    E[(2,3)] = -6
#    E[4] = -4.3
#    E[(4,3)] = -4.2
#    E[(5,6)] = -4.11 
#    E[(4,5)] = -4. 
#    E[(1,2)] = 
#    E[(1,4)] = 
#    E[(1,6)] = 

def test1():
    from numpy import exp
    E3 = 0.;
    E1 = 2.;
    E2 = 3.;
    E12 = 4.1;
    E23 = 5.2;
    
    Ei = dict()
    Ei[1] = E1
    Ei[2] = E2
    Ei[3] = E3
    Eij = dict()
    Eij[(1,2)] = E12
    Eij[(2,3)] = E23
    
    T = 0.1
    spect = MSTSpectralDecomposition(Ei, Eij, T=T)
    m, eval, evec = get_eigs(Ei, Eij, T=T)
    
    if True:
        print "\ntesting eigenvectors"
        k = 2
        v = spect.eigenvectors[:,k].copy()
        lam = -spect.eigenvalues[k]
        v /= np.linalg.norm(v)
        print lam, v
        print lam * v
        print m.dot(v)
        
        print "\ntesting exact eigenvectors"
        v = evec[:,k].copy()
        lam = eval[k]
        v /= np.linalg.norm(v)
        print lam, v
        print lam * v
        print m.dot(v)
        
#        print m


def test2():
    np.random.seed(1)
    Ei, Eij = make_random_energies_complete(4)
    print Ei.items()
    print Eij.items()
    T = .02
    spect = MSTSpectralDecomposition(Ei, Eij, T=T)
    
    m, eval, evec = get_eigs(Ei, Eij, T=T)

    if True:
        print "\ntesting eigenvectors"
        k = 1
        v = spect.eigenvectors[:,k].copy()
        lam = spect.eigenvalues[k]
        v /= np.linalg.norm(v)
        print lam, v
        print lam * v
        print m.dot(v)
        
        print "\ntesting exact eigenvectors"
        v = evec[:,k].copy()
        lam = eval[k]
        v /= np.linalg.norm(v)
        print lam, v
        print lam * v
        print m.dot(v)

def cond(m):
    u, s, v = np.linalg.svd(m)
    print "singular values"
    print s
    print s**2
    
def test_precond1():
    np.random.seed(1)
    Ei = dict()
    Ei[3] = 0.
    Ei[1] = 2.
    Ei[2] = 3.
    Eij = dict()
    Eij[(1,2)] = 4.1
    Eij[(2,3)] = 5.2

    T = .1
    precond = MSTPreconditioning(Ei, Eij, T=T)
    
def test_precond2():
    np.random.seed(3)
    Ei, Eij = make_random_energies_complete(8)
    T = .05
    precond = MSTPreconditioning(Ei, Eij, T=T)
    
    


if __name__ == "__main__":
    from tests.test_preconditioning import make_random_energies_complete, get_eigs

#     test2()
    test_precond1()
