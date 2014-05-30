"""
methods for preconditioning the rate matrix to improve the performance 
of the linear algebra solvers
"""
from itertools import izip

import numpy as np
import scipy
import networkx as nx


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
        
    
    def matricize_eigenvectors(self, evecs):
        nodes = sorted(self.Ei.iterkeys())
        node2index = dict(((n, i) for i, n in enumerate(nodes)))
        eigenvectors = np.zeros([len(nodes), len(nodes)])
        eigenvectors[:,0] = 1. / np.sqrt(len(nodes))
        for k, evec in evecs.iteritems():
            for node, v in evec.iteritems():
                i = node2index[node]
                eigenvectors[i,k] = v
            eigenvectors[:,k] /= np.linalg.norm(eigenvectors[:,k])
        self.node_list = nodes
        self.node2index = node2index
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
        print "equilibrium occupation probability"
        print self.pi
        
    
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
        sinks = set()
        sinks.add(s1)
        barrier_function[s1] = 0.
        escape_function[s1] = 0.
        
        self.compute_barrier_function(mst, s1, barrier_function, escape_function)
            
        for k in xrange(1, len(self.Ei)):
            print ""
            # find the next sink
            s = max(escape_function.iterkeys(), key=lambda i:escape_function[i])
#            sold = slist[-1]
#            print slist
            print "current sink", s
            print "past sinks", sinks
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
            
            use_alternate_evec = False
            if use_alternate_evec:
                evecs[k] = evec_alternate
            
            if True:
                print "evecs[",k,"] = ", evecs[k]
                print "alternate      ", evec_alternate
            
            # recompute barrier_function and escape_function
            self.recompute_barrier_function(mst, s, barrier_function, escape_function)
            
            sinks.add(s)
        
        # process eigenvalues
        print delta
#        evals = np.cumsum(delta)
        self.eigenvalues = np.exp(-delta / self.T)
        self.eigenvalues[0] = 0.
        
        # process eigenvectors
        self.eigenvectors = self.matricize_eigenvectors(evecs)
        
        # make the equilibrium occupation probabilities
        self._make_equilibrium_occupation_probabilities()
        
        # test to see if the eigenvectors are orthonormal
        print "\ntesting to see if the eigenvectors are normalized"
        for k in xrange(len(self.Ei)):
            v = self.eigenvectors[:,k]
            val = np.dot(v, v*self.pi)
            self.eigenvectors[:,k] /= np.sqrt(val)
            print "k", k, ":", np.dot(v, v*self.pi)
        print "testing to see if the eigenvectors are orthogonal"
        for k1 in xrange(len(self.Ei)):
            for k2 in xrange(k1):
                v1 = self.eigenvectors[:,k1]
                v2 = self.eigenvectors[:,k2]
                print "k1 k2", k1, k2, ":", np.dot(v1, v2*self.pi)
            
        
        print "\neigenvalues", self.eigenvalues

        print "eigenvectors"
        print self.eigenvectors

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
            print v.shape, mat.shape, pi.shape
            mat += -evals[k] * np.outer(v, v * pi)
            mat_inv += 1./evals[k] * np.outer(v, v * pi)
        print "approximate K matrix"
        print mat
        print "approximate inv K", 
        mat_inv
        print "approximate K * Kinv"
        mat_inv.dot(mat)

        if True:       
            m, eval, evec = get_eigs(self.spect.Ei, self.spect.Eij, T=self.spect.T)
            print "exact K matrix"
            print m
            print "condition number of m", np.linalg.cond(m[:-1, :-1])

        Kcond = np.dot(mat_inv[:-1,:-1], m[:-1,:-1])
        print "condition number of conditiond matrix", np.linalg.cond(Kcond)
            
 
#
# only testing below here
#
from tests.test_preconditioning import make_random_energies_complete, get_eigs

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
        k = 0
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
        
#        print m


def test2():
    np.random.seed(1)
    Ei, Eij = make_random_energies_complete(8)
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
#    np.random.seed(1)
    Ei, Eij = make_random_energies_complete(3)
    T = .1
    precond = MSTPreconditioning(Ei, Eij, T=T)
    
    


if __name__ == "__main__":
#    test2()
    test_precond1()
