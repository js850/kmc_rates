"""
methods for preconditioning the rate matrix to improve the performance 
of the linear algebra solvers
"""
from itertools import izip

import numpy as np
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
            if E > barrier_function[parent]:
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
            
    def get_connected_sink(self, mst, s1, barrier_function):
        s2 = None
        for parent, child in nx.bfs_edges(mst, s1):
            if barrier_function[child] == 0.:
                s2 = child
                break
        if s2 is None:
            print nx.connected_components(mst)
        assert s2 is not None
        return s2
    
    def find_cutting_edge(self, mst, s1, barrier_function):
        u1 = barrier_function[s1]
        for parent, child in nx.bfs_edges(mst, s1):
            print "u[parent], u[child]", s1, ":", parent, child, barrier_function[parent], barrier_function[child], u1
            if barrier_function[child] < u1:
                assert barrier_function[parent] == u1
                i, j = child, parent
                break
        E = self.get_edge_energy(i, j)
#        path = nx.shortest_path(mst, s1, s2)
#        print "path", path
#        E, i, j = max( [(self.get_edge_energy(i,j), i, j) for i, j in izip(path, path[1:])] )
        print "cutting edge", i, j, E
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
        return eigenvectors
    
    def _print_component(self, mst, s):
        for parent, child in nx.bfs_edges(mst, s):
            last = child
        # now do depth first search from child
        print last,
        for parent, child in nx.dfs_edges(mst, child):
            print "-", child,
        print ""
    
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
        for i, ui in barrier_function.iteritems():
            print "barrier_function escape_function", i, barrier_function, escape_function[i]
            
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
                s2 = self.get_connected_sink(mst, s, barrier_function)
            
            for i, ui in barrier_function.iteritems(): print "  barrier_function[",i,"] = ", ui
            
            # find the new cutting edge
            Epq, p, q = self.find_cutting_edge(mst, s, barrier_function)
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
        print "eigenvalues", self.eigenvalues
        
        # process eigenvectors
        self.eigenvectors = self.matricize_eigenvectors(evecs)
        print "eigenvectors"
        print self.eigenvectors
        
        

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
    T = 1.
    beta = 1./T
    E3 = 0.;
    E1 = 2.;
    E2 = 3.;
    E12 = 4.1;
    E23 = 5.2;
    k = dict()
    k[(1,2)] = exp(-beta*(E12-E1));
    k[(1,3)] = 0;
    k[(1,1)] = -k[(1,3)] - k[(1,2)];
    k[(2,1)] = exp(-beta*(E12-E2));
    k[(2,3)] = exp(-beta*(E23-E2));
    k[(2,2)] = -k[(2,1)] - k[(2,3)];
    k[(3,1)] = 0;
    k[(3,2)] = exp(-beta*(E23-E3));
    k[(3,3)] = -k[(3,1)] -k[(3,2)];
    
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

def test2():
    np.random.seed(0)
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
    

if __name__ == "__main__":
    test2()
