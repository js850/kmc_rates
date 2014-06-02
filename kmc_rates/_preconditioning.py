"""
methods for preconditioning the rate matrix to improve the performance 
of the linear algebra solvers
"""
from itertools import izip
import sys

import numpy as np
import scipy
import scipy.misc
import networkx as nx
import pylab as plt

def print_matrix(m, fmt="%9.3g"):
    print np.array_str(m, precision=3, max_line_width=300)
    return
    np.savetxt(sys.stdout, m, fmt=fmt)

class ViewMSTSpectralDecomp(object):
    def __init__(self, spect):
        self.spect = spect
        self.run()

    def draw_nodes(self, nodes, color='k'):
        xypos = [self.pos[n] for n in nodes]
        x = [xy[0] for xy in xypos]
        y = [xy[1] for xy in xypos]
        plt.scatter(x, y, c=color, linewidths=0, s=60)

    def draw_edges(self, edges, color='k', lw=1.):
        for u, v in edges:
            xy = np.array([self.pos[u], self.pos[v]])
            plt.plot(xy[:,0], xy[:,1], color, lw=lw)
    
    def draw(self, k):
        plt.clf()
        
        m=k+1
        self.draw_nodes(self.spect.slist[:m], color='b')
        self.draw_nodes([self.spect.slist[m]], color='r')
        if m < len(self.spect.slist):
            self.draw_nodes(self.spect.slist[m+1:])
        
        # draw edges
        self.draw_edges(self.spect.cut_edges[:k], lw=.3, color='k--')
        self.draw_edges([self.spect.cut_edges[k]], color='r--')
        if k < len(self.spect.slist):
            self.draw_edges(self.spect.cut_edges[(k+1):])

        # draw edges lightly that have been removed
        plt.show()
      
    def run(self):
        self.mst = self.spect._make_minimum_spanning_tree()
        self.pos = nx.spring_layout(self.mst)
        self.removed_edges = []
        
        for k in xrange(len(self.spect.Ei)-1):
            self.draw(k)
            
            self.mst.remove_edge(*self.spect.cut_edges[k])
            self.removed_edges = self.spect.cut_edges[:(k+1)]
        

class MSTSpectralDecompositionOld(object):
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
        mst = nx.minimum_spanning_tree(self.graph, weight='E')
        return mst
    
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
        
    
    def compute_kth_eigenvector_components(self, mst, s, p, q):
        """compute the non zero components of the kth eigenvector"""
        
        S1 = dict(((i, 1.) for i in nx.node_connected_component(mst, p)))
        S2 = dict(((i, 1.) for i in nx.node_connected_component(mst, q)))
        if s in S2:
            S1, S2 = S2, S1
        assert s in S1
        assert s not in S2
        return S1, S2

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
#         m = m / self.pi[np.newaxis,:]
        m = m * self.pi[:,np.newaxis]
#         u1, u2 = self.node_list[1], self.node_list[2]
#         print "symmetrized?", m[1,2], m[2,1], np.exp(-(self.get_edge_energy(u1, u2) - self.Ei[u2])/self.T) / self.pi[1]
#         print "symmetrized?", m[1,2], m[2,1], np.exp(-(self.get_edge_energy(u1, u2) - self.Ei[u2])/self.T) * self.pi[2]
        r = eval * evec - np.dot(m, evec * self.pi)
        p = r
        alpha = np.dot(r,r) / np.dot(p, np.dot(m, p))
        evec += alpha * p
        print "alpha", alpha, "r*r", r.dot(r), "p*m*p", np.dot(p, np.dot(m, p))
        if True:
            rnew = eval * evec - np.dot(m, evec * self.pi)
            assert np.linalg.norm(rnew) <= np.linalg.norm(r), "%g %g" % (np.linalg.norm(rnew), np.linalg.norm(r))

    def improve_eigenvectors_by_orthogonalization(self):
        # graeme schmidt orthogonalization
        from tests.test_preconditioning import make_rate_matrix
        for k in xrange(len(self.Ei)):
            v = self.eigenvectors[:,k]
            for j in xrange(k):
                u = self.eigenvectors[:,j]
                proj = np.sum(u * v * self.pi) / np.sum(u * u * self.pi)
                v -= u * proj

    def normalize_eigenvectors(self, eigenvectors):
        # normalize them
        for k in xrange(len(self.Ei)):
            v = eigenvectors[:,k]
            eigenvectors[:,k] /= np.sqrt(np.dot(v, v * self.pi))
    
    def _make_node_to_index(self):
        if hasattr(self, "node2index"):
            return
        nodes = sorted(self.Ei.iterkeys())
        node2index = dict(((n, i) for i, n in enumerate(nodes)))
        self.node_list = nodes
        self.node2index = node2index
        self._make_equilibrium_occupation_probabilities()
        
    
    def matricize_eigenvectors(self, evecs, slist):
        self._make_node_to_index()
        
        eigenvectors = np.zeros([len(self.node_list), len(self.node_list)])
        eigenvectors[:,0] = 1.
        for k, evec in evecs.iteritems():
            for node, v in evec.iteritems():
                i = self.node2index[node]
                eigenvectors[i,k] = v
        
        print "eigenvectors before any fiddling"
        print_matrix(eigenvectors)
        self.normalize_eigenvectors(eigenvectors)

        return eigenvectors
    
    def matricize_eigenvectors_new(self, S1list, S2list):
        self._make_node_to_index()

        eigenvectors = np.zeros([len(self.node_list), len(self.node_list)])
        eigenvectors[:,0] = 1.
        
        
        for k in xrange(1, len(self.Ei)):
            S1 = set(S1list[k].iterkeys())
            S2 = set(S2list[k].iterkeys())
            
            isum_pi1 = 1./sum((self.pi_dict[node] for node in S1))
            isum_pi2 = 1./sum((self.pi_dict[node] for node in S2))
            C1 =  isum_pi1 / np.sqrt(isum_pi1 + isum_pi2)
            C2 = -isum_pi2 / np.sqrt(isum_pi1 + isum_pi2)
            for node in S1:
                i = self.node2index[node]
                eigenvectors[i,k] = C1
            for node in S2:
                i = self.node2index[node]
                eigenvectors[i,k] = C2
        
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
        self.mst = self.compute_mst()
        mst = self.mst
        print self.mst.number_of_edges()
        
        # the barrier function gives the minimax barrier separating
        # a node from a sink node
        barrier_function = dict()
        # the escape function is usually u[i] - E[i]
        escape_function = dict() # escape function

        delta = np.zeros(len(self.Ei)) # used to compute eigenvalues
        evecs = dict() # eigenvectors
        evecs2 = dict() # eigenvectors
        
        k = 0
        s1 = self.sorted_indices[0]
        sinks = set([s1])
        slist = [s1]
        self.cut_edges = []
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
            evecs[k], evecs2[k] = self.compute_kth_eigenvector_components(mst, s, p, q)
            
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
            self.cut_edges.append((p,q))

        self.slist = slist
        print ""
                
        # process eigenvalues
#        evals = np.cumsum(delta)
        self.eigenvalues = -np.exp(-delta / self.T)
        self.eigenvalues[0] = 0.
        
        # process eigenvectors
#        self.eigenvectors = self.matricize_eigenvectors(evecs, slist)
        self.eigenvectors = self.matricize_eigenvectors_new(evecs, evecs2)

        self._test_orthogonality()

        if False:  
            print "\ntrying to improve eigenvectors by CG minimization"      
            for k in xrange(len(self.eigenvalues)):
                self.improve_eigenvector(self.eigenvectors[:,k], self.eigenvalues[k])
            self.normalize_eigenvectors(self.eigenvectors)
            self._test_orthogonality()
        
        if False:  
            print "\ntrying to improve eigenvectors by graeme schmidt orthogonalization"      
            self.improve_eigenvectors_by_orthogonalization()
            self.normalize_eigenvectors(self.eigenvectors)
            self._test_orthogonality()
        
    
        
        print "\neigenvalues"
        print_matrix(self.eigenvalues)

        print "eigenvectors (normalized to norm==1)"
        v = self.eigenvectors / (np.sqrt(np.sum(self.eigenvectors**2, axis=0)))[np.newaxis,:]
        print_matrix(v)
        print "eigenvectors (normalized to sum(v*v*pi) = 1)"
        print_matrix(self.eigenvectors)

class MSTSpectralDecomposition(object):
    """spectral decomposition of a graph using the minimum spanning tree in the limit of small T"""
    def __init__(self, Ei, Eij, T=.1, verbose=True):
        self.Ei = Ei
        self.Eij = Eij
        self.T = T
        self.verbose = verbose
        self.nnodes = len(self.Ei)
        self._make_node_to_index()
        self._make_equilibrium_occupation_probabilities()
        self.run()
 
    def _make_node_to_index(self):
        if hasattr(self, "node2index"):
            return
        nodes = sorted(self.Ei.iterkeys())
        node2index = dict(((n, i) for i, n in enumerate(nodes)))
        self.node_list = nodes
        self.node2index = node2index

    def _make_equilibrium_occupation_probabilities(self):
        log_pi = np.array([-self.Ei[node] / self.T for node in self.node_list])
        # normalize log_pi
        log_sum_pi = scipy.misc.logsumexp(log_pi)
        log_pi -= log_sum_pi
        self.pi = np.exp(log_pi)
        self.log_pi = log_pi
        
        self.log_pi_dict = dict((self.node_list[i], lpi) for i, lpi in enumerate(self.log_pi))
        
        print "equilibrium occupation probability"
        print_matrix(self.pi)

    def _get_edge_energy(self, u, v):
        try:
            return self.Eij[(u,v)]
        except KeyError:
            return self.Eij[(v,u)]

    def _make_minimum_spanning_tree(self):
        self.graph = nx.Graph()
        for i, E in self.Ei.iteritems():
            self.graph.add_node(i, E=E)
        for (i, j), E in self.Eij.iteritems():
            self.graph.add_edge(i, j, E=E)
        mst = nx.minimum_spanning_tree(self.graph, weight='E')
        return mst

    def _compute_barrier_function(self, mst, s, barrier_function, escape_function):
        """compute the barrier function and escape function
        
        compute the maximum barrier separating each node from the sink node s
        """
        barrier_function[s] = 0.
        escape_function[s] = 0. 
        for parent, child in nx.bfs_edges(mst, s):
            Ets = self._get_edge_energy(parent, child)
            if parent == s or Ets > barrier_function[parent]:
                barrier_function[child] = Ets
            else:
                barrier_function[child] = barrier_function[parent]
            escape_function[child] = barrier_function[child] - self.Ei[child]

    def _recompute_barrier_function(self, mst, s, barrier_function, escape_function):
        """recompute the barrier and escape function"""
        bar_new = dict()
        esc_new = dict()
        # compute the new barrier function for all nodes connected to s
        self._compute_barrier_function(mst, s, bar_new, esc_new)
        # if the new barrier function is less than the old then replace the old
        for i, bi in bar_new.iteritems():
            if bi < barrier_function[i]:
                barrier_function[i] = bi
            escape_function[i] = min(escape_function[i], barrier_function[i] - self.Ei[i])
        barrier_function[s] = 0.
        escape_function[s] = 0.

    def _get_connected_sink(self, mst, s, sinks):
        s2 = None
        for parent, child in nx.bfs_edges(mst, s):
            if child in sinks:
                s2 = child
                break
        assert s2 is not None
        return s2


    def _find_cutting_edge(self, mst, s1, barrier_function, sinks):
        """Find the cutting edge
        
        This is the edge with the maximum barrier in the path connecting s to 
        another sink
        """
        s2 = self._get_connected_sink(mst, s1, sinks)
        path = nx.shortest_path(mst, s1, s2)
        if self.verbose:
            print "path", path
        Ets, i, j = max(((self._get_edge_energy(i,j), i, j) for i, j in izip(path, path[1:])))
        assert barrier_function[i] == barrier_function[s1]
#         assert barrier_function[]
        return Ets, i, j

    def _make_kth_eigenvector(self, mst, k, p, q, s):
        """make the current eigenvector
        
        (p, q) is the cutting edge, and is already removed.  Assume that
        p is the node that is still connected to the new sink s.
        The eigenvector v[i] will be C1 for all nodes i connected to p, and
        C2 for all nodes i connected to q.  It will be zero elsewhere.
        C1 and C2 are determined by normalization and orthogonalization. In
        particular, we force orthogonalization with respect to the 0th 
        eigenvector, which is just np.ones(nnodes) 
        """
        S1 = nx.node_connected_component(mst, p)
        S2 = nx.node_connected_component(mst, q)
        if s in S2:
            S1, S2 = S2, S1
        assert s in S1 and s not in S2
        
        # compute the inverse sum of pi for all nodes in S1 and S2
        isum_pi1 = - scipy.misc.logsumexp([self.log_pi_dict[node] for node in S1])
        isum_pi2 = - scipy.misc.logsumexp([self.log_pi_dict[node] for node in S2])
        norm = scipy.misc.logsumexp([isum_pi1, isum_pi2])
        logC1 = isum_pi1 - norm / 2.
        logC2 = isum_pi2 - norm / 2.
        is1 = np.array([self.node2index[node] for node in S1])
        is2 = np.array([self.node2index[node] for node in S2])
        
        self.eigenvectors[:,k] = 0.
        self.eigenvectors[is1, k] =  np.exp(logC1)
        self.eigenvectors[is2, k] = -np.exp(logC2)
        
        self.log_eigenvectors[:,k] = 1.
        self.log_eigenvectors[is1, k] = logC1
        self.log_eigenvectors[is2, k] = logC2
        
        self.eigenvector_indicator[:,k] = 0
        self.eigenvector_indicator[is1,k] =  1
        self.eigenvector_indicator[is2,k] = -1

    def _test_eigenvectors(self): 
        # test to see if the eigenvectors are orthonormal
        print "\ntesting to see if the eigenvectors are normalized sum(v*v*pi) == 1"
        for k in xrange(len(self.Ei)):
            v = self.eigenvectors[:,k]
            r = np.sum(v * v * self.pi)
            if abs(r-1) > .1:
                print "k", k, ":", r
        print "testing to see if the eigenvectors are orthogonal sum(u*v*pi) == 0"
        for k1 in xrange(len(self.Ei)):
            for k2 in xrange(k1):
                v1 = self.eigenvectors[:,k1]
                v2 = self.eigenvectors[:,k2]
                r = np.sum(v1 * v2 * self.pi)
                if np.abs(r) > .01:
                    print "k1 k2", k1, k2, ":", r

        from tests.test_preconditioning import make_rate_matrix
        m = make_rate_matrix(self.Ei, self.Eij, T=self.T)
        print "testing to see if the the eigenvectors satisfy the eigenvalue equation dot(v*pi, dot(m, v)) = lambda"
        for k in xrange(len(self.eigenvalues)):
            v = self.eigenvectors[:,k]
            l = np.dot(v*self.pi, np.dot(m, v))
            print "k", k, ":", l, "=?=", self.eigenvalues[k], "normalized |diff|", np.abs((l-self.eigenvalues[k])/l)



    def run(self):
        self.mst = self._make_minimum_spanning_tree()
        mst = self.mst

        # the barrier function gives the minimax barrier separating
        # a node from a sink node
        barrier_function = dict()
        # the escape function is usually barrier_function[i] - E[i]
        escape_function = dict() # escape function

        self.log_eigenvalues = np.zeros(self.nnodes) # used to compute eigenvalues
        self.eigenvalues = np.zeros(self.nnodes) # used to compute eigenvalues
        self.eigenvectors = np.zeros([self.nnodes, self.nnodes])
        self.eigenvectors[:,0] = 1.
        self.eigenvector_indicator = np.zeros(self.eigenvectors.shape)
        self.log_eigenvectors = np.zeros(self.eigenvectors.shape)
        
        k = 0
        
        # the first sink is the node with the lowest energy
        s1 = min((((E, n) for n, E in self.Ei.iteritems())))[1]
        sinks = set([s1])
        self.slist = [s1]
        self.cut_edges = []
        
        # compute the initial values of the barrier function and escape function
        self._compute_barrier_function(mst, s1, barrier_function, escape_function)

        for k in xrange(1, len(self.Ei)):
            # the sink is the node with the maximum value of escape function
            s = max(((v, i) for i, v in escape_function.iteritems()))[1]
            
            if self.verbose:
                print "\niteration", k
                print "current sink", s
                print "past sinks", self.slist
                print "total # edges", mst.number_of_edges()
                print "size of component connected to current sink", len(nx.node_connected_component(mst, s))
            
                for i, ui in barrier_function.iteritems(): 
                    print "  barrier_function[",i,"] = ", ui
            
            # find the new cutting edge
            Epq, p, q = self._find_cutting_edge(mst, s, barrier_function, sinks)
            if self.verbose:
                s2 = self._get_connected_sink(mst, s, sinks)
                print "p q Epq", p, q, Epq
                print "s Es   ", s, self.Ei[s], "sold Esold", s2, self.Ei[s2] 
            
            # save the eigenvalue
            self.log_eigenvalues[k] = -(Epq - self.Ei[s]) / self.T
            self.eigenvalues[k] = - np.exp(self.log_eigenvalues[k])
            
            # remove the cutting edge
            mst.remove_edge(p, q)
            if self.verbose:
                print "removing edge (",p,",",q,")"
            
            # compute the eigenvector
            self._make_kth_eigenvector(mst, k, p, q, s)
            
            # recompute barrier_function and escape_function
            self._recompute_barrier_function(mst, s, barrier_function, escape_function)
            
            sinks.add(s)
            self.slist.append(s)
            self.cut_edges.append((p,q))

        if self.verbose:
            self._test_eigenvectors()
            
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

    def logaddexp(self, v1, v2):
        if np.isnan(v1):
            return v2
        else:
            return np.logaddexp(v1, v2)
        

    def make_K_inv(self):
        log_kinv_plus = np.zeros([self.spect.nnodes, self.spect.nnodes])
        log_kinv_minus = np.zeros([self.spect.nnodes, self.spect.nnodes])
        levec = self.spect.log_eigenvectors
        ind = self.spect.eigenvector_indicator
        leval = self.spect.log_eigenvalues
        for i in xrange(self.spect.nnodes):
            for j in xrange(self.spect.nnodes):
                for k in xrange(1, self.spect.nnodes):
                    indicator = -1 * ind[i,k] * ind[j,k]
                    if indicator == 0:
                        continue
                    val = - leval[k] + levec[i,k] + levec[j,k] + self.spect.log_pi[j]
                    if indicator > 0:
                        log_kinv_plus[i,j] = self.logaddexp(log_kinv_plus[i,j], val)
                    else:
                        log_kinv_minus[i,j] = self.logaddexp(log_kinv_minus[i,j], val)
        
        # now merge the positive and negative components of log_kinv
        max_vals = np.where(log_kinv_plus > log_kinv_minus, log_kinv_plus, log_kinv_minus)
        kinv = np.exp(log_kinv_plus - max_vals) - np.exp(log_kinv_minus - max_vals)
        kinv *= np.exp(max_vals)
        
        return kinv
                
    
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
            mat_inv += 1./evals[k] * np.outer(v, v * pi)
        print "approximate K matrix"
        print_matrix(mat)
        print "approximate inv K"
        print_matrix(mat_inv)
        
        if True:
            mat_inv_alt = self.make_K_inv()
            print "approximate inv K make from log sums"
            print_matrix(mat_inv_alt)
            mat_inv = mat_inv_alt
        

        if True:       
            m, eval, evec = get_eigs(self.spect.Ei, self.spect.Eij, T=self.spect.T)
            print "exact K matrix"
            print_matrix(m)
            print "condition number of m", np.linalg.cond(m[:-1, :-1])
        


        Kcond = np.dot(mat_inv, m)
        print "conditioned matrix"
        print_matrix(Kcond)
        print "reduced conditioned matrix"
        print_matrix(np.dot(mat_inv[1:,1:], m[1:,1:]))
        print "condition number of conditioned matrix", np.linalg.cond(Kcond[:-1,:-1])
       
        viewer = ViewMSTSpectralDecomp(self.spect)

 
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
    
    viewer = ViewMSTSpectralDecomp(spect)

        
#        print m


def test2():
    np.random.seed(1)
    Ei, Eij = make_random_energies_complete(15)
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
    
    viewer = ViewMSTSpectralDecomp(spect)

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


def test_precond2(n=8, T=.01):
    np.random.seed(5)
    Ei, Eij = make_random_energies_complete(n)
    precond = MSTPreconditioning(Ei, Eij, T=T)
    
    


if __name__ == "__main__":
    from tests.test_preconditioning import make_random_energies_complete, get_eigs

#     test1()
    test_precond2(n=5, T=.05)
