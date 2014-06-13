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
import mpmath

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

class ViewMSTSpectralDecompDgraph(object):
    def __init__(self, spect):
        self.spect = spect
        self.run()

    def draw_nodes(self, nodes, color='k'):
        xypos = [self.pos[n] for n in nodes]
        x = [xy[0] for xy in xypos]
        y = [xy[1] for xy in xypos]
        plt.scatter(x, y, c=color, linewidths=0, s=60)

    def draw_edges_old(self, edges, color='k', lw=1.):
        for u, v in edges:
            xy = np.array([self.pos[u], self.pos[v]])
            plt.plot(xy[:,0], xy[:,1], color, lw=lw)
    
    def draw_edges(self, edges, color='k', lw=1.):
        from matplotlib.collections import LineCollection
        line_segments = []
        line_colors = []
        eoffset = .1
#        line_segments, line_colors = self.dgraph._get_line_segments(self.dgraph.tree_graph, eoffset=.1)
        for u, v in edges:
            if self.spect.Ei[u] < self.spect.Ei[v]:
                u, v = v, u
#            tree = self.dgraph_get_ts_tree(u, v)
            trees = [self.dgraph_get_leaf(u)]#, self.dgraph_get_leaf(v)]
#            parents = [t.parent for t in trees]
#            trees += parents
            for tree in trees:
                self.dgraph._get_line_segment_single(line_segments, line_colors, tree, eoffset)
#            tree = self.dgraph_get_leaf(u)
#            tree = self.dgraph_get_leaf(v)
#            self.dgraph._get_line_segment_single(line_segments, line_colors, tree, eoffset)
#            parents = tree.parent
#            self.dgraph._get_line_segment_single(line_segments, line_colors, tree, eoffset)
#        linecollection = LineCollection([ [(x[0],y[0]), (x[1],y[1])] for x,y in line_segments])
#        linecollection.set_linewidth(lw)
#        linecollection.set_color(color)
#        ax = plt.gca()
#        ax.add_collection(linecollection)

        for x, y in line_segments:
            plt.plot(x, y, color, lw=lw)
            

    
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
    
    def draw_old(self, k):
        m=k+1
        groups = [self.spect.slist[:m]]
        groups.append([self.spect.slist[m]])
        colors = ["red", "black"]
            
#        self.dgraph.color_by_group(groups, colors=colors)
#        self.dgraph.plot(show_minima=True)
#        self.dgraph.show()    

    
    def dgraph_get_leaf(self, nodeid):
        for leaf in self.dgraph.tree_graph.leaf_iterator():
            if leaf.data["minimum"].nodeid == nodeid:
                return leaf
    
    def dgraph_get_ts_tree(self, u, v):
        from pele.utils.disconnectivity_graph import TreeLeastCommonAncestor
        leafu = self.dgraph_get_leaf(u)
        leafv = self.dgraph_get_leaf(v)
        lca = TreeLeastCommonAncestor([leafu, leafv])
        return lca.least_common_ancestor
    
    def make_dgraph(self, mst):
        from pele.storage import Database
        from pele.utils.disconnectivity_graph import DisconnectivityGraph, database2graph
        db = Database()
        self.node2minimum = dict()
        for u, v, data in mst.edges_iter(data=True):
            Ets = data["energy"]
            mu = db.addMinimum(self.spect.Ei[u], [0])
            mv = db.addMinimum(self.spect.Ei[v], [0])
            mu.nodeid = u
            mv.nodeid = v
            self.node2minimum[u] = mu
            self.node2minimum[v] = mv
            db.addTransitionState(Ets, [0], mu, mv)
        
        
        dgraph = DisconnectivityGraph(database2graph(db), nlevels=self.spect.nnodes+10)
        dgraph.calculate()
        
        line_segments = dgraph._get_line_segments(dgraph.tree_graph)
        
        # get the minima layout
        xpos, minima = dgraph.get_minima_layout()
        ypos = [m.energy for m in minima]
        self.pos = dict(((m.nodeid, (x, y)) for m, x, y in izip(minima, xpos, ypos)))
        print "positions", self.pos
        
        # get the transition states
        return dgraph
        
    
    def run(self):
        self.mst = self.spect._make_minimum_spanning_tree()
        self.dgraph = self.make_dgraph(self.mst)
#        self.pos = nx.spring_layout(self.mst)
        self.removed_edges = []
        
        for k in xrange(len(self.spect.Ei)-1):
            self.draw(k)
            
#            self.mst.remove_edge(*self.spect.cut_edges[k])
            self.removed_edges = self.spect.cut_edges[:(k+1)]

def _normalize_eigenvectors(eigenvectors, pi):
    for k in xrange(eigenvectors.shape[0]):
        v = eigenvectors[:,k]
        eigenvectors[:,k] /= np.sqrt(np.dot(v, v * pi))


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

class _EigenvectorRepresentation(object):
    def __init__(self, group1, group2, log_pi, mppi):
        self.group1 = np.array(group1)
        self.group1.sort()
        self.group2 = np.array(group2)
        self.group2.sort()
        self.nnodes = len(log_pi)

        self.log_pi = log_pi
        self.mppi = mppi
        
        isum_pi1 = - scipy.misc.logsumexp(self.log_pi[self.group1])
        isum_pi2 = - scipy.misc.logsumexp(self.log_pi[self.group2])
        log_norm = np.logaddexp(isum_pi1, isum_pi2) / 2
        self.logC1 = isum_pi1 - log_norm
        self.logC2 = isum_pi2 - log_norm
        
        self.log_terms1 = self.log_pi[self.group1] - log_norm
        self.log_terms2 = self.log_pi[self.group2] - log_norm
        self.indicator = np.zeros(self.nnodes)
        self.indicator[self.group1] = 1
        self.indicator[self.group2] = -1

    
    def make_log_eigenvector(self):
        log_evec = np.zeros(self.nnodes)
        log_evec[self.group1] = self.logC1
        log_evec[self.group2] = self.logC2
        return log_evec, self.indicator
    
    def make_eigenvector(self):
        evec = np.zeros(self.nnodes)
        evec[self.group1] = np.exp(self.logC1)
        evec[self.group2] = -np.exp(self.logC2)
        return evec
    
    def get_log_terms(self, index):
        ind = self.indicator[index]
        if index == 0:
            return np.array([0]), 0
        elif index > 0:
            return self.log_terms1, ind
        else:
            return self.log_terms2, ind
    
    def make_mp_eigenvector(self):
        isum_pi1 = mpmath.mpf(1) / sum((self.mppi[i] for i in self.group1))
        isum_pi2 = mpmath.mpf(1) / sum((self.mppi[i] for i in self.group2))
        norm = mpmath.sqrt(isum_pi1 + isum_pi2)
        C1 = isum_pi1 / norm
        C2 = isum_pi2 / norm
        
        v = mpmath.matrix([0.]*self.nnodes)
        for i in self.group1:
            v[i] = C1
        for i in self.group2:
            v[i] = -C2
        return v
            
        

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

    def _make_equilibrium_occupation_probabilities_mp(self):
        self.mppi = mpmath.matrix([mpmath.exp(mpmath.mpf(-self.Ei[node]) / self.T) 
                         for node in self.node_list])
        self.mppi /= sum(self.mppi)
        print "mp pi"#, mp.dps, mpmath
        print self.mppi
        
    def _make_equilibrium_occupation_probabilities(self):
        self._make_equilibrium_occupation_probabilities_mp()
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
            self.graph.add_node(i, energy=E)
        for (i, j), E in self.Eij.iteritems():
            self.graph.add_edge(i, j, energy=E)
        mst = nx.minimum_spanning_tree(self.graph, weight='energy')
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
        
        is1 = np.array([self.node2index[node] for node in S1])
        is2 = np.array([self.node2index[node] for node in S2])
        evec_repr = _EigenvectorRepresentation(is1, is2, self.log_pi, self.mppi)
        
        
        # compute the inverse sum of pi for all nodes in S1 and S2
#         isum_pi1 = - scipy.misc.logsumexp([self.log_pi_dict[node] for node in S1])
#         isum_pi2 = - scipy.misc.logsumexp([self.log_pi_dict[node] for node in S2])
#         norm = scipy.misc.logsumexp([isum_pi1, isum_pi2])
#         logC1 = isum_pi1 - norm / 2.
#         logC2 = isum_pi2 - norm / 2.
        
        self.eigenvectors[:,k] = evec_repr.make_eigenvector()
#         self.eigenvectors[is1, k] =  np.exp(logC1)
#         self.eigenvectors[is2, k] = -np.exp(logC2)
        
        self.log_eigenvectors[:,k], self.eigenvector_indicator[:,k] = evec_repr.make_log_eigenvector()
        self.eigenvector_full_representations[k] = (evec_repr)
        self.eigevectors_mp[:,k] = evec_repr.make_mp_eigenvector()
#         self.log_eigenvectors[is1, k] = logC1
#         self.log_eigenvectors[is2, k] = logC2
        
#         self.eigenvector_indicator[:,k] = 0
#         self.eigenvector_indicator[is1,k] =  1
#         self.eigenvector_indicator[is2,k] = -1

    def _test_eigenvectors(self, evals=None, evecs=None, m=None, pi=None):
        if evals is None:
            evals = self.eigenvalues
        if evecs is None:
            evecs = self.eigenvectors
        if pi is None:
            pi = self.pi
        if m is None:
            from tests.test_preconditioning import make_rate_matrix
            m = make_rate_matrix(self.Ei, self.Eij, T=self.T)
        # test to see if the eigenvectors are orthonormal
        print "\ntesting to see if the eigenvectors are normalized sum(v*v*pi) == 1"
        for k in xrange(len(self.Ei)):
            v = evecs[:,k]
            r = np.sum(v * v * pi)
            if abs(r-1) > .1:
                print "k", k, ":", r
        print "testing to see if the eigenvectors are orthogonal sum(u*v*pi) == 0"
        for k1 in xrange(len(self.Ei)):
            for k2 in xrange(k1):
                v1 = evecs[:,k1]
                v2 = evecs[:,k2]
                r = np.sum(v1 * v2 * pi)
                if np.abs(r) > .01:
                    print "k1 k2", k1, k2, ":", r

        print "testing to see if the the eigenvectors satisfy the eigenvalue equation dot(v*pi, dot(m, v)) = lambda"
        for k in xrange(len(evals)):
            v = evecs[:,k]
            l = np.dot(v*pi, np.dot(m, v))
            print "k", k, ":", l, "=?=", evals[k], "normalized |diff|", np.abs((l-evals[k])/l)

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
        self.eigenvector_full_representations = [None] * self.nnodes
        self.eigevectors_mp = mpmath.matrix(self.eigenvectors)

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
            print "eigenvalue", self.eigenvalues[k], Epq, self.Ei[s] 
            
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
            print ""
            print "eigenvectors (mpmath)"
            print self.eigevectors_mp

def subtract_exp(pos, neg):
    if np.isnan(pos):
        p = 0.
    else:
        p = np.exp(pos)
    if np.isnan(neg):
        n = 0.
    else:
        n = np.exp(neg)
    return p - n
    
def logaddexp(v1, v2):
    if np.isnan(v1):
        return v2
    else:
        return np.logaddexp(v1, v2)

# class LogAccumulator(object):
#     log_pos = np.NAN
#     log_neg = np.NAN
#     
#     def insert(self, log_val, indicator):
#         if indicator == 0:
#             return
#         if indicator > 0:
#             self.log_pos = logaddexp(self.log_pos, log_val)
#         else:
#             self.log_neg = logaddexp(self.log_neg, log_val)
#     
#     def get_result(self):
#         if np.isnan(self.log_pos):
#             p = 0.
#         else:
#             p = np.exp(self.log_pos)
#         if np.isnan(self.log_neg):
#             n = 0.
#         else:
#             n = np.exp(self.log_neg)
#         return p - n
from _precision_summation import MSum
LogAccumulator = MSum

            

class MSTPreconditioning(object):
    def __init__(self, Ei, Eij, sinks, T=0.1, verbose=True):
        self.verbose = verbose
        self.sinks = set(sinks)
        assert len(sinks) == 1
        Ei = self.prepare_Ei(Ei, sinks)
        self.spect = MSTSpectralDecomposition(Ei, Eij, T=T)
        self.make_log_rate_matrix()
    
    def prepare_Ei(self, Ei, sinks):
        self.Ei_orig = Ei.copy()
        Emin = min(Ei.itervalues())
        for s in sinks:
            Ei[s] = Emin - 20
        return Ei

        

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
                        log_kinv_plus[i,j] = logaddexp(log_kinv_plus[i,j], val)
                    else:
                        log_kinv_minus[i,j] = logaddexp(log_kinv_minus[i,j], val)
        
        # now merge the positive and negative components of log_kinv
        max_vals = np.where(log_kinv_plus > log_kinv_minus, log_kinv_plus, log_kinv_minus)
        kinv = np.exp(log_kinv_plus - max_vals) - np.exp(log_kinv_minus - max_vals)
        kinv *= np.exp(max_vals)
        
        return kinv
    
    def make_log_rate_matrix(self):
        indicator = np.zeros(self.spect.eigenvectors.shape)
        logK = np.zeros(indicator.shape)
        K = np.zeros(logK.shape)
        for (i,j), Eij in self.spect.Eij.iteritems():
            logK[i,j] = -(Eij - self.spect.Ei[i]) / self.spect.T
            logK[j,i] = -(Eij - self.spect.Ei[j]) / self.spect.T
            indicator[i,j] = 1
            indicator[j,i] = 1
            K[i,j] = np.exp(logK[i,j])
            K[j,i] = np.exp(logK[j,i])
        for i in xrange(self.spect.nnodes):
            ind = np.where(indicator[i,:]>0)
            log_sumK = scipy.misc.logsumexp(logK[i,ind])
            indicator[i,i] = -1
            logK[i,i] = log_sumK
            K[i,i] = -np.exp(logK[i,i])
        self.logK = logK
        self.logK_indicator = indicator
        self.K = K



    def make_preconditioned_submatrix(self, iremove=None):
        if iremove is None:
            assert len(self.sinks) == 1
            iremove = iter(self.sinks).next()
        ilist = range(iremove) + range(iremove+1,self.spect.nnodes)
        assert len(ilist) == self.spect.nnodes - 1
        nnew = self.spect.nnodes - 1
        M = np.zeros([nnew, nnew])
        levecs = self.spect.log_eigenvectors
        ind = self.spect.eigenvector_indicator
        levals = self.spect.log_eigenvalues

        for inew in xrange(nnew): # index i
            i = ilist[inew]
            for jnew in xrange(nnew): # index j
                j = ilist[jnew]
                acc = LogAccumulator()
                for n in xrange(1,self.spect.nnodes): # eigenvalue number
                    for k in xrange(self.spect.nnodes): # matrix multiplication index
                        if k == iremove:
                            continue
                        indicator = -1 * ind[i,n] * ind[k,n] * self.logK_indicator[k,j]
                        logval = (- levals[n]
                                  + levecs[i,n]
                                  + levecs[k,n]
                                  + self.spect.log_pi[k]
                                  + self.logK[k,j]
                                  )
                        acc.insert(logval, indicator)
                M[inew,jnew] = acc.get_result()
#        print M
        return M

    def make_mfpt_right_vector(self):
        """the pseudo inverse of K * ones()"""
        ilist = sorted(set(xrange(self.spect.nnodes)) - self.sinks)
        nnew = len(ilist)
        right = np.zeros(nnew)
        ind = self.spect.eigenvector_indicator
        levals = self.spect.log_eigenvalues
        levecs = self.spect.log_eigenvectors
        for inew in xrange(nnew):
            i = ilist[inew]
            sum_logval_pos = np.nan
            sum_logval_neg = np.nan
            for n in xrange(1,self.spect.nnodes):
                for k in ilist:
                    indicator = -1 * ind[i,n] * ind[k,n]
                    if indicator == 0:
                        continue
                    logval = (- levals[n]
                              + levecs[i,n]
                              + levecs[k,n]
                              + self.spect.log_pi[k]
                              )
                    if indicator > 0:
                        sum_logval_pos = logaddexp(sum_logval_pos, logval)
                    else:
                        sum_logval_neg = logaddexp(sum_logval_neg, logval)
            if np.isnan(sum_logval_pos):
                vp = 0
            else:
                vp = np.exp(sum_logval_pos)
            if np.isnan(sum_logval_neg):
                vn = 0
            else:
                vn = np.exp(sum_logval_neg)
#                print sum_logval_pos, sum_logval_neg, vp, vn
            right[inew] = vp - vn
        return right
                    
        
    
    def make_submatrix(self, K, i):
        n = K.shape[0]
        M = np.zeros([n-1,n-1])
        M[:i,:i] = K[:i,:i]
        M[i:,i:] = K[i+1:,i+1:]
        return M

    def compute_mfp_times(self):
        M = self.make_preconditioned_submatrix()
        b = self.make_mfpt_right_vector()
        
        print "\nconditioned matrix"
        print_matrix(M)
        print "condition number of conditioned matrix", np.linalg.cond(M)
        times = np.linalg.solve(M, b)
        print "mean first passage times"
        print_matrix(times)
    
    def test_eigenvector_orthonormality(self, n, m, iremove):
        ind = self.spect.eigenvector_indicator
        levecs = self.spect.log_eigenvectors
        lpi = self.spect.log_pi
        acc = LogAccumulator()
        for i in xrange(self.spect.nnodes):
            if i == iremove:
                continue
            indicator = ind[i,n] * ind[i,m]
            val = levecs[i,n] + levecs[i,m] + lpi[i]
            acc.insert(val, indicator)
        return acc.get_result()

    def test_eigenvector_equation(self, n, m, iremove):
        ind = self.spect.eigenvector_indicator
        levecs = self.spect.log_eigenvectors
        lpi = self.spect.log_pi
        acc = LogAccumulator()
        from _precision_summation import MSum
        acc = MSum()
        for i in xrange(self.spect.nnodes):
            if i == iremove:
                continue
            for j in xrange(self.spect.nnodes):
                if j == iremove: 
                    continue
                indicator = ind[i,n] * ind[j,m] * self.logK_indicator[i,j]
                val = levecs[i,n] + lpi[i] + self.logK[i,j] + levecs[j,m]
                acc.insert(val, indicator)
                if n == m == 1:
                    pass
        return acc.get_result()

    
    def test_submatrix_eigenvectors(self):
        assert len(self.sinks) == 1
        iremove = iter(self.sinks).next()
#        nnew = self.spect.nnodes - 1
#        m = self.K.copy()
#        Kii = m[iremove,iremove]
#        phi_ii = self.spect.eigenvectors[iremove,1]
#        pi_ii = self.spect.pi[iremove]
#        
#        m = self.make_submatrix(m, iremove)
#        evecs = np.zeros([self.spect.nnodes-1, self.spect.nnodes])
#        evecs[:iremove,:] = self.spect.eigenvectors[:iremove,:]
#        evecs[iremove:,:] = self.spect.eigenvectors[iremove+1:,:]
#        evals = self.spect.eigenvalues.copy()
#        pi = np.zeros(self.spect.nnodes-1)
#        pi[:iremove] = self.spect.pi[:iremove]
#        pi[iremove:] = self.spect.pi[iremove+1:]
#        
#        norm = np.sqrt(1. - phi_ii**2 * pi_ii)
#        print "normalizing eigenvector 1 by", norm
#        print "actual norm of eigenvector 1 is", np.sum(evecs[:,1]**2 * pi)
#        evecs[:,1] /= norm
        
#        self.spect._test_eigenvectors(evals=evals, evecs=evecs, m=m, pi=pi)
#        print "norm of eigenvector 1 should be", 1. - self.spect.eigenvectors[iremove,1]**2 * self.spect.pi[iremove]
        
#        print "full eigenvector 1 dot full eigenvector 2"
#        print np.sum(self.spect.eigenvectors[:,1] * self.spect.eigenvectors[:,2] * self.spect.pi)
#        print "reduced eigenvector 1 dot reduced eigenvector 2"
#        print np.sum(evecs[:,1] * evecs[:,2] * pi)
        
        print "\ntesting eigenvectors are normal"
        for n1 in xrange(1,self.spect.nnodes):
            r = self.test_eigenvector_orthonormality(n1, n1, iremove)
            if np.abs(r-1) > .1:
                print "  ", n1, ":", r 

        print "testing eigenvectors are orthogonal"
        for n1 in xrange(1,self.spect.nnodes):
            for m1 in xrange(1,self.spect.nnodes):
                if n1 == m1:
                    continue
                r = self.test_eigenvector_orthonormality(n1, m1, iremove)
                if np.abs(r) > 1e-2:
                    print "  ", n1, m1, ":", r
        print "testing eigenvector equation (n==m)"
        for n1 in xrange(1,self.spect.nnodes):
            lam = self.test_eigenvector_equation(n1, n1, iremove)
            lam_full = self.test_eigenvector_equation(n1, n1, -1)
            print "  ", n1, ":", lam, "expected", self.spect.eigenvalues[n1], "full", lam_full
        print "testing eigenvector equation (n!=m)"
        for n1 in xrange(1,self.spect.nnodes):
            for m1 in xrange(1,self.spect.nnodes):
                if n1 == m1:
                    continue
                r = self.test_eigenvector_equation(n1, m1, iremove)
                if np.abs(r) > 1e-2:
                    print "  ", n1, m1, ":", r
        pass   
        
#        print "\ntesting eigenvector equation with eigenvector 1"
#        print "K * evec[:,1]"
#        v = np.dot(m, evecs[:,1])
#        print_matrix(v)
#        print "eval[1] * evec[:,1]"
#        print_matrix(evals[1] * evecs[:,1])
#        eval1new = norm**-2 *(evals[1] - phi_ii**2 * pi[iremove] * (2*evals[1] - Kii))
#        print "using new estimate of eigenvalue", eval1new, "(compared with", evals[1],")"
#        print_matrix(eval1new * evecs[:,1])
#        print "is v=dot(m, evecs[:,1]) orthogonal to all other eigenvectors?"
#        print "first normalize v", np.sum(v**2 * pi)
#        v /= np.sum(v**2 * pi)
#        for k in xrange(1,self.spect.nnodes):
#            print "  ", k, ":", np.dot(v * pi, evecs[:,k])
##        print "now with respect to pi_new: first normalize v", np.sum(v**2 * pi)
##        v /= np.sum(v**2 * pi)
##        for k in xrange(1,self.spect.nnodes):
##            print "  ", k, ":", np.dot(v * pi, evecs[:,k])
#        
#        
#        print "\neigenvalues of submatrix K"
#        print_matrix(np.linalg.eigvals(m))
#        print "eigenvalues of full K"
#        print_matrix(self.spect.eigenvalues)
#        print "eigenvalues of full K"
#        print_matrix(get_eigs(self.spect.Ei, self.spect.Eij, self.spect.T)[1])
#        print ""
    
    def test_result(self, K, Kinv, iremove):
        i = iremove
        Ksub = self.make_submatrix(K, iremove)
        Kinv_sub = self.make_submatrix(Kinv, iremove)
        red = Kinv_sub.dot(Ksub)
        cond = np.linalg.cond(red)
        Ksub_cond = self.make_preconditioned_submatrix(iremove)
        cond_approx = np.linalg.cond(Ksub_cond)
        if self.spect.pi[iremove] > .1 or cond < 1e5 or cond_approx < 1e5:
            print "\nremoving the row and column", iremove
            print "pi of this node", self.spect.pi[iremove]
            print "(submatrix(pseudo inverse K) * submatrix(K)) (exact)"
            print_matrix(red)
            print "condition number of unconditioned submatrix", np.linalg.cond(Ksub)
            print "eigenvalues of submatrix(K)"
            print_matrix(np.linalg.eigvals(Ksub))
            print "eigenvalues of (submatrix(pseudo inverse K) * submatrix(K)) (exact)"
            print_matrix(np.linalg.eigvals(red))
            print "singular values of (submatrix(pseudo inverse K) * submatrix(K)) (exact)"
            print_matrix(np.linalg.svd(red)[1])
            print "condition number:", cond
            print "******"
            print "(submatrix(pseudo inverse K) * submatrix(K)) (approx)"
            print_matrix(Ksub_cond)
            print "eigenvalues of (submatrix(pseudo inverse K) * submatrix(K)) (approx)"
            print_matrix(np.linalg.eigvals(Ksub_cond))
            print "singular values of (submatrix(pseudo inverse K) * submatrix(K)) (approx)"
            print_matrix(np.linalg.svd(Ksub_cond)[1])
            print "condition number:", cond_approx
    
    def run_tests(self):
        n = len(self.spect.Ei)
        evals = self.spect.eigenvalues
        evecs = self.spect.eigenvectors
        mat = np.zeros([n, n])
        mat_inv = np.zeros([n, n])
        pi = self.spect.pi
        
        if self.verbose:
            m, xevals, xevecs = get_eigs(self.spect.Ei, self.spect.Eij, T=self.spect.T)
            xmat = np.zeros(mat.shape)
            xmat_inv = np.zeros(mat.shape)
            _normalize_eigenvectors(xevecs, pi)
            for k in xrange(1,n):
                v = xevecs[:,k]
                xmat += xevals[k] * np.outer(v, v * pi)
                xmat_inv += 1./xevals[k] * np.outer(v, v * pi)

        
#        for k in xrange(1,n):
#            v = evecs[:,k]
#            mat += evals[k] * np.outer(v, v * pi)
#            mat_inv += 1./evals[k] * np.outer(v, v * pi)

        if True:
            mat = self.K
            mat_inv = self.make_K_inv()
        
        if self.verbose:
            print ""
            print "approximate K matrix"
            print_matrix(mat)
            print "K (from exact eigenvectors)"
            print_matrix(xmat)
            print "exact K matrix"
            print_matrix(m)
            print ""
            print "pseudo inverse K (from approximate eigenvectors)"
            print_matrix(mat_inv)
#            print "pseudo inverse K (from log sums)"
#            print_matrix(mat_inv_alt)
            print "pseudo inverse K (from exact eigenvectors)"
            print_matrix(xmat_inv)
            
        
#        mat_inv = mat_inv_alt
        

#        print "pseudo inverse K * K (from exact eigenvectors)"
#        print_matrix(xmat_inv * m)
#        Kcond = np.dot(mat_inv, m)
#        print "conditioned matrix (K_inv * K)"
#        print_matrix(Kcond)
        if True:
            print "\n(pseudo inverse K) * K (approx)"
            print_matrix(mat_inv.dot(m))
            print "(pseudo inverse K) * K (exact)"
            print_matrix(xmat_inv.dot(m))
            print "\nremoving the first row and column"
            print "submatrix(pseudo inverse K) * submatrix(K) (exact)"
            sl = slice(1,None)
            print_matrix(np.dot(xmat_inv[sl,sl], m[sl,sl]))
            print "eigenvalues of submatrix(K)"
            print_matrix(np.linalg.eigvals(m[sl,sl]))
            print "eigenvalues of (submatrix(pseudo inverse K) * submatrix(K)) (exact)"
            print_matrix(np.linalg.eigvals(xmat_inv[sl,sl].dot(m[sl,sl])))
            print "condition number:", np.linalg.cond(xmat_inv[sl,sl].dot(m[sl,sl]))
        
        if False:
            print "\nnow removing the last row and column"
            sl = slice(None,-1)
            print "eigenvalues of submatrix(K)"
            print_matrix(np.linalg.eigvals(m[sl,sl]))
            print "eigenvalues of (submatrix(pseudo inverse K) * submatrix(K)) (exact)"
            print_matrix(np.linalg.eigvals(xmat_inv[sl,sl].dot(m[sl,sl])))
            print "condition number:", np.linalg.cond(xmat_inv[sl,sl].dot(m[sl,sl]))
            print ""

        for i in xrange(self.spect.nnodes):
            self.test_result(m, xmat_inv, i)
        
        if False:
            Kcond = np.dot(mat_inv[:-1,:-1], m[:-1,:-1])
            Kcond_full = np.dot(mat_inv, m)
    #        print_matrix(Kcond)
            if self.verbose:
                print "condition number of submatrix", np.linalg.cond(m[:-1,:-1])
            print "condition number of conditioned submatrix", np.linalg.cond(Kcond)
            print "condition number of conditioned submatrix", np.linalg.cond(Kcond_full[:-1,:-1])
       
        if False:
            print "\ninvert the pseudo inverse"
            print_matrix(np.linalg.inv(xmat_inv))
            print "rank 1 update: remove the last row and column"
            A = m[:-1,:-1]
            print "the eigenvalues of the sub matrix"
            print np.linalg.eigvals(A)
            print np.linalg.eigvals(m[1:,1:])
            Ainv = xmat_inv[:-1,:-1] - np.outer(m[:-1,-1], m[-1,:-1]) / m[-1,-1]
            print_matrix(A)
            print_matrix(Ainv)
            print_matrix(A.dot(Ainv))
            
            
        
        if False:
            print "\ntest if it is a pseudo inverse"
            print_matrix(mat.dot(mat_inv).dot(mat))
            print_matrix(mat)
            print "\nare they the same?"
            print_matrix(mat_inv.dot(mat).dot(mat_inv))
            print_matrix(mat_inv)
            print "\nis symmetric?"
            print_matrix(mat.dot(mat_inv))
            print "\nis symmetric?"
            print_matrix(mat_inv.dot(mat))
        
        self.test_submatrix_eigenvectors()

#         viewer = ViewMSTSpectralDecomp(self.spect)
        
 
#
# only testing below here
#


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
    seed = 5699024 #n=6
    seed = 1492541 #n=5
    seed = 8833179 #n=4
    seed = 10
#     seed = np.random.randint(10000000)
    print "seed", seed
    np.random.seed(seed)
    Ei, Eij = make_random_energies_complete(n)
    print "energies", Ei
    precond = MSTPreconditioning(Ei, Eij, [1], T=T)   
    precond.run_tests()
    precond.compute_mfp_times()
    

if __name__ == "__main__":
    from tests.test_preconditioning import make_random_energies_complete, get_eigs
    mpmath.mp.dps = 20
#     test1()
    test_precond2(n=9, T=.1)
