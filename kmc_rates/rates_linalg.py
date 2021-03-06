import time
from itertools import izip
from collections import defaultdict

import numpy as np
import scipy.sparse
import scipy.sparse.linalg
import networkx as nx

class LinalgError(Exception):
    """raised when a the linear algebra solver fails"""

def reduce_rates(rates, B, A=None):
    B = set(B)
    if A is not None:
        A = set(A)
        if A.intersection(B):
            raise Exception("A and B share", len(A.intersection(B)), "nodes")
    graph = nx.Graph()
    graph.add_edges_from(rates.iterkeys())
    
    # remove nodes not connected to B
    # TODO: this only works if B is fully connected
    connected_nodes = nx.node_connected_component(graph, iter(B).next())
    connected_nodes = set(connected_nodes)
    all_nodes = set(graph.nodes())
    if len(connected_nodes) != len(all_nodes):
        print "removing", len(all_nodes) - len(connected_nodes), "nodes that are not connected to B"
    
        rates = dict((uv, rate) for uv, rate in rates.iteritems()
                          if uv[0] in connected_nodes
                          )
        
        if B - connected_nodes:
            raise Exception("the nodes in B are not all connected")
        
        if A is not None:
            if A - connected_nodes:
                raise Exception("the A nodes are not all connected to the B nodes")

    return rates

def compute_sum_out_rates(rates):
    rates_list = defaultdict(list)
    for uv, rate in rates.iteritems():
        rates_list[uv[0]].append(rate)
    
    #sum rates more precisely
#    print "recomputing the sum of the rates more precisely"
    sum_out_rates = dict()
    for u, urates in rates_list.iteritems():
        urates.sort()
        sumrate = sum(urates)
        if False:
            import decimal
            urates_dec = map(decimal.Decimal, urates)
            sumrate = float(sum(urates_dec))
        sum_out_rates[u] = sumrate
    return sum_out_rates


class EstimateRates(object):
    """make a rough estimate of the rates at very low temperatures
    
    This is aimed to be used as an initial guess for the rate solver.
    Equivalently it can be used to precondition the matrix in the system of
    linear equations
    
    Notes
    -----
    The idea is that at very low temperatures a trajectory will first roll down hill
    to the lowest energy minimum it can get to.  The folding rate will be given by 
    the time to get from the lowest energy down hill minimum to the target structure.
    """
    def __init__(self, rate_constants, Peq, B):
        self.union_find = nx.utils.UnionFind()
        self.rate_constants = rate_constants
        self.Peq = Peq
        self.B = set(iter(B))
        
        self.run()
    
    def minimum_spaning_edges(self, edges):
        # assume the edges are already sorted appropriately
        subtrees = self.union_find
        for u,v in edges:
            uroot, vroot = subtrees[u], subtrees[v]
            if uroot != vroot:
                yield u, v, uroot, vroot
                subtrees.union(u,v)

    
    def run(self):
        # add B to the union find and make sure they're all connected
        b = iter(self.B).next()
        for x in self.B:
            if x != b:
                self.union_find.union(x, b)
        
        # sort edges with smallest free energy to the left.
        # the free energy is proportional to k_uv * Peq[u]
        edges = [((u,v), k*self.Peq[u]) for (u,v), k in 
                 self.rate_constants.iteritems() if u<v]
        edges.sort(key=lambda uvk: -uvk[1])
        assert edges[0][1] > edges[1][1]
        edges = [uv for uv, k in edges]
        
        rate_estimates = dict()
        graph = nx.Graph()
        for u, v, uroot, vroot in self.minimum_spaning_edges(edges):
            broot = self.union_find[b]
            
            if uroot == broot or vroot == broot:
                if uroot == broot:
                    new_node = v
                if vroot == broot:
                    new_node = u
                
                # all nodes connected to new_node are newly connected to the B nodes
                if new_node not in graph:
                    graph.add_node(new_node)
                nodes = nx.node_connected_component(graph, new_node)
                # The edge u, v is the transition state with the smallest free energy 
                # that connects these nodes to B.  This transition state is the bottleneck  
                ts_free_energy = self.rate_constants[(u,v)] * self.Peq[u]
#                print len(nodes), "nodes connecting to B"
                # find the node in this group with the smallest free energy (largest probability)
                P = max([self.Peq[x] for x in nodes])
                for x in nodes:
                    assert x not in rate_estimates
                    rate_estimates[x] = ts_free_energy / P
#                    rate_estimates[x] = ts_free_energy / self.Peq[x]
                
                graph.remove_nodes_from(nodes)
            else:
                graph.add_edge(u,v)
                    

        self.rate_estimates = rate_estimates
                
                
            
        
        
    


class CommittorLinalg(object):
    """compute committor probabilites using sparse linear algebra"""
    def __init__(self, rates, A, B, debug=False, weights=None):
        self.rates = rates
        self.A = set(A)
        self.B = set(B)
        self.nodes = set()
        for u, v in self.rates.iterkeys():
            self.nodes.add(u)
            self.nodes.add(v)
        
        self.time_solve = 0.
        
    def make_matrix(self):
        intermediates = self.nodes - self.A - self.B
        
        node_list = list(intermediates)
        n = len(node_list)
        matrix = scipy.sparse.dok_matrix((n,n))
        right_side = np.zeros(n)
        node2i = dict([(node,i) for i, node in enumerate(node_list)])
        
        
        for uv, rate in self.rates.iteritems():
            u, v = uv
#            v, u = uv

            if u in intermediates:
                iu = node2i[u]
                matrix[iu,iu] -= rate

                if v in intermediates:
                    matrix[iu, node2i[v]] = rate
        
                if v in self.B:
                    right_side[iu] -= rate
        
        self.node_list = node_list
        self.node2i = node2i
        self.matrix =  matrix.tocsr()
        self.right_side = right_side
        
    def compute_committors(self):
        self.make_matrix()
        if self.right_side.size == 1:
            # some versions of scipy can't handle matrices of size 1
            committors = np.array([self.right_side[0] / self.matrix[0,0]])
        else:
            t0 = time.clock()
            committors = scipy.sparse.linalg.spsolve(self.matrix, self.right_side)
            self.time_solve += time.clock() - t0
        eps = 1e-10
        if np.any(committors < -eps) or np.any(committors > 1+eps):
            qmax = committors.max()
            qmin = committors.min()
            raise LinalgError("The committors are not all between 0 and 1.  max=%.18g, min=%.18g" % (qmax, qmin))
        self.committor_dict = dict(((node, c) for node, c in izip(self.node_list, committors)))
#        self.committors = committors
#        print "committors", committors
        return self.committor_dict
    

class MfptLinalgSparse(object):
    """compute mean first passage times using sparse linear algebra"""
    def __init__(self, rates, B, sum_out_rates=None, check_graph=True):
        self.rates = rates
        self.B = set(B)
        if check_graph:
            self.rates = reduce_rates(self.rates, B)
    
        if False:
            self._make_subgroups()
        
        self.nodes = set()
        for u, v in self.rates.iterkeys():
            self.nodes.add(u)
            self.nodes.add(v)
        
        if sum_out_rates is None:
            self.sum_out_rates = compute_sum_out_rates(self.rates)
        else:
            self.sum_out_rates = sum_out_rates
        
        self.mfpt_dict = dict()
        self.time_solve = 0.
    
    def _make_subgroups(self):
        graph = nx.Graph()
        graph.add_edges_from(filter(lambda uv: uv[0] not in self.B and uv[1] not in self.B,
                                    self.rates.iterkeys()))
        cc = nx.connected_components(graph)
        self.subgroups = [set(c) for c in cc]
        print len(self.subgroups), "subgroups"
        if len(self.subgroups) <= 10:
            print "subgroup sizes", [len(c) for c in self.subgroups]
        else:
            print "subgroup sizes", [len(c) for c in self.subgroups[:10]], "..."
        
    def make_matrix(self, intermediates):
        assert not self.B.intersection(intermediates)
        
        node_list = list(intermediates)
        n = len(node_list)
        matrix = scipy.sparse.dok_matrix((n,n))
        node2i = dict([(node,i) for i, node in enumerate(node_list)])
        
        for iu, u in enumerate(node_list):
            matrix[iu,iu] = -self.sum_out_rates[u]
        
        for uv, rate in self.rates.iteritems():
            u, v = uv
            if u in intermediates and v in intermediates: 
                ui = node2i[u]
                vi = node2i[v]
                assert ui != vi
                matrix[ui,vi] = rate
        
        self.node_list = node_list
        self.node2i = node2i
        self.matrix =  matrix.tocsr()
    
    def compute_mfpt(self, use_umfpack=True, cg=False, mfpt_estimate=None):
        if not hasattr(self, "matrix"):
            self.make_matrix(self.nodes - self.B)
        t0 = time.clock()
        if cg:
            x0 = np.array([mfpt_estimate[u] for u in self.node_list])
            times, info = scipy.sparse.linalg.cgs(self.matrix, -np.ones(self.matrix.shape[0]),
                                                  x0=x0)
            print "time to solve using conjugate gradient", time.clock() - t0
        else:
            times = scipy.sparse.linalg.spsolve(self.matrix, -np.ones(self.matrix.shape[0]),
                                                use_umfpack=use_umfpack)
        self.time_solve += time.clock() - t0
        self.mfpt_dict = dict(((node, time) for node, time in izip(self.node_list, times)))
        if np.any(times < 0):
            raise LinalgError("error the mean first passage times are not all greater than zero")
        return self.mfpt_dict
    
    def compute_mfpt_symmetric(self, Peq):
        """make the matrix symmetric by multiplying both sides of the equation by Peq
        
            sum_v Peq_u k_uv = -Peq_u
        """
        intermediates = self.nodes - self.B
        
        node_list = list(intermediates)
        n = len(node_list)
        matrix = scipy.sparse.dok_matrix((n,n))
        node2i = dict([(node,i) for i, node in enumerate(node_list)])
        
        right_side = -np.array([Peq[u] for u in node_list])
        
        for iu, u in enumerate(node_list):
            matrix[iu,iu] = -self.sum_out_rates[u] * Peq[u]
        
        for uv, rate in self.rates.iteritems():
            u, v = uv
            if u in intermediates and v in intermediates: 
                ui = node2i[u]
                vi = node2i[v]
                assert ui != vi
                matrix[ui,vi] = rate * Peq[u]
        
        node_list = node_list
        node2i = node2i
        matrix =  matrix.tocsc()
        
        t0 = time.clock()
        from scikits.sparse.cholmod import cholesky
        factor = cholesky(matrix)
        times = factor(right_side)
#        times = scikits.sparse.spsolve(matrix, right_side,
#                                            use_umfpack=True)
        print "time solving symmetric linalg", time.clock() - t0
        self.time_solve += time.clock() - t0
        self.mfpt_dict = dict(((node, time) for node, time in izip(node_list, times)))
        if np.any(times < 0):
            raise LinalgError("error the mean first passage times are not all greater than zero")
        return self.mfpt_dict

    def compute_mfpt_from_estimate(self, mfpt_estimates):
        intermediates = self.nodes - self.B
        
        node_list = list(intermediates)
        n = len(node_list)
        matrix = scipy.sparse.dok_matrix((n,n))
        node2i = dict([(node,i) for i, node in enumerate(node_list)])
        
        right_side = -np.ones(len(node_list))
        
        for iu, u in enumerate(node_list):
            matrix[iu,iu] = -self.sum_out_rates[u] * mfpt_estimates[u]
        
        for uv, rate in self.rates.iteritems():
            u, v = uv
            if u in intermediates and v in intermediates: 
                ui = node2i[u]
                vi = node2i[v]
                assert ui != vi
                matrix[ui,vi] = rate * mfpt_estimates[v]
        
        matrix_max = np.max(matrix.values())
        print "matrix max value", np.max(matrix.values())
        print "matrix min value", np.min(matrix.values())
        print "matrix min abs value", np.min([np.abs(v) for v in matrix.values()])
#        for ij, v in matrix.iteritems():
#            matrix[ij] = v / matrix_max
#        print "new matrix max value", np.max(matrix.values())
#        print "new matrix min value", np.min(matrix.values())

            
        
        matrix =  matrix.tocsc()
        
        t0 = time.clock()
        cg = True
        if cg:
            times, info = scipy.sparse.linalg.cgs(matrix, right_side)
            print "time to solve using conjugate gradient", time.clock() - t0
        else:
            times = scipy.sparse.linalg.spsolve(matrix, right_side,
                                                use_umfpack=True)
            print "time solving symmetric linalg", time.clock() - t0
        self.time_solve += time.clock() - t0
        self.mfpt_dict = dict(((u, time * mfpt_estimates[u]) for u, time in izip(node_list, times)))
        if np.any(times < 0):
            raise LinalgError("error the mean first passage times are not all greater than zero")
        return self.mfpt_dict

    def compute_mfpt_subgroups(self, use_umfpack=True):
        for group in self.subgroups:
            self.make_matrix(group)
            t0 = time.clock()
            times = scipy.sparse.linalg.spsolve(self.matrix, -np.ones(self.matrix.shape[0]),
                                                use_umfpack=use_umfpack)
            self.time_solve += time.clock() - t0
            for node, time in izip(self.node_list, times):
                self.mfpt_dict[node] = time


class TwoStateRates(object):
    """compute committors and several different rates between two groups"""
    def __init__(self, rate_constants, A, B, weights=None, check_rates=True):
        if check_rates:
            self.rate_constants = reduce_rates(rate_constants, B, A=A)
        else:
            self.rate_constants = rate_constants
        self.A = A
        self.B = B
        self.weights = weights
        if self.weights is None:
            self.weights = dict([(a, 1.) for a in self.A])
        
        self.sum_out_rates = compute_sum_out_rates(self.rate_constants)
        self.mfpt_computer = MfptLinalgSparse(self.rate_constants, self.B,
                                              sum_out_rates=self.sum_out_rates, 
                                              check_graph=False)
        
        

    def get_rate_AB(self):
        """return the rate from A to B
        
        the rate is the inverse mean first passage time averaged over the nodes in A
        """
        rate = sum((self.weights[a] / self.mfptimes[a] for a in self.A))
        norm = sum((self.weights[a] for a in self.A))
        
        return rate / norm
    
    def get_alternate_committors(self, Agroup, Bgroup):
        """for each node a in A, compute the probability that it ends up in 
        B before coming back to itself or reaching another node in A.
        """
        a_committors = dict([(a, 0.) for a in Agroup])
        for uv, rate in self.rate_constants.iteritems():
            u, v = uv
            if u in Agroup and v not in Agroup:
                if v in Bgroup:
                    vcomm = 1.
                else:
                    vcomm = self.committor_dict[v]
                a_committors[u] += rate * vcomm
        for a in Agroup:
            a_committors[a] /= self.sum_out_rates[a]
        return a_committors

    def get_rate_AB_SS(self):
        """
        return the steady state rate from A to B
        """
        # for each node a in A, compute the probability that it ends up in
        # B before coming back to itself or reaching another node in A.
        a_committors = self.get_alternate_committors(self.A, self.B)
        
        rate = sum((self.weights[a] * a_committors[a] * self.sum_out_rates[a] for a in self.A))
        norm = sum((self.weights[a] for a in self.A))
        
        return rate / norm

    def get_committor(self, x):
        """return the probability that a trajectory starting from x reaches B before A"""
        if x in self.A:
            return 0.
        elif x in self.B:
            return 1.
        else:
            return self.committor_dict[x]
    
    def compute_rates(self, use_umfpack=True, subgroups=False, symmetric=False, Peq=None,
                      cg=False, mfpt_estimates=None):
        """compute the mean first passage times."""
        if mfpt_estimates is not None and not cg:
            self.mfptimes = self.mfpt_computer.compute_mfpt_from_estimate(mfpt_estimates)
        elif symmetric:
            assert Peq is not None
            self.mfptimes = self.mfpt_computer.compute_mfpt_symmetric(Peq)
            return
        elif subgroups:
            self.mfptimes = self.mfpt_computer.compute_mfpt_subgroups(use_umfpack=use_umfpack)
        else:
            self.mfptimes = self.mfpt_computer.compute_mfpt(use_umfpack=use_umfpack, cg=cg,
                                                            mfpt_estimate=mfpt_estimates)
    
    def compute_committors(self):
        """compute the committors""" 
        self.committor_computer = CommittorLinalg(self.rate_constants, self.A, self.B)
        self.committor_dict = self.committor_computer.compute_committors()
        return self.committor_dict

