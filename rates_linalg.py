import time
from itertools import izip
from collections import defaultdict

import numpy as np
import scipy.sparse
import scipy.sparse.linalg
import networkx as nx
from scipy.weave.catalog import intermediate_dir


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
    
    def compute_mfpt(self, use_umfpack=True, cg=False):
        if not hasattr(self, "matrix"):
            self.make_matrix(self.nodes - self.B)
        t0 = time.clock()
        if cg:
            times, info = scipy.sparse.linalg.cgs(self.matrix, -np.ones(self.matrix.shape[0]))
            print "time to solve using conjugate gradient", time.clock() - t0
        else:
            times = scipy.sparse.linalg.spsolve(self.matrix, -np.ones(self.matrix.shape[0]),
                                                use_umfpack=use_umfpack)
        self.time_solve += time.clock() - t0
        self.mfpt_dict = dict(((node, time) for node, time in izip(self.node_list, times)))
        if np.any(times < 0):
            raise RuntimeError("error the mean first passage times are not all greater than zero")
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
            raise RuntimeError("error the mean first passage times are not all greater than zero")
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
    
    def get_rate_AB_SS(self):
        """
        return the steady state rate from A to B
        """
        # for each node a in A, compute the probability that it ends up in
        # B before coming back to itself or reaching another node in A.
        a_committors = dict([(a, 0.) for a in self.A])
        for uv, rate in self.rate_constants.iteritems():
            u, v = uv
            if u in self.A and v not in self.A:
                if v in self.B:
                    vcomm = 1.
                else:
                    vcomm = self.committor_dict[v]
                a_committors[u] += rate * vcomm
        for a, c in a_committors.iteritems():
            a_committors[a] /= self.sum_out_rates[a]
        # the sum_out_rates cancels here, we can remove it
        
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
                      cg=False):
        """compute the mean first passage times."""
        if symmetric:
            assert Peq is not None
            self.mfptimes = self.mfpt_computer.compute_mfpt_symmetric(Peq)
            return
        if subgroups:
            self.mfptimes = self.mfpt_computer.compute_mfpt_subgroups(use_umfpack=use_umfpack)
        else:
            self.mfptimes = self.mfpt_computer.compute_mfpt(use_umfpack=use_umfpack, cg=cg)
    
    def compute_committors(self):
        """compute the committors""" 
        self.committor_computer = CommittorLinalg(self.rate_constants, self.A, self.B)
        self.committor_dict = self.committor_computer.compute_committors()
        
