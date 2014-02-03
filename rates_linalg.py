import numpy as np
import scipy.sparse
import scipy.sparse.linalg
from itertools import izip
from collections import defaultdict
import networkx as nx
import time


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
    
    def compute_mfpt(self, use_umfpack=True):
        if not hasattr(self, "matrix"):
            self.make_matrix(self.nodes - self.B)
        t0 = time.clock()
        times = scipy.sparse.linalg.spsolve(self.matrix, -np.ones(self.matrix.shape[0]),
                                            use_umfpack=use_umfpack)
        self.time_solve += time.clock() - t0
        self.mfpt_dict = dict(((node, time) for node, time in izip(self.node_list, times)))
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
    
    def compute_rates(self, use_umfpack=True, subgroups=False):
        """compute the mean first passage times."""
        if subgroups:
            self.mfptimes = self.mfpt_computer.compute_mfpt_subgroups(use_umfpack=use_umfpack)
        else:
            self.mfptimes = self.mfpt_computer.compute_mfpt(use_umfpack=use_umfpack)
    
    def compute_committors(self):
        """compute the committors""" 
        self.committor_computer = CommittorLinalg(self.rate_constants, self.A, self.B)
        self.committor_dict = self.committor_computer.compute_committors()
        
             

def test():
#    from test_graph_transformation import _three_state_graph, _MakeRandomGraph
    from kmc_rates import GraphReduction, kmcgraph_from_rates
    from kmc import KineticMonteCarlo
    nnodes = 4
    rates = dict()
    for i in xrange(0,nnodes):
        for j in xrange(0,nnodes):
            if i != j:
                rates[(i,j)] = np.random.rand()        
                rates[(i,j)] = float(i+1)
#                rates[(i,j)] = 1.
    rates.pop((0,3))
#    rates.pop((3,0))
    print "rates", rates
    A=[0]
    B=[3]
    graph = kmcgraph_from_rates(rates)
#    import networkx as nx
#    import matplotlib.pyplot as plt
#    nx.draw(graph)
#    plt.draw()
#    plt.show()
#    maker = _MakeRandomGraph(nnodes=8, nedges=60, node_set=A+B)
#    graph = maker.run()
    
    c = CommittorLinalg(rates.copy(), A, B)
    cp = c.compute_committors()
    print "linalg committors"
    print cp
    
    
    from kmc_rates import GraphReduction
    red = GraphReduction(graph.copy(), A, B)
    cprob = red.compute_committor_probabilities(set(graph.nodes()) - set(A) - set(B))
    print "NGT committors"
    print cprob
    
#    print cp[0]
#    print cprob[B[0]+1]
    
    mfpt_comp = MfptLinalg(rates, B)
    mfptimes = mfpt_comp.compute_mfpt()
    print "linalg rates", 1./ mfptimes
    
    mfpt_comps = MfptLinalgSparse(rates, B)
    mfptimess = mfpt_comps.compute_mfpt()
    print "sparse linalg rates", 1./ mfptimess
    
    red = GraphReduction(graph.copy(), A, B)
    red.compute_rates()
    print "NGT rate AB", A, "->", B, red.get_rate_AB()
    print "NGT rate BA", red.get_rate_BA()
    print "NGT MFPT AB", 1./red.get_rate_AB()
    
#    print "NGT SS rate AB", red.get_SS_rate_AB()
    
    kmc = KineticMonteCarlo(graph.copy())
    tkmc = kmc.mean_first_passage_time(A[0], B, niter=10000)
    print "kmc mfp time", tkmc
    print "kmc rate", 1./tkmc

def make_sparse_network(nnodes, nts):
    rates = dict()
    for u in xrange(1,nnodes):
        v = np.random.randint(u)
        rates[(u,v)] = np.random.rand()
        rates[(v,u)] = np.random.rand()
    
    nodes = range(nnodes)
    while len(rates) < 2*nts:
        np.random.shuffle(nodes)
        u,v = nodes[:2]
        rates[(u,v)] = np.random.rand()
        rates[(v,u)] = np.random.rand()
    
    return rates
        

def benchmarks():
    from test_kmc import _MakeRandomGraph
    from kmc_rates import kmcgraph_from_rates, GraphReduction
    from matplotlib import pyplot as plt
    plt.ion()
    import sys
    import time
    nlist = 1.4**np.array(np.arange(15,40))
    nlist = [int(n) for n in nlist]
    print nlist
    tsplist = []
    plt.figure()
    plt.show()
    for n in nlist:
#        maker = _MakeRandomGraph(n, n, node_set=set(range(n)))
#        rates = maker.make_rates()    
        rates = make_sparse_network(n, n*1.1)
#        graph = kmcgraph_from_rates(rates)
        
        
        A = [0]
        B = [1]
        
        t0 = time.clock()
#        red = MfptLinalg(rates, B)
#        mfpt = red.compute_mfpt()
#        print 1./mfpt[0]
        
        red = MfptLinalgSparse(rates, B)
        red.make_matrix()
        t1 = time.clock()
        mfpt = red.compute_mfpt()
        t2 = time.clock()
        print 1./mfpt[0]
        tsplist.append(t2-t1)
        
#        red = GraphReduction(graph, A, B)
#        red.compute_rates()
#        print red.get_rate_AB()
#        t3 = time.clock()
        
        
        print n, len(rates) / 2, ": times", t1-t0, t2-t1 #, t3-t2
        sys.stdout.flush()
        
        if len(tsplist) > 3:
            plt.clf()
            plt.loglog(nlist[:len(tsplist)], tsplist, '-.')
            plt.draw()
#            plt.show()

        if True:
            out = np.zeros([len(tsplist),2])
            out[:,0] = nlist[:len(tsplist)]
            out[:,1] = tsplist
            np.savetxt("bench.dat", out)
    
    raw_input("press key")
#    vals = tsplist.items()
#    vals.sort(key=lambda v:v[0])
#    nlist = [n for n,t in vals]
#    times = [t for n,t in vals]
#    print nlist
#    print times
        

def benchmark_ngt():
    from test_kmc import _MakeRandomGraph
    from kmc_rates import kmcgraph_from_rates, GraphReduction
    from matplotlib import pyplot as plt
    plt.ion()
    import sys
    import time
    nlist = 1.4**np.array(np.arange(8,40))
    nlist = [int(n) for n in nlist]
    print nlist
    tsplist = []
    plt.figure()
    plt.show()
    for n in nlist:
#        maker = _MakeRandomGraph(n, n, node_set=set(range(n)))
#        rates = maker.make_rates()    
        rates = make_sparse_network(n, n*1.1)
        graph = kmcgraph_from_rates(rates)
        
        
        A = [0]
        B = [1]
        
        red = GraphReduction(graph, A, B)
        t0 = time.clock()
        red.compute_rates()
        print red.get_rate_AB()
        t1 = time.clock()
        tsplist.append(t1-t0)
        
        
        print n, len(rates) / 2, ": times", t1-t0
        sys.stdout.flush()
        
        if len(tsplist) > 5:
            if max(tsplist) > 0.1:
                print tsplist
                plt.clf()
                plt.loglog(nlist[:len(tsplist)], tsplist, '-.')
                plt.draw()

        if True:
            out = np.zeros([len(tsplist),2])
            out[:,0] = nlist[:len(tsplist)]
            out[:,1] = tsplist
            np.savetxt("bench.ngt.dat", out)
    
    raw_input("press key")
#    vals = tsplist.items()
#    vals.sort(key=lambda v:v[0])
#    nlist = [n for n,t in vals]
#    times = [t for n,t in vals]
#    print nlist
#    print times
        

        

if __name__ == "__main__":
#    test()
#    benchmarks()
    benchmark_ngt()
        