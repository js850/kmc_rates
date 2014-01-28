
import numpy as np
import scipy.sparse
import scipy.sparse.linalg
from numpy import bench
from itertools import izip


class CommittorLinalg(object):
    def __init__(self, rates, A, B, debug=False, weights=None):
        self.rates = rates
        self.A = set(A)
        self.B = set(B)
        self.nodes = set()
        for u, v in self.rates.iterkeys():
            self.nodes.add(u)
            self.nodes.add(v)
        self.intermediates = self.nodes - self.A - self.B
        
    def make_matrix(self):
        intermediates = self.nodes - self.A - self.B
        
        nodes = list(intermediates)
        n = len(nodes)
        matrix = np.zeros([n,n])
        right_side = np.zeros(n)
        node2i = dict([(node,i) for i, node in enumerate(nodes)])
        
        
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
        
        self.node_list = nodes
        self.node2i = node2i
        self.matrix =  matrix
        print "matrix", self.matrix
        self.right_side = right_side
        print self.matrix
        
    def compute_committors(self):
        self.make_matrix()
        committors = np.linalg.solve(self.matrix, self.right_side)
        self.committors = committors
#        print "committors", committors
        return committors
    

class MfptLinalgSparse(object):
    def __init__(self, rates, B):
        self.rates = rates
        self.B = set(B)
        self.initialize()
        
    def initialize(self):
        import networkx as nx
        graph = nx.Graph()
        graph.add_edges_from(self.rates.iterkeys())
        
        # remove nodes not connected to B
        # TODO: this only works if B is fully connected
        connected_nodes = nx.node_connected_component(graph, iter(self.B).next())
        connected_nodes = set(connected_nodes)
        all_nodes = set(graph.nodes())
        self.nodes = set(connected_nodes)
        if len(connected_nodes) != len(all_nodes):
            print "removing", len(all_nodes) - len(connected_nodes), "nodes that are not connected to B"
        
            self.rates = dict((uv, rate) for uv, rate in self.rates.iteritems()
                              if uv[0] in connected_nodes
                              )
        
#        # now remove the B nodes from the graph and see if it is split into multiple parts
#        graph = nx.Graph()
#        graph.add_edges_from(filter(lambda uv: uv[0] not in self.B and uv[1] not in self.B,
#                                    self.rates.iterkeys()))
#        graph.add_nodes_from(self.nodes)
#        cc = nx.connected_components(graph)
#        self.independent_sets = [set(c) for c in cc]
#        print len(self.independent_sets), "independent sets"
#        print [len(c) for c in self.independent_sets]
#        
#        # compute the sums of the rates out of each node.
#        # These will be the negative of the diagonal of the rate matrix        
#        self.sum_out_rates = dict()
#        for uv, rate in self.rates.iteritems():
#            u = uv[0]
#            try:
#                self.sum_out_rates[u] += rate
#            except KeyError:
#                self.sum_out_rates[u] = rate
        
        
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
#        print "matrix", matrix
        self.matrix =  matrix.tocsr()
    
    def _compute_mfpt(self, nodes):
        self.make_matrix(nodes)
    
    def compute_mfpt(self):
        if not hasattr(self, "matrix"):
            self.make_matrix(self.nodes - self.B)
        times = scipy.sparse.linalg.spsolve(self.matrix, -np.ones(self.matrix.shape[0]))
        self.time_dict = dict(((node, time) for node, time in izip(self.node_list, times)))
        if np.any(times < 0):
            print "error the mean first passage times are not all greater than zero"
        return times

              

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
        