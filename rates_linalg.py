
import numpy as np
import scipy.sparse
import scipy.sparse.linalg
from numpy import bench


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
    
#    def compute_rate(self):
#        rate = 0.
#        for u in self.nodes:
#            for v in self.B:
#                if u != v:
#                    if u in self.A:
#                        com = 0.
#                    else:
#                        com = self.committors[self.node2i[u]]
#                    rate += self.rates[(u,v)] * (1. - com)
#        print "rate AB", rate
#        return rate
#
#    def compute_eq_prob_pysal(self):
#        from kmc_rates import kmcgraph_from_rates
#        from pysal.spatial_dynamics.ergodic import steady_state, fmpt
#        node_list = list(self.nodes)
#        graph = kmcgraph_from_rates(self.rates)
#        Tmat = np.zeros([len(self.nodes)]*2)
#        self.node2i_full = dict([(node, i) for i, node in enumerate(node_list)])
#        for u, v, data in graph.out_edges_iter(data=True):
#            Tmat[self.node2i_full[u], self.node2i_full[v]] = data["P"]
#        
#        pi = steady_state(Tmat)
#        self.eq_prob = dict()
#        for i, node in enumerate(node_list):
#            self.eq_prob[node] = np.real(pi[i])
#        print "Tmat",
#        print Tmat
#        print "pysal fmpt", fmpt(Tmat)
#        print "tau[", node_list[2],"]", graph.node[node_list[2]]["tau"]
#        print "pysal eq prob", np.real(pi)
#        
#
#    def compute_eq_prob(self):
#        rin = dict([(u,0.) for u in self.nodes])
#        rout = dict([(u,0.) for u in self.nodes])
#        for uv, rate in self.rates.iteritems():
#            u, v = uv
#            rout[u] += rate
#            rin[v] += rate
#        
#        self.eq_prob = dict()
#        for u in self.nodes:
#            self.eq_prob[u] = rin[u] / rout[u]
#        
#        norm = sum(self.eq_prob.itervalues())
#        for u in self.nodes:
#            self.eq_prob[u] /= norm
#        
#        #testing
#        u, v = 2,3
#        print "test detailed balance"
#        print self.eq_prob[u] * self.rates[(u,v)], self.eq_prob[v] * self.rates[(v,u)]
#        
#            
#    def compute_rate2(self):
#        self.compute_eq_prob_pysal()
#        weights = dict([(node, 1.) for node in self.nodes])
#        weights = self.eq_prob
#        
#        rate = 0.
#        norm = 0.
#        for u in self.nodes:
#            if u in self.B: continue
#            for v in self.B:
#                if u != v:
#                    if (u,v) in self.rates:
#                        if u in self.A:
#                            com = 0.
#                        else:
#                            com = self.committors[self.node2i[u]]
#                        print "r", u,v, weights[u], self.rates[(u,v)], com
#                        rate += weights[u] * self.rates[(u,v)] * (1. - com)
#                        if u in self.A:
#                            norm += weights[u]
##        rate /= norm
#        print "rate AB", rate, norm, rate / norm
#        return rate
#
#    def compute_rate3(self):
#        self.A = set()
#        self.make_matrix()
#        times = np.linalg.solve(self.matrix, -np.ones(self.matrix.shape[0]))
#        return 1./times
        
class MfptLinalg(object):
    def __init__(self, rates, B):
        self.rates = rates
        self.B = set(B)
        self.nodes = set()
        for u, v in self.rates.iterkeys():
            self.nodes.add(u)
            self.nodes.add(v)
        
    def make_matrix(self):
        intermediates = self.nodes - self.B
        
        nodes = list(intermediates)
        n = len(nodes)
        matrix = np.zeros([n,n])
        node2i = dict([(node,i) for i, node in enumerate(nodes)])
        
        
        for uv, rate in self.rates.iteritems():
            u, v = uv
#            v, u = uv

            if u in intermediates:
                iu = node2i[u]
                matrix[iu,iu] -= rate

                if v in intermediates:
                    matrix[iu, node2i[v]] = rate
        
        
        self.node_list = nodes
        self.node2i = node2i
        self.matrix =  matrix
#        print "matrix", self.matrix
#        print self.matrix
    
    def compute_mfpt(self):
        self.make_matrix()
        times = np.linalg.solve(self.matrix, -np.ones(self.matrix.shape[0]))
        return times

class MfptLinalgSparse(object):
    def __init__(self, rates, B):
        self.rates = rates
        self.B = set(B)
        self.nodes = set()
        for u, v in self.rates.iterkeys():
            self.nodes.add(u)
            self.nodes.add(v)
        
    def make_matrix(self):
        intermediates = self.nodes - self.B
        
        nodes = list(intermediates)
        n = len(nodes)
        matrix = scipy.sparse.dok_matrix((n,n))
        node2i = dict([(node,i) for i, node in enumerate(nodes)])
        
        
        for uv, rate in self.rates.iteritems():
            u, v = uv
#            v, u = uv

            if u in intermediates:
                iu = node2i[u]
                matrix[iu,iu] -= rate

                if v in intermediates:
                    matrix[iu, node2i[v]] = rate
        
        
        self.node_list = nodes
        self.node2i = node2i
#        print "matrix", matrix
        self.matrix =  matrix.tocsr()
    
    def compute_mfpt(self):
        if not hasattr(self, "matrix"):
            self.make_matrix()
        times = scipy.sparse.linalg.spsolve(self.matrix, -np.ones(self.matrix.shape[0]))
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
    nlist = 1.4**np.array(np.arange(15,30))
    nlist = [int(n) for n in nlist]
    print nlist
    tsplist = []
    plt.figure()
    plt.show()
    for n in nlist:
#        maker = _MakeRandomGraph(n, n, node_set=set(range(n)))
#        rates = maker.make_rates()    
        rates = make_sparse_network(n, 2*n)
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
        

        
        

if __name__ == "__main__":
#    test()
    benchmarks()
        