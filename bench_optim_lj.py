from pele.utils.optim_compatibility import OptimDBConverter
from pele.rates import RateCalculation
from pele.storage import Database, Minimum
import networkx as nx
import numpy as np
import scipy.sparse.linalg
import time
from rates_linalg import TwoStateRates, EstimateRates
from pele.utils.disconnectivity_graph import Tree
import sys
tstart = time.clock()

class TreeUnionFind(object):
    def __init__(self):
        self.trees = dict()

    def _root_tree(self, tree):
        if tree.parent is None:
            return tree
        else:
            return self._root_tree(tree.parent)

    def root_tree(self, u):
        try:
            tree = self.trees[u]
        except KeyError:
            tree = Tree()
            tree.data["minimum"] = u
            self.trees[u] = tree
            return tree
        return self._root_tree(tree)
    
    def __getitem__(self, u):
        return self.root_node(u)
    
    def root_node(self, u):
        return self.root_tree(u).data["minimum"]
    
    def union(self, u, v):
        uroot = self.root_tree(u)
        vroot = self.root_tree(v)
        if uroot != vroot:
            if v.energy > u.energy:
                vroot.add_branch(uroot)
            else:
                uroot.add_branch(vroot)
            return True
        else:
            return False

def find_nearest_neighbors(db, a, nmax=10, rates=None, Peq=None):
    if rates is None:
        edges = [ (ts.minimum1, ts.minimum2) for ts in 
                 db.transition_states(order_energy=True)]
    else:
        edges = [((u,v), k*Peq[u]) for (u,v), k in 
                 rates.iteritems() if u<v]
        edges.sort(key=lambda uvk: -uvk[1])
        assert edges[0][1] > edges[1][1]
        edges = [uv for uv, k in edges]
        
    subtrees = TreeUnionFind()
    subtrees[a]
    for u,v in edges:
        uroot = subtrees[u]
        vroot = subtrees[v]
        if uroot != vroot:
            subtrees.union(u,v)
            aroot = subtrees[a]
            uroot = subtrees[u]
#            graph.add_edge(u,v)
            if uroot == aroot:
                aroot_tree = subtrees.trees[aroot]
                if aroot_tree.number_of_subtrees() > nmax:
                    print "found a tree of sufficient size"
                    return [t.data["minimum"] for t in aroot_tree.get_all_trees()]
    aroot_tree = subtrees.root_tree(u)
    print "finished without finding a large tree"
    return [t.data["minimum"] for t in aroot_tree.get_all_trees()]

def find_nearest_neighbors2(db, a, nmax=10, rates=None, Peq=None):
    if rates is None:
        edges = [ (ts.minimum1, ts.minimum2) for ts in 
                 db.transition_states(order_energy=True)]
    else:
        edges = [((u,v), k*Peq[u]) for (u,v), k in 
                 rates.iteritems() if u<v]
        edges.sort(key=lambda uvk: -uvk[1])
        assert edges[0][1] > edges[1][1]
        edges = [uv for uv, k in edges]

    subtrees = nx.utils.UnionFind()
    subtrees[a]
    graph = nx.Graph()
    graph.add_node(a)
    for u,v in edges:
        uroot = subtrees[u]
        vroot = subtrees[v]
        if uroot != vroot:
            subtrees.union(u,v)
            graph.add_edge(u,v)
#            graph.add_edge(u,v)
            if subtrees[u] == subtrees[a]:
                cc = nx.node_connected_component(graph, a)
                if len(cc) >= nmax:
                    print "found a tree of sufficient size"
                    return [m for m in cc]
    print "finishing without having found a tree of sufficient size"
    cc = nx.node_connected_component(graph, a)
    return [m for m in cc]
    

def read_minA(fname, db):
    """load data from min.A or min.B"""
    with open(fname) as fin:
        ids = []
        for i, line in enumerate(fin):
            if i == 0:
                nminima = int(line.split()[0])
            else:
                sline = line.split()
                ids += map(int, sline)
    
    assert nminima == len(ids)
    print len(ids), "minima read from file:", fname
    return [db.getMinimum(mid) for mid in ids]


db = Database()
direc = "lj38/20000.minima"
if db.number_of_minima() == 0:
    converter = OptimDBConverter(db, mindata=direc+"/min.data", 
             tsdata=direc+"/ts.data")
    converter.pointsmin_data = None
    converter.pointsts_data = None
    converter.ReadMinDataFast()
    converter.ReadTSdataFast()




A = read_minA(direc+"/min.A", db)
B = read_minA(direc+"/min.B", db)
print len(A), "A minima"
print len(B), "B minima"

        
        

print A[0].energy
print B[0].energy

#rate_constants = reduce_rates(rate_constants, B, A=A)


print "computing rates from minima data"
T = .2
print "temperature", T
pele_rates = RateCalculation(db.transition_states(), A, B, T=T, use_fvib=True)
pele_rates._make_kmc_graph()
rates = pele_rates.rate_constants
weights = pele_rates._get_equilibrium_occupation_probabilities()
rate_norm = np.exp(-pele_rates.max_log_rate)
#Peq = pele_rates._get_equilibrium_occupation_probabilities()


print "max rate constant", max(rates.itervalues())
print "min rate constant", min(rates.itervalues())

if True:
    amin = sorted(A, key=lambda m:m.energy)[0]
    bmin = sorted(B, key=lambda m:m.energy)[0]
#    graph = nx.Graph()
#    Emin = db.minima()[0].energy
#    for ts in db.transition_states():
#        graph.add_edge(ts.minimum1, ts.minimum2, weight=ts.energy - Emin)
#    newA = set([amin])
#    dist, path = nx.single_source_dijkstra(graph, amin, cutoff=10)
#    dist = sorted([(d,u) for u, d in dist.iteritems()])
#    for d, v in dist:
#        if len(newA) >= 10:
#            break
#        newA.add(v)
#    newB = set([bmin])
#    dist, path = nx.single_source_dijkstra(graph, bmin, cutoff=10)
#    dist = sorted([(d,u) for u, d in dist.iteritems()])
#    for d, v in dist:
#        if len(newB) >= 10:
#            break
#        newB.add(v)
        
#    newA = find_nearest_neighbors2(db, amin, nmax=30)
#    newB = find_nearest_neighbors2(db, bmin, nmax=100)
    newA = find_nearest_neighbors2(db, amin, nmax=30, rates=rates, Peq=weights)
    newB = find_nearest_neighbors2(db, bmin, nmax=100, rates=rates, Peq=weights)
    print "new A", len(newA)
    print [m._id for m in newA]
    print "new B", len(newB)
    print [m._id for m in newB]
    exit()


tsr = TwoStateRates(rates, A, B, weights=weights)
print "after removing unconnected we have", len(tsr.mfpt_computer.nodes), "nodes and", len(tsr.mfpt_computer.rates), "rates"
if False:
    tsr.compute_rates(symmetric=True, Peq=weights)
else:
    tsr.compute_rates(use_umfpack=False, cg=False)
kAB = tsr.get_rate_AB()
mfpt = tsr.mfpt_computer.mfpt_dict
print "rate AB", kAB * np.exp(pele_rates.max_log_rate)
print "sparse linalg finished in", tsr.mfpt_computer.time_solve, "seconds"
print "max mfpt time", max(mfpt.itervalues())
print "min mfpt time", min(mfpt.itervalues())
print ""
sys.stdout.flush()

if True:
    print "\nestimating rates"
    lin = tsr
    estimator = EstimateRates(lin.rate_constants, weights, B)
    print "rate estimate", estimator.rate_estimates[A[0]] / rate_norm
    estimates = dict()
    for u, t in lin.mfpt_computer.mfpt_dict.iteritems():
        kcalc = 1. / t / rate_norm
        kest = estimator.rate_estimates[u] / rate_norm
        print u._id, kest, kcalc, kest / kcalc
        estimates[u] = (kcalc, kest)
    print "max rate", max([kcalc for kcalc, kest in estimates.values()])
    print "min rate", min([kcalc for kcalc, kest in estimates.values()])
    print "max overestimate ratio", max([kest / kcalc for kcalc, kest in estimates.values()])
    print "min underestimate ratio", min([kest / kcalc for kcalc, kest in estimates.values()])



if True:
    print "computing committors and steady state rate constants"
    tsr.compute_committors()
    rABss = tsr.get_rate_AB_SS() * np.exp(pele_rates.max_log_rate)
    print "steady state rate", rABss

print "total duration", time.clock() - tstart 