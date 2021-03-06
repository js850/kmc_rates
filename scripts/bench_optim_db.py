from pele.utils.optim_compatibility import OptimDBConverter
from pele.rates import RateCalculation
from pele.storage import Database
import networkx as nx
import numpy as np
import scipy.sparse.linalg
import time

from rates_linalg import MfptLinalgSparse, TwoStateRates, EstimateRates

#db = Database("db.cf.sqlite")
db = Database()
if db.number_of_minima() == 0:
    converter = OptimDBConverter(db, mindata="min.data", 
             tsdata="ts.data")
    converter.pointsmin_data = None
    converter.pointsts_data = None
    converter.ReadMinDataFast()
    converter.ReadTSdataFast()



m1 = db.getMinimum(1)
m2 = db.getMinimum(3)
print "m1, m2", m1._id, m2._id
A = [m1]
B = [m2]
print "energy of mA", m1.energy
print "energy of mB", m2.energy

# remove disconnected components
# graph = nx.Graph()
# for ts in db.transition_states():
#     m1 = ts.minimum1
#     m2 = ts.minimum2
#     graph.add_edge(m1, m2)
# cc = nx.node_connected_component(graph, A[0])
# nodes = set(cc)
# if len(nodes) != db.number_of_minima():
#     print "discarding", db.number_of_minima() - len(nodes), "minima not connected to A"
#     print len(nodes), "minima remaining"
# assert B[0] in nodes
# transition_states = filter(lambda ts: ts.minimum1 in nodes and ts.minimum2 in nodes, 
#                            db.transition_states())

T=1.
print "temperature", T
if True:
    tstates = list(db.transition_states())
    print len(tstates)
pele_rates = RateCalculation(db.transition_states(), A, B, T=T, use_fvib=True)
pele_rates._make_kmc_graph()
rates = pele_rates.rate_constants
rate_norm = np.exp(-pele_rates.max_log_rate)
Peq = pele_rates._get_equilibrium_occupation_probabilities()

if True:
    print "saving rate constants to"
    with open("test_rate_consts", "w") as fout:
        for (u,v), k in sorted(rates.iteritems(), key=lambda m: (m[0][0]._id, m[0][1]._id)):
            fout.write("%6d %6d %s\n" % (u._id,v._id,k/rate_norm))

if False:
    print "making all rates 1"
    for uv in rates.iterkeys():
        rates[uv] = 1.
        u, v = uv
#        assert u != v
        assert (v,u) in rates

if False:
    print "remove self transition states"
    nodes = set([u for u,v in rates.iterkeys() if u==v])
    for u in nodes:
        rates.pop((u,u))
    
        


print "max rate constant", max(rates.itervalues())
print "min rate constant", min(rates.itervalues())

#if True:
#    print "applying a minimum rate to the rate matrix. breaking the rate matrix"
#    rmin = 1e-3
#    for uv, r in rates.iteritems():
#        if r < rmin:
#            rates[uv] = rmin



lin = MfptLinalgSparse(rates, B)
if True:
    lin.make_matrix(lin.nodes - lin.B)
    t0 = time.clock()
    mfpt = lin.compute_mfpt(use_umfpack=True)
else:
    t0 = time.clock()
    lin.compute_mfpt_subgroups(use_umfpack=True)
    mfpt = lin.time_dict
t1 = time.clock()
times = lin.mfpt_dict
print "mean first passage times"
print times[A[0]] / np.exp(pele_rates.max_log_rate)
print "rates"
print 1./times[A[0]] * np.exp(pele_rates.max_log_rate)
print "sparse linalg finished in", t1-t0, "seconds"
print "max time", max(mfpt.itervalues())
print "min time", min(mfpt.itervalues())

if True:
    print "\nestimating rates"
    estimator = EstimateRates(lin.rates, Peq, B)
    print "rate estimate", estimator.rate_estimates[A[0]] / rate_norm
    estimates = dict()
    for u, t in lin.mfpt_dict.iteritems():
        kcalc = 1. / t / rate_norm
        kest = estimator.rate_estimates[u] / rate_norm
#        print u._id, kcalc, kest, kcalc / kest
        estimates[u] = (kcalc, kest)
    print "max rate", max([kcalc for kcalc, kest in estimates.values()])
    print "min rate", min([kcalc for kcalc, kest in estimates.values()])
    print "max ratio", max([kest / kcalc for kcalc, kest in estimates.values()])
    print "min ratio", min([kest / kcalc for kcalc, kest in estimates.values()])

if True:
    print "computing rates using symmetric method"
    lin.compute_mfpt_symmetric(Peq)
    print "rate symetric", 1./lin.mfpt_dict[A[0]] / rate_norm
            
if False:
    print "computing rates using conjugant gradient method"
    lin.compute_mfpt(cg=True)
    print "rate cg", 1./lin.mfpt_dict[A[0]] / rate_norm


if True:
    print "computing committors and steady state rate constants"
    tsr = TwoStateRates(rates, A, B)
    tsr.compute_rates()
    rAB = tsr.get_rate_AB() * np.exp(pele_rates.max_log_rate)
    print "rate AB", rAB
    tsr.compute_committors()
    rABss = tsr.get_rate_AB_SS() * np.exp(pele_rates.max_log_rate)
    print "steady state rate", rABss

if False:
    print "saving the graph structure"
    nodes = list(lin.nodes)
    node2i = dict([(node, i) for i, node in enumerate(nodes)])
    node2i = dict([(node, node._id) for node in nodes])
    with open("error.graph.id", "w") as fout:
        fout.write("# %d\n" % (node2i[iter(B).next()]))
        for uv in lin.rates.iterkeys():
            i = node2i[uv[0]]
            j = node2i[uv[1]]
            if i < j:
                fout.write("%d %d\n" % (i,j))


if False:
    print "\nrecomputing rates"
    rates2 = dict([((uv[0]._id, uv[1]._id),rate) for uv, rate in lin.rates.iteritems()])
    lin2 = MfptLinalgSparse(rates2, [B[0]._id])
    mfpt = lin2.compute_mfpt()
    print "new max time", mfpt.max()
    print "new min time", mfpt.min()
    import pickle
    pickle.dump(rates2, open("error.graph.id2", "w"))
    

if False:
    print "\ncomputing for reduced matrix"
    lin.make_matrix(lin.independent_sets[0] - lin.B)
    mfpt = lin.compute_mfpt()
    print 1./ mfpt[lin.node2i[A[0]]] * np.exp(pele_rates.max_log_rate)


if False:
    print "diagonalizing the matrix"
    w, v = scipy.sparse.linalg.eigs(lin.matrix)
    print "min max eigenvalues"
    print w.max(), w.min()
    import numpy as np
    np.savetxt("bench.eigs.txt", sorted(np.real(w)))

if False:
    print "computing the LU decomposition"
    lu = inverse = scipy.sparse.linalg.splu(lin.matrix)
    newtimes = lu.solve(-np.ones(lin.matrix.shape[0]))
    newtimes /= np.exp(pele_rates.max_log_rate)
    t = newtimes[lin.node2i[A[0]]] 
    print "LU times", t, 1./t
    print "max min LU times"
    print newtimes.max(), newtimes.min()

if False:
    print "smallest magnitude eigenvalue"
    u, v = scipy.sparse.linalg.eigs(lin.matrix, which="SM", k=1, tol=1e-6)
    print u
                             
#for t in mfpt:
#    print t, 1./t

#if True:
#    print pele_rates.compute_rates()