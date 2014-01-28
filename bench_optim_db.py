from pele.utils.optim_compatibility import OptimDBConverter
from pele.rates import RateCalculation
from pele.storage import Database
import networkx as nx
import numpy as np
import scipy.sparse.linalg
import time

from rates_linalg import MfptLinalgSparse

db = Database("db.40000.sqlite")
#if db.number_of_minima() == 0:
#    converter = OptimDBConverter(db, mindata="min.data.40000", 
#             tsdata="ts.data.40000")
#    converter.pointsmin_data = None
#    converter.pointsts_data = None
#    converter.ReadMindata()
#    converter.ReadTSdata()



m1 = db.getMinimum(18)
m2 = db.getMinimum(1)
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


pele_rates = RateCalculation(db.transition_states(), A, B, T=.592, use_fvib=True)
pele_rates._make_kmc_graph()

rates = pele_rates.rate_constants
print "max rate constant", max(rates.itervalues())
print "min rate constant", min(rates.itervalues())

#if True:
#    print "applying a minimum rate to the rate matrix. breaking the rate matrix"
#    rmin = 1e-3
#    for uv, r in rates.iteritems():
#        if r < rmin:
#            rates[uv] = rmin



lin = MfptLinalgSparse(rates, B)
lin.make_matrix(lin.nodes - lin.B)
t0 = time.clock()
mfpt = lin.compute_mfpt()
t1 = time.clock()
times = lin.time_dict
print "mean first passage times"
print times[A[0]] / np.exp(pele_rates.max_log_rate)
print "rates"
print 1./times[A[0]] * np.exp(pele_rates.max_log_rate)
print "sparse linalg finished in", t1-t0, "seconds"

print "max time", mfpt.max()
print "min time", mfpt.min()

if True:
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

if True:
    print "computing the LU decomposition"
    lu = inverse = scipy.sparse.linalg.splu(lin.matrix)
    newtimes = lu.solve(-np.ones(lin.matrix.shape[0]))
    newtimes /= np.exp(pele_rates.max_log_rate)
    t = newtimes[lin.node2i[A[0]]] 
    print "LU times", t, 1./t
    print "max min LU times"
    print newtimes.max(), newtimes.min()

if True:
    print "smallest magnitude eigenvalue"
    u, v = scipy.sparse.linalg.eigs(lin.matrix, which="SM", k=1, tol=1e-6)
    print u
                             
#for t in mfpt:
#    print t, 1./t

#if True:
#    print pele_rates.compute_rates()