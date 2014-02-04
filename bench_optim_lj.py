from pele.utils.optim_compatibility import OptimDBConverter
from pele.rates import RateCalculation
from pele.storage import Database, Minimum
import networkx as nx
import numpy as np
import scipy.sparse.linalg
import time
from rates_linalg import MfptLinalgSparse, TwoStateRates
import sys
tstart = time.clock()

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

print "computing rates from minima data"
T = 2.5
print "temperature", T
pele_rates = RateCalculation(db.transition_states(), A, B, T=T, use_fvib=True)
pele_rates._make_kmc_graph()
rates = pele_rates.rate_constants
weights = pele_rates._get_equilibrium_occupation_probabilities()

print "max rate constant", max(rates.itervalues())
print "min rate constant", min(rates.itervalues())



tsr = TwoStateRates(rates, A, B, weights=weights)
print "after removing unconnected we have", len(tsr.mfpt_computer.nodes), "nodes and", len(tsr.mfpt_computer.rates), "rates"
if False:
    tsr.compute_rates(symmetric=True, Peq=weights)
else:
    tsr.compute_rates(use_umfpack=False, cg=True)
kAB = tsr.get_rate_AB()
mfpt = tsr.mfpt_computer.mfpt_dict
print "rate AB", kAB * np.exp(pele_rates.max_log_rate)
print "sparse linalg finished in", tsr.mfpt_computer.time_solve, "seconds"
print "max mfpt time", max(mfpt.itervalues())
print "min mfpt time", min(mfpt.itervalues())
print ""
sys.stdout.flush()



if True:
    print "computing committors and steady state rate constants"
    tsr.compute_committors()
    rABss = tsr.get_rate_AB_SS() * np.exp(pele_rates.max_log_rate)
    print "steady state rate", rABss

print "total duration", time.clock() - tstart 