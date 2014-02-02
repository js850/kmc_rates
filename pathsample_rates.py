import argparse
import sys
import os
import time
import numpy as np
import networkx as nx

from pele.utils.optim_compatibility import OptimDBConverter
from pele.storage import Database
from pele.rates import RateCalculation

from rates_linalg import TwoStateRates, reduce_rates

def _check_AB(connected_components, A, AB="A"):
    Aclist = [A.intersection(c) for c in connected_components]
    Aclist = filter(lambda c:len(c)>0, Aclist)
    Aconn = set()
    for c in Aclist: Aconn.update(c)
    if Aconn != A:
        print "the following", AB, "nodes are not connected at all"
        print [m._id for m in A - Aconn]
    if len(Aclist) > 1:
        print "the following groups of", AB, "minima are connected within the group but not between groups"
        for c in Aclist:
            print len(c), "minima:", [m._id for m in c]

def analyze_graph_error(rates, A, B):
    A = set(A)
    B = set(B)
    
    if A.intersection(B):
        print "the following minima are in both A and B"
        print [m._id for m in A.intersection(B)]
    
    graph = nx.Graph()
    graph.add_edges_from(rates.iterkeys())
    
    # remove nodes not connected to B
    # TODO: this only works if B is fully connected
    cclist = nx.connected_components(graph)
    cclist = [set(c) for c in cclist]
    
    _check_AB(cclist, A, AB="A")
    _check_AB(cclist, B, AB="B")
    
     

def read_minA(fname, db):
    with open(fname) as fin:
        ids = []
        for i, line in enumerate(fin):
            if i == 0:
                nminima = int(line.split()[0])
            else:
                sline = line.split()
                ids += map(int, sline)
    
    assert nminima == len(ids)
    return [db.getMinimum(n) for n in ids]

def load_database(directory):
    print "reading from directory:", os.path.abspath(directory)
    db = Database()
    converter = OptimDBConverter(db, mindata=directory+"/min.data", 
             tsdata=directory+"/ts.data")
    converter.pointsmin_data = None
    converter.pointsts_data = None
    converter.ReadMinDataFast()
    converter.ReadTSdataFast()

    A = read_minA(directory+"/min.A", db)
    B = read_minA(directory+"/min.B", db)
    print len(A), "A minima read from min.A"
    print len(B), "B minima read from min.B"
    
    return db, A, B

def make_rates(database, A, B, T):
    print "computing rates from transition states"
    pele_rates = RateCalculation(database.transition_states(), A, B, T=T, use_fvib=True)
    pele_rates._make_kmc_graph()
    rates = pele_rates.rate_constants
    print "computing equilibrium occupation probabilities"
    Peq = pele_rates._get_equilibrium_occupation_probabilities()
    
    return rates, Peq, np.exp(-pele_rates.max_log_rate)


def main():
    tstart =  time.clock()
    parser = argparse.ArgumentParser(description="compute rates from a pathsample database")

    parser.add_argument("-d", type=str, default=".",
                        help="directory with the min.data, ts.data files")
    parser.add_argument("-T", type=float, default=1., help="Temperature")
    args = parser.parse_args()
    print "temperature", args.T
    
    # load the database
    db, A, B = load_database(args.d)
    # compute rate constants
    rate_constants, Peq, knorm = make_rates(db, A, B, args.T)
    print "checking and reducing the graph structure"
    try:
        rate_constants = reduce_rates(rate_constants, B, A=A)
    except Exception:
        analyze_graph_error(rate_constants, A, B)
        raise
        
    
    calculator = TwoStateRates(rate_constants, A, B, weights=Peq, check_rates=False)
    calculator.compute_rates(use_umfpack=True)
    kAB = calculator.get_rate_AB()
    
    print "k(B<-A)", kAB / knorm
    
    if True:
        fname = "out.mfpt"
        print "saving mean first passage times for all minima to reach B to file", fname
        mfpt = sorted([(m._id, t) for m, t in 
                       calculator.mfpt_computer.mfpt_dict.iteritems()])
        with open(fname, "w") as fout:
            for i, t in mfpt:
                fout.write("%d %.12g\n" % (i, t * knorm))

    print "computing committor probabilities"
    sys.stdout.flush()
    calculator.compute_committors()
    print "kSS(B<-A)", calculator.get_rate_AB_SS() / knorm
    
    if True:
        fname = "out.committors"
        print "saving committor probabilities for all minima to file", fname
        coms = sorted([(m._id, c) for m, c in 
                       calculator.committor_computer.committor_dict.iteritems()])
        with open(fname, "w") as fout:
            for i, c in coms:
                fout.write("%d %.12g\n" % (i, c))
    
    time_solve = calculator.mfpt_computer.time_solve + calculator.committor_computer.time_solve
    print "time spent solving linear equations", time_solve, "seconds"
    print "total time", time.clock() - tstart
     


if __name__ == "__main__":
    main()