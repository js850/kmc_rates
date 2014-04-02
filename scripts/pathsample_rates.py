"""Compute rates from a pathsample database.  The rates are computed
by solving a set of linear equations using a sparse solver.

Requirements
------------
Python2.7 and the packages numpy, scipy, and networkx.  
"""
import argparse
import sys
import os
import time
from itertools import izip
import numpy as np
import networkx as nx

from kmc_rates.rates_linalg import TwoStateRates, reduce_rates
from kmc_rates.utils import read_minA, make_rates, analyze_graph_error



description="""Compute rates from a pathsample database.  The transition
states and minima data are read from min.data and ts.data.  The product
and reactant states are read from min.A and min.B.  The rates are computed
by solving a set of linear equations using a sparse solver.
"""

def main():
    """the main loop"""
    tstart =  time.clock()
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("-d", type=str, default=".",
                        help="directory with the min.data, ts.data files")
    parser.add_argument("-T", type=float, default=1., help="Temperature")
    args = parser.parse_args()
    print "temperature", args.T
    directory = args.d
    print "reading from directory:", os.path.abspath(directory)
    
    A = read_minA(directory+"/min.A")
    B = read_minA(directory+"/min.B")
    rate_constants, Peq, knorm = make_rates(args.d, args.T)
            

    if True:
        fname = "out.rate_consts"
        print "saving rate constants to", fname
        with open(fname, "w") as fout:
            for (u,v), k in sorted(rate_constants.iteritems()):
                fout.write("%6d %6d %s\n" % (u,v,k/knorm))

    print "checking and reducing the graph structure"
    try:
        rate_constants = reduce_rates(rate_constants, B, A=A)
    except Exception:
        analyze_graph_error(rate_constants, A, B)
        raise
    
    print "computing mean first passage times"
    calculator = TwoStateRates(rate_constants, A, B, weights=Peq, check_rates=False)
    calculator.compute_rates(use_umfpack=True)
    print "computing rates"
    kAB = calculator.get_rate_AB()
    print "k(B<-A)", kAB / knorm
    
    if True:
        fname = "out.mfpt"
        print "saving mean first passage times for all minima to reach B to file", fname
        mfpt = sorted([(m, t) for m, t in 
                       calculator.mfpt_computer.mfpt_dict.iteritems()])
        with open(fname, "w") as fout:
            for i, t in mfpt:
                fout.write("%d %.12g\n" % (i, t * knorm))

    print "computing committor probabilities"
    sys.stdout.flush()
    calculator.compute_committors()
    print "computing steady state rates"
    kSS = calculator.get_rate_AB_SS() / knorm
    print "kSS(B<-A)", kSS
    
    if True:
        fname = "out.committors"
        print "saving committor probabilities for all minima to file", fname
        coms = sorted([(m, c) for m, c in 
                       calculator.committor_computer.committor_dict.iteritems()])
        with open(fname, "w") as fout:
            for i, c in coms:
                fout.write("%d %.12g\n" % (i, c))
    
    time_solve = calculator.mfpt_computer.time_solve + calculator.committor_computer.time_solve
    print "time spent solving linear equations", time_solve, "seconds"
    print "total time", time.clock() - tstart
     


if __name__ == "__main__":
    main()
