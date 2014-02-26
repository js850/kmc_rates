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

from pathsample_rates import make_rates, read_minA, reduce_rates, analyze_graph_error
from ngt_wrapper import NGT




description="""Compute rates from a pathsample database.  The transition
states and minima data are read from min.data and ts.data.  The product
and reactant states are read from min.A and min.B.  The rates are computed
by solving a set of linear equations using a sparse solver.
"""

def write_dat_files(rates, A, B, Peq):
    with open("in.A", "w") as fout:
        for a in A:
            fout.write("%d\n" % a)
    with open("in.B", "w") as fout:
        for b in B:
            fout.write("%d\n" % b)
    with open("in.Peq", "w") as fout:
        for u, p in Peq.iteritems():
            fout.write("%d %.16g\n" % (u, p))
    with open("in.rates", "w") as fout:
        for (u, v), p in rates.iteritems():
            fout.write("%d %d %.16g\n" % (u, v, p))
        

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
    
    write_dat_files(rate_constants, A, B, Peq)
    
    print "computing mean first passage times"
    calculator = NGT(rate_constants, A, B, weights=Peq)#, check_rates=False)
    calculator.compute()
    print "computing rates"
    kAB = calculator.get_rate_AB()
    print "k(B<-A)", kAB / knorm
    kBA = calculator.get_rate_BA()
    print "k(A<-B)", kBA / knorm
    
#    if True:
#        fname = "out.mfpt"
#        print "saving mean first passage times for all minima to reach B to file", fname
#        mfpt = sorted([(m, t) for m, t in 
#                       calculator.mfpt_computer.mfpt_dict.iteritems()])
#        with open(fname, "w") as fout:
#            for i, t in mfpt:
#                fout.write("%d %.12g\n" % (i, t * knorm))
#
#    print "computing committor probabilities"
#    sys.stdout.flush()
#    calculator.compute_committors()
#    print "computing steady state rates"
#    kSS = calculator.get_rate_AB_SS() / knorm
#    print "kSS(B<-A)", kSS
#    
#    if True:
#        fname = "out.committors"
#        print "saving committor probabilities for all minima to file", fname
#        coms = sorted([(m, c) for m, c in 
#                       calculator.committor_computer.committor_dict.iteritems()])
#        with open(fname, "w") as fout:
#            for i, c in coms:
#                fout.write("%d %.12g\n" % (i, c))
#    
    time_solve = calculator.time_solve #+ calculator.committor_computer.time_solve
    print "time spent solving linear equations", time_solve, "seconds"
    print "total time", time.clock() - tstart
     


if __name__ == "__main__":
    main()
