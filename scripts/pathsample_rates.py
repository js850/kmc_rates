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

from kmc_rates.rates_linalg import TwoStateRates, reduce_rates
from kmc_rates.utils import read_minA, make_rates, analyze_graph_error



description="""Compute rates from a pathsample database.  The transition
states and minima data are read from min.data and ts.data.  The product
and reactant states are read from min.A and min.B.  The rates are computed
by solving a set of linear equations using a sparse solver.
"""

def run(directory, T, A, B, out_prefix, tstart, reverse=False):
    if not reverse:
        source = "A"
        destination = "B"
        direction = "A->B"
    else:
        source = "B"
        destination = "A"
        direction = "B->A"

    rate_constants, Peq, knorm = make_rates(directory, T)
            

    if True:
        fname = "{}.rate_consts".format(out_prefix)
        print "saving rate constants to:", fname
        with open(fname, "w") as fout:
            fout.write("#starting_minimum ending_minimum rate_constant\n")
            for (u,v), k in sorted(rate_constants.iteritems()):
                fout.write("{} {} {}\n".format(u, v, k/knorm))

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
    print "k({}) {}".format(direction, kAB / knorm)
    
    if True:
        fname = "{}.rates".format(out_prefix)
        print "saving rates and mean first passage times for all minima to reach {} to file {}".format(destination, fname)
        mfpt = sorted([(m, t) for m, t in 
                       calculator.mfpt_computer.mfpt_dict.iteritems()])
        with open(fname, "w") as fout:
            fout.write("#Rates and mean first passage times for each node to reach {}\n".format(destination))
            fout.write("#The rate is just the inverse of the mean first passage time\n")
            fout.write("#minimum_index rate mean_first_passage_time\n")
            for i, t in mfpt:
                mt = t * knorm
                fout.write("{index} {rate} {mfpt}\n".format(index=i, rate=1./mt,
                                                            mfpt=mt))

    print "computing committor probabilities"
    sys.stdout.flush()
    calculator.compute_committors()
    print "computing steady state rate"
    kSS = calculator.get_rate_AB_SS() / knorm
    print "kSS({}) {}".format(direction, kSS)
    
    # print the committors
    # get the alternate definition of committors for the nodes in A and B
    Acomm = calculator.get_alternate_committors(A, B)
    Bcomm = calculator.get_alternate_committors(B, A)
    all_committors = calculator.committor_computer.committor_dict.copy()
    all_committors.update(Acomm)
    all_committors.update(Bcomm)
    
    fname = "{}.committors".format(out_prefix)
    print "saving committor probabilities for all minima to file", fname
    coms = sorted([(m, c) for m, c in all_committors.iteritems()])
    with open(fname, "w") as fout:
        fout.write("#probability a trajectory starting from each node ends up "
                   "in {B} before returning to {A}\n".format(B=destination, A=source))
        fout.write("#minimum_index committor_probability\n")
        for i, c in coms:
            fout.write("{} {}\n".format(i, c))
    
    time_solve = calculator.mfpt_computer.time_solve + calculator.committor_computer.time_solve
    print "time spent solving linear equations", time_solve, "seconds"
    print "total time", time.clock() - tstart


def main():
    """the main loop"""
    tstart =  time.clock()
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("-d", type=str, default=".",
                        help="directory with the min.data, ts.data files")
    parser.add_argument("-T", type=float, default=1., help="Temperature")
    parser.add_argument("--prefix", type=str, default="out", help="output prefix")
    args = parser.parse_args()
    
    print "temperature", args.T
    directory = args.d
    print "reading from directory:", os.path.abspath(directory)
    out_prefix = args.prefix
    
    A = read_minA(directory+"/min.A")
    B = read_minA(directory+"/min.B")
    
    print "\ncomputing rates from A to B"
    run(directory, args.T, A, B, out_prefix + ".AtoB", tstart, reverse=False)
    print "\ncomputing rates from B to A"
    run(directory, args.T, B, A, out_prefix + ".BtoA", tstart, reverse=True)
    
    
     


if __name__ == "__main__":
    main()
