import argparse
import sys
import os
import time
import numpy as np

from pele.utils.optim_compatibility import OptimDBConverter
from pele.storage import Database
from pele.rates import RateCalculation

from rates_linalg import TwoStateRates

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
    pele_rates = RateCalculation(database.transition_states(), A, B, T=T, use_fvib=True)
    pele_rates._make_kmc_graph()
    rates = pele_rates.rate_constants
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
    
    calculator = TwoStateRates(rate_constants, A, B, weights=Peq)
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