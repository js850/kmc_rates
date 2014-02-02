import argparse
import sys
import os
import time
from itertools import izip
import numpy as np
import networkx as nx


from pele.utils.optim_compatibility import OptimDBConverter
from pele.storage import Database
from pele.rates import RateCalculation

from rates_linalg import TwoStateRates, reduce_rates


def log_sum2(a, b):
    """
    return log( exp(a) + exp(b) )
    """
    if a > b:
        return a + np.log(1.0 + np.exp(-a + b) )
    else:
        return b + np.log(1.0 + np.exp(a - b) )

def compute_log_rates(m1vals, tsdata, T):
    menergy = m1vals[:,0]
    mfvib = m1vals[:,1]
    mpgorder = m1vals[:,2]
    tsenergy = tsdata[:,0]
    tsfvib = tsdata[:,1]
    tspgorder = tsdata[:,2]
    return (np.log(mpgorder / (2. * np.pi * tspgorder))
            + (mfvib - tsfvib)/2. 
            - (tsenergy - menergy) / T)



class RatesFromPathsampleDB(object):
    def __init__(self, mindata_file, tsdata_file, T):
        self.T = T
        self.run(mindata_file, tsdata_file)
    
    def _log_equilibrium_occupation_probabilities(self, mindata):
        """return the log equilibrium occupation probability
        
        This is computed from the harmonic superposition approximation.  Some
        constants that are the same for all minima have been left out.
        """
        energy = mindata[:,0]
        fvib = mindata[:,1]
        pgorder = mindata[:,2]
        return (-self.T * energy - np.log(pgorder) - 0.5 * fvib)
    
    def compute_rates(self, mindata, tsdata, mindices):
        """return the array of rates over the transition states
        """
        m1vals = mindata[mindices,:]
        
        assert m1vals.shape[0] == tsdata.shape[0]
        return compute_log_rates(m1vals, tsdata, self.T)
    
    def make_rates_dict(self, m1_indices, m2_indices, log_kuv_list, log_kvu_list):
        """turn the arrays into a dictionary of rates"""
        log_rates = dict()
        
        
        for u, v, log_kuv, log_kvu in izip(m1_indices, m2_indices, 
                                           log_kuv_list, log_kvu_list):
            if u == v:
                # don't add transition states from a minimum to itself
                continue
            if (u,v) in log_rates:
                log_rates[(u,v)] = log_sum2(log_rates[(u,v)], log_kuv)
                log_rates[(v,u)] = log_sum2(log_rates[(v,u)], log_kvu)
            else:
                log_rates[(u,v)] = log_kuv
                log_rates[(v,u)] = log_kvu
        
        rates = dict(( (uv,np.exp(k)) for uv, k in log_rates.iteritems() ))
        return rates
                
    
    def run(self, mindataf, tsdataf):
        # mindata format
        # energy fvib pgorder I1 I2 I3
        print "reading minima data from file:", mindataf
        mindata = np.genfromtxt(mindataf)
#         mindata = mindata[:,:3]
        
        # tsdata format
        # energy fvib pgorder min1 min2 I1 I2 I3
        print "reading transition state data from file:", tsdataf
        tsdata = np.genfromtxt(tsdataf)
#         tsdata = tsdata[:,:5]
        
        m1_indices = tsdata[:,3].astype(int) - 1
        m2_indices = tsdata[:,4].astype(int) - 1
        
        print "computing rate constants"
        log_k12 = self.compute_rates(mindata, tsdata, m1_indices)
        log_k21 = self.compute_rates(mindata, tsdata, m2_indices)

        max_log_rate = max(log_k12.max(), log_k21.max())
        log_k12 -= max_log_rate
        log_k21 -= max_log_rate
        self.log_rate_subtracted = max_log_rate
        self.rate_norm = np.exp(-max_log_rate)
        
        self.rate_constants = self.make_rates_dict(m1_indices, m2_indices, log_k12, log_k21)
        
        print "computing equilibrium occupation probabilities"
        log_Peq = self._log_equilibrium_occupation_probabilities(mindata)
        Peq = np.exp(log_Peq - log_Peq.max())
        self.Peq = dict(( (u, P) for u, P in enumerate(Peq) ))
        

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
    
     

def read_minA(fname, db=None):
    with open(fname) as fin:
        ids = []
        for i, line in enumerate(fin):
            if i == 0:
                nminima = int(line.split()[0])
            else:
                sline = line.split()
                ids += map(int, sline)
    
    assert nminima == len(ids)
    if db is None:
        return ids
    else:
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

def make_rates_pele(database, A, B, T):
    print "computing rates from transition states"
    pele_rates = RateCalculation(database.transition_states(), A, B, T=T, use_fvib=True)
    pele_rates._make_kmc_graph()
    rates = pele_rates.rate_constants
    print "computing equilibrium occupation probabilities"
    Peq = pele_rates._get_equilibrium_occupation_probabilities()
    
    return rates, Peq, np.exp(-pele_rates.max_log_rate)

def make_rates(directory, T):
    print "reading from directory:", os.path.abspath(directory)
    mindata = directory + "/min.data"
    tsdata = directory + "/ts.data"
    generator = RatesFromPathsampleDB(mindata, tsdata, T)
    
    return generator.rate_constants, generator.Peq, generator.rate_norm


def main():
    tstart =  time.clock()
    parser = argparse.ArgumentParser(description="compute rates from a pathsample database")

    parser.add_argument("-d", type=str, default=".",
                        help="directory with the min.data, ts.data files")
    parser.add_argument("-T", type=float, default=1., help="Temperature")
    args = parser.parse_args()
    print "temperature", args.T
    directory = args.d
    
#     if True:
#         reader = RatesFromPathsampleDB(args.T)
#         reader.parse_files("min.data.cf", "ts.data.cf")
#         exit()
    
    if True:
        A = read_minA(directory+"/min.A")
        B = read_minA(directory+"/min.B")
        print len(A), "A minima read from min.A"
        print len(B), "B minima read from min.B"
        rate_constants, Peq, knorm = make_rates(args.d, args.T)
    else:
        # load the database
        db, A, B = load_database(args.d)
        # compute rate constants
        rate_constants, Peq, knorm = make_rates_pele(db, A, B, args.T)

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
        mfpt = sorted([(m, t) for m, t in 
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