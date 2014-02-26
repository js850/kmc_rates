from kmc_rates import GraphReduction, KineticMonteCarlo
from kmc_rates._kmc_rates import kmcgraph_from_rates


def main():
    nnodes = 3
    # make matrix of transition rates with three states
    transition_matrix = [ [0., 1., 1.], 
                          [1., 0., 1.], 
                          [1., 1., 0.] ]
    print "Computing rates and committor probabilities for transition matrix"
    print transition_matrix
    # an easy way to set up the graph is to make
    # a dictionary of rates
    rates = dict()
    for i in range(nnodes):
        for j in range(nnodes):
            if i != j:
                rates[(i,j)] = transition_matrix[i][j]


    # we want to compute rates from node 0 to node 1
    A = [0]
    B = [1]
    x = 0
    
    # set up the class which will compute rates for us
    reducer = GraphReduction(rates, A, B)
    
    
    # now to check the values do a kinetic monte carlo run
    # use the utility function to build the rate graph with the correct formatting
    kmc_graph = kmcgraph_from_rates(rates)
    kmc = KineticMonteCarlo(kmc_graph)
    
    niterlist = [10, 100, 1000, 10000]

    print ""
    print "computing the probability for a trajectory to start at", x, "and reach", B, "before reaching", A
    PxB = reducer.compute_committor_probability(x)
    print "the committor probability computed by graph transformation is", PxB
    for niter in niterlist:
        PxB_KMC = kmc.committor_probability(x, A, B, niter=niter)
        print "the KMC committor probability averaged over %8d trajectories is %s. abs(PxB-PxB_KMC) = %s" % (niter, PxB_KMC, abs(PxB-PxB_KMC))


    print ""
    print "computing the transition rate from nodes", A, "to nodes", B
    print "the exact rate for this network is 1.0"
    reducer.compute_rates()
    rAB = reducer.get_rate_AB()
    print "the transition rate computed by graph transformation is", rAB
    for niter in niterlist:
        rAB_KMC = kmc.mean_rate(A, B, niter=niter)
        print "the KMC rate averaged over %8d trajectories is %s. abs(rAB-rAB_KMC) = %s" % (niter, rAB_KMC, abs(rAB-rAB_KMC))

if __name__ == "__main__":
    main()
