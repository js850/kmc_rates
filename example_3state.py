from kmc_rates import GraphReduction, kmcgraph_from_rates
from kmc import KineticMonteCarlo


def main():
    nnodes = 3
    # make matrix of transition rates with three states
    transition_matrix = [ [0., 1., 1.], 
                          [1., 0., 1.], 
                          [1., 1., 0.] ]
    print "Computing rates for transition matrix"
    print transition_matrix
    # an easy way to set up the graph is to make
    # a dictionary of rates
    rates = dict()
    for i in range(nnodes):
        for j in range(nnodes):
            if i != j:
                rates[(i,j)] = transition_matrix[i][j]

    # use the utility function to build the rate graph with the correct formatting
    rate_graph = kmcgraph_from_rates(rates)
    rate_graph_backup = rate_graph.copy()

    # we want to compute rates from node 0 to node 1
    A = [0]
    B = [1]
    
    # set up the class which will compute rates for us
    reducer = GraphReduction(rate_graph, A, B)
    
    # do the calculation
    print "computing the transition rate from nodes", A, "to nodes", B
    reducer.compute_rates()
    rAB = reducer.get_rate_AB()
    
    print "the exact rate for this network is 1.0"
    print "the transition rate computed by graph transformation is", rAB
    
    # now to check the values do a kinetic monte carlo run
    kmc = KineticMonteCarlo(rate_graph_backup)
    
    niterlist = [10, 100, 1000, 10000]
    for niter in niterlist:
        rAB_KMC = kmc.mean_rate(A, B, niter=niter)
        print "the KMC rate averaged over %8d iterations is %s. abs(rAB-rAB_KMC) = %s" % (niter, rAB_KMC, abs(rAB-rAB_KMC))

if __name__ == "__main__":
    main()
