import numpy as np

from kmc_rates import GraphReduction, graph_from_rates
from kmc import KineticMonteCarlo

def readme_example():
    from kmc_rates import GraphReduction, graph_from_rates
    nnodes = 4
    # create a dictionary of transition rates
    rates = dict()
    for i in range(nnodes):
        for j in range(nnodes):
            if i != j:
                rates[(i,j)] = np.random.rand()
    
    # set up the calculation of the transition rate from node 0 to node 1
    A = [0]
    B = [1]
    kmc_graph = graph_from_rates(rates)
    reducer = GraphReduction(kmc_graph, A, B)
    reducer.compute_rates()
    rAB = reducer.get_rate_AB()
    print "the transition rate from nodes", A, "to nodes", B, "is", rAB
    
    

def main():
    nnodes = 6
    # the graph need not be made from a transition matrix, but it's an 
    # easy way to make a random graph ensuring that everything is connected
    transition_matrix = np.random.uniform(0,1,[nnodes, nnodes])
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
    kmc_graph = graph_from_rates(rates)
    kmc_graph_backup = kmc_graph.copy()

    # we want to compute rates from node 0 to node 1
    A = [0,1]
    B = [2,3]
    
    # set up the class which will compute rates for us
    reducer = GraphReduction(kmc_graph, A, B)
    reducer.compute_rates()
    rAB = reducer.get_rate_AB()
    
    # do the calculation
    print "computing the transition rate from nodes", A, "to nodes", B
    reducer = GraphReduction(kmc_graph, A, B)
    reducer.compute_rates()
    
    print "the transition rate computed by graph transformation is", rAB
    
    # now to check the values do a kinetic monte carlo run
    kmc = KineticMonteCarlo(kmc_graph_backup)
    
    niterlist = [10, 100, 1000, 10000]
    for niter in niterlist:
        rAB_KMC = kmc.mean_rate(A, B, niter=niter)
        print "the KMC rate averaged over %8d iterations is %s. abs(rAB-rAB_KMC) = %s" % (niter, rAB_KMC, abs(rAB-rAB_KMC))

if __name__ == "__main__":
#     readme_example()
    main()
