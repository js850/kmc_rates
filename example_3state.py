from kmc_rates import GraphReduction, graph_from_rates
from kmc import KineticMonteCarlo


def main():
    # make matrix of transition rates with three states
    transition_matrix = [ [0., 1., 1.], 
                          [1., 0., 1.], 
                          [1., 1., 0.] ]
    # an easy way to set up the graph is to make
    # a dictionary of rates
    rates = dict()
    for i in range(3):
        for j in range(3):
            if i != j:
                rates[(i,j)] = transition_matrix[i][j]

    # use the utility function to build the rate graph with the correct formatting
    rate_graph = graph_from_rates(rates)
    rate_graph_backup = rate_graph.copy()

    # we want to compute rates from node 0 to node 1
    A = [0]
    B = [1]
    
    # set up the class which will compute rates for us
    reducer = GraphReduction(rate_graph, A, B)
    
    # do the calculation
    rAB, rBA = reducer.compute_rates()
    
    print "For this transition network with 3 nodes the exact rate constant should be 1.0"
    print "the transition rate computed by graph transformation is", rAB
    
    # now to check the values do a kinetic monte carlo run
    kmc = KineticMonteCarlo(rate_graph_backup)
    rAB100 = kmc.mean_rate(A, B, niter=100)
    print "the KMC rate averaged over 100 iterations is", rAB100
    rAB1000 = kmc.mean_rate(A, B, niter=1000)
    print "the KMC rate averaged over 1000 iterations is", rAB1000
    rAB10000 = kmc.mean_rate(A, B, niter=10000)
    print "the KMC rate averaged over 10000 iterations is", rAB10000

if __name__ == "__main__":
    main()