import networkx as nx
import numpy as np

from kmc_rates import GraphReduction, kmcgraph_from_rates, KineticMonteCarlo

def main(plot=True):
    # make a graph representing a lattice in two dimensions
    dim = [10,10]
    grid_graph = nx.grid_graph(dim, periodic=False)
    A = [(0,0)]
    B = [(dim[0]-1, dim[1]-1)]
    
    # make some random rates for each of the edges
    rates = dict()
    for u, v in grid_graph.edges_iter():
        rates[(u,v)] = np.random.rand()
        rates[(v,u)] = np.random.rand()
        
    
    # set up the graph reduction object
    reducer = GraphReduction(rates, A, B)
    
    
    # make the kmc graph from the rates
    kmc_graph = kmcgraph_from_rates(rates)
    com_prob = reducer.compute_committor_probabilities(kmc_graph.nodes())

    if plot:    
        # put it into matrix form and plot
        P = np.zeros(dim)
        for node, p in com_prob.iteritems():
            x, y = node
            P[x,y] = p
        import matplotlib.pyplot as plt
        
        plt.imshow(P, cmap="BrBG", vmin=0, vmax=1)
        plt.title("probability of a trajectory reaching the lower right corner\n before reaching the upper left")
        plt.colorbar()
        plt.show()

if __name__ == "__main__":
    main()