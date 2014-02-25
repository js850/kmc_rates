from ngt_wrapper import NGT

from kmc_rates import GraphReduction, kmcgraph_from_rates

print "hello"

def main():
    nnodes = 100
    # an easy way to set up the graph is to make
    # a dictionary of rates
    rates = dict()
    for i in range(nnodes):
        for j in range(i+1,nnodes):
            rates[(i,j)] = float(i+j) / (i+1)
            rates[(j,i)] = float(i+j) / (j+1)

    # we want to compute rates from node 0 to node 1
    A = [0 ,1, 2]
    B = [3, 4]
    x = 2
    
    weights=dict([(a,1.) for a in A+B])
    weights[2] = 5e3


    ngt = NGT(rates, A, B)
    ngt.compute()
    kAB = ngt.get_rate_AB()
    kBA = ngt.get_rate_BA()
    print "rate AB", kAB
    print "rate BA", kBA
    
    graph = kmcgraph_from_rates(rates)
    pyngt = GraphReduction(graph, A, B)
    pyngt.compute_rates()
    print "pyngt rate AB", pyngt.get_rate_AB()
    print "pyngt rate AB", pyngt.get_rate_BA()

if __name__ == "__main__":
    main()