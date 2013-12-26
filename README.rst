Kinetic Monte Carlo Rates
=========================

This package contains tools to quickly (and exactly) compute transition rates,
first passage times and and commitor probabilities in a transition network.
The resulting rates are exact, in the sense of Kinetic Monte Carlo, but the
analysis can be orders of magnitude faster than doing a Kinetic Monte Carlo
simulation.

getting started
---------------
This package uses external packages numpy and networkx. Those can be installed
from the python package index.
The following few lines of code creates a random graph with 4 nodes and 
then uses the graph renormalization method to compute the exact rates.::

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
    rAB, rBA = reducer.compute_rates()
    print "the transition rate from nodes", A, "to nodes", B, "is", rAB

See the examples files for more examples, including comparisons to kinetic
Monte Carlo simulations,


Description of the method
-------------------------
The rates are computed using the New Graph Transformation (NGT) method of
described in the paper

Calculating rate constants and committor probabilities for transition networks
by graph transformation
David Wales (2009) J. Chem. Phys., 130, 204111 
http://dx.doi.org/10.1063/1.3133782

The method uses a graph renormalization method (renormalization in the sense of
renormalization group theory) to compute exact Kinetic Monte Carlo rates and
first passage probabilities from a reactant group A to a product group B.  Each
node `u` has an attribute `tau_u` which is the waiting time at that node.  Each
edge `u -> v` has an associated transition probability and `P_uv`.  An
important feature of this algorithm is that each node has a loop edge pointing
back to itself and associated probability `P_uu` which gives the self-transitio
probability.  In the typical case the self-transition probabilities will all be
zero initially, but will take non zero values after renormalization.  The
transition probabilities always satisfy `sum_v P_uv = 1`.

graph transformation
++++++++++++++++++++

The algorithm is most easily described if we first assume that A, and B each
contain only one node (called `a`, and `b`).  In the algorithm, nodes are
iteratively removed from the graph until the only two remaining nodes are a and
`b`.  Upon removing node `x`, in order to preserve the transition times and
probabilities, the properties of the neighbors of `x` are all updated.  For
each neighbor, `u`, of `x`, the transition times are updated according to

    tau_u -> tau_u + P_ux * tau_x / (1 - P_xx)

Similarly, for each pair, `u` and `v`, of neighbors of `x`
the transition probabilities are updated according to 

    P_uv -> P_uv + P_ux * P_xv / (1 - P_xx)

Note that the self-transition probabilities `P_uu` are also updated according to the
above equation.

Once the graph is reduced to only the two nodes, `a`, and `b`, 
the probability `P_ab` is interpreted as the commitor probability from `a` to `b`.  
That is, the probability that a trajectory starting at `a` will end up at `b` before returning to `a`.  
Similarly, the mean first passage time from `a` to `b` is simply
`tau_a / P_ab`.  Note that because the probabilies sum to 1 the mean first passage time can also
be written `tau_a / (1-P_aa)`.
The transtion rate from `a` to `b`
is simply the inverse of the mean first passage time.  The rates and
probabilites from `b -> a` are read from the resuling graph in the same way.
The above interpretations are exact in the sense that a Kinetic Monte Carlo
simulation will give the same result.


if B has more than one element
++++++++++++++++++++++++++++++

If there is more than one element in `B` the calculation of rates from `a -> B`
is nearly as simple.  Following the same procedure described above, all nodes
except those in `A` or `B` are iteratively removed.  
The commitor probability from `a` to `B` is then the sum over the transition probabilities
from `a` to `b` for each element `b` in `B`.  This can also be written as

    1 - P_aa

The mean first passage time from `a` to `B` is given by

    T_aB = tau_a / (1 - P_aa)

if A has more than one element
++++++++++++++++++++++++++++++

In this, the most general case, when both A and B have more than one element,
the transition rate from `A` to `B` must be computed as an average over
the inverse mean first passage time for each
element `a` in `A`. That is

  k_AB = average( 1 / T_aB )

The computation is done in two phases.  In the first phase the intermediate
nodes (those not in `A` or in `B`) are all removed from the graph.  In the
second phase we first make a backup copy of the graph.  Then for each node `a`
in `A` we remove from the graph all nodes in `A` (except `a`). This allows us
to compute commitor probabilities and mean first passage times (`T_aB`) from
`a` to `B` as described in the preceding section.

If the nodes `a` are not all equally likely to be occupied, then the above
average can be a weighted average where each node is weighted according to its
equilibrium occupation probabilities.

The rates `B -> A` can be computed in a similar manner
