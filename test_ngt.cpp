#include "graph.hpp"
#include "breadth_first_search.hpp"
#include "grid_graph.hpp"
#include "ngt.hpp"
#include <iostream>

using std::cout;
using namespace graph_ns;
using graph_ns::Graph;
using graph_ns::node_id;
using graph_ns::edge_ptr;
using graph_ns::node_ptr;



int main()
{
    std::vector<node_id> A, B;
    A.push_back(0);
    B.push_back(1);

    typedef std::pair<node_id, node_id> pair_t;
    NGT::rate_map_t rates;
    rates[pair_t(0,1)] = 1.;
    rates[pair_t(1,0)] = 1.;
    rates[pair_t(0,2)] = 1.;
    rates[pair_t(2,0)] = 1.;
    rates[pair_t(1,2)] = 1.;
    rates[pair_t(2,1)] = 1.;

    NGT ngt(rates, A, B);
    ngt.debug = true;
    node_ptr a = ngt._graph->get_node(0);
    node_ptr b = ngt._graph->get_node(1);
    edge_ptr ab = a->get_successor_edge(b);
    edge_ptr ba = b->get_successor_edge(a);

    double tau_a = ngt.get_tau(a);
    double Paa = ngt.get_node_P(a);
    double Pab = ngt.get_P(ab);
    cout << "a tau " << tau_a << "\n";
    cout << "a P " << Paa << "\n";
    cout << "ab P " << Pab << "\n";

    ngt.remove_intermediates();


    tau_a = ngt.get_tau(a);
    Paa = ngt.get_node_P(a);
    Pab = ngt.get_P(ab);
    cout << "tau_a " << tau_a << "\n";
    cout << "Paa " << Paa << "\n";
    cout << "Pab " << Pab << "\n";
    cout << "rate A -> B " << ngt.get_rate_AB() << "\n";
    cout << "rate B -> A " << ngt.get_rate_BA() << "\n";
}
