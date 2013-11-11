#include "graph.hpp"
#include <iostream>

using std::cout;
using graph_ns::Graph;
using graph_ns::node_id;

int main()
{
    Graph G = Graph();
    node_id u1 = G.add_node();
    node_id u2 = G.add_node();
    G.add_edge(u1, u2);

    cout << G.number_of_nodes() << "\n";
    cout << G.number_of_edges() << "\n";
    
}
