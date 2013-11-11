#include "graph.hpp"
#include "breadth_first_search.hpp"
#include <iostream>

using std::cout;
using graph_ns::Graph;
using graph_ns::node_id;
using graph_ns::BreadthFirstSearch;
using graph_ns::edge_ptr;
using graph_ns::node_ptr;

int main()
{
    Graph G = Graph();
    G.add_nodes(10);
    G.add_edge(0, 1);
    G.add_edge(1, 2);
    G.add_edge(2, 3);
    G.add_edge(2, 4);
    G.add_edge(4, 1);

    cout << G.number_of_nodes() << "\n";
    cout << G.number_of_edges() << "\n";

    BreadthFirstSearch bfs(&G, 0);
    edge_ptr e = bfs.get_next_edge();
    while (e != NULL){
        node_id u = e->tail()->id();
        node_id v = e->head()->id();
        cout << u << " -> " << v << "\n";
        e = bfs.get_next_edge();
    }
    
}
