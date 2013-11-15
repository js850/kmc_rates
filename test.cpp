#include "graph.hpp"
#include "breadth_first_search.hpp"
#include "grid_graph.hpp"
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

    auto builder = graph_ns::BuildGridGraph2d(G, 3, 4, true);
    builder.build_graph();


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
