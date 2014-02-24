#include "graph.hpp"
#include "breadth_first_search.hpp"
#include "grid_graph.hpp"
#include <iostream>

using std::cout;
using namespace graph_ns;
using graph_ns::Graph;
using graph_ns::node_id;
using graph_ns::edge_ptr;
using graph_ns::node_ptr;

class print_visitor : public BFSVisitor
{
public:
    void examining_node(node_t u, Graph & g){
        cout << "    examining node " << u->id() << "\n";
    }
    void discover_node(node_t u, Graph & g){
        cout << "discovered node " << u->id() << "\n";
    }
    void tree_edge(edge_t e, Graph & g){
        cout << "discovered edge " << e->tail()->id() << " -> " << e->head()->id() << "\n";
    }
    void examine_edge(edge_t e, Graph & g){
        cout << "    examining edge " << e->tail()->id() << " -> " << e->head()->id() << "\n";
    }
};


int main()
{
    Graph G = Graph();

    graph_ns::BuildGridGraph2d builder = graph_ns::BuildGridGraph2d(G, 3, 4, true);
    builder.build_graph();


    cout << G.number_of_nodes() << "\n";
    cout << G.number_of_edges() << "\n";

    BreadthFirstSearchEdges bfs(&G, 0);
    edge_ptr e = bfs.get_next_edge();
    while (e != NULL){
        node_id u = e->tail()->id();
        node_id v = e->head()->id();
        cout << u << " -> " << v << "\n";
        e = bfs.get_next_edge();
    }

    cout << "Now runing bfs visitor\n";
    print_visitor visitor;
    BreadthFirstSearch<print_visitor> bfs1(G, 0, visitor);
    bfs1.run();
    
}
