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
    Graph G;

    graph_ns::BuildGridGraph2d builder(G, 3, 4, true);
    builder.build_graph();


    cout << "number of nodes " << G.number_of_nodes() << "\n";
    cout << "number of edges " << G.number_of_edges() << "\n";

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

    cout << "deleting node 0\n";
    G.remove_node(0);
    cout << "number of nodes " << G.number_of_nodes() << "\n";
    cout << "number of edges " << G.number_of_edges() << "\n";
    cout << "deleting node 1\n";
    G.remove_node(1);
    cout << "number of nodes " << G.number_of_nodes() << "\n";
    cout << "number of edges " << G.number_of_edges() << "\n";

    cout << "\n\nmaking a copy of the graph\n";
    Graph G2(G);
    BreadthFirstSearchEdges bfs2(&G2, 2);
    cout << "number of nodes " << G2.number_of_nodes() << "\n";
    cout << "number of edges " << G2.number_of_edges() << "\n";
    e = bfs2.get_next_edge();
    while (e != NULL){
        node_id u = e->tail()->id();
        node_id v = e->head()->id();
        cout << u << " -> " << v << "\n";
        e = bfs2.get_next_edge();
    }


    
}
