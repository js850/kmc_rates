#ifndef _BREADTH_FIRST_SEARCH_HPP_
#define _BREADTH_FIRST_SEARCH_HPP_

#include <cstdlib>
#include <iostream>
#include <list>
#include <queue>
#include <assert.h>

#include "graph.hpp"

namespace graph_ns
{
class BreadthFirstSearch
{
    Graph * graph_;
    std::vector<color_type> node_color;
    std::queue<Node *> node_queue_;

    // these hold the current state of the iteration
    Node::out_edge_list::iterator edge_iter_;
    Node::out_edge_list * edge_list;

    unsigned long Nedges;
    unsigned long Nnodes;

public:
    BreadthFirstSearch(Graph * graph, node_id n):
        graph_(graph),
        node_color(graph_->number_of_nodes(), color_white),
        edge_list(NULL),
        Nedges(0)
    {
        // add the first node to the queue
        assert(n < graph_->number_of_nodes());
        node_color[n] = color_grey;
        Node * node = graph_->get_node(n);
        node_queue_.push(node);
        ++Nnodes;
    }

    Edge * get_next_edge()
    {
        if ( node_queue_.empty() ){
            return NULL;
        } else if (Nedges == 0){
            Node *u = node_queue_.front();
            edge_list = &(u->get_out_edges());
            edge_iter_ = edge_list->begin();
        } else if (edge_iter_ == edge_list->end()){
            // done with this node.  time to go on to the next node.
            // pop the old node and color it black
            Node * uold = node_queue_.front();
            node_queue_.pop();
            assert(node_color[uold->id()] == color_grey);
            node_color[uold->id()] = color_black;

            // get the next node from the front of the queue
            if (node_queue_.empty()){
                edge_iter_ = edge_list->end();
            } else {
                Node * u = node_queue_.front();
                edge_list = &(u->get_out_edges());
                edge_iter_ = edge_list->begin();
            }
        } else {
            // got to next out edge from this node
            ++edge_iter_;
        }

        if (edge_iter_ == edge_list->end()){
            return get_next_edge();
        } else {
            // examine the head of the current edge if it's color is white then
            // color it grey and add it to the queue
            Node * v = (*edge_iter_)->head();
            assert(v!=NULL);
            assert(v->id() < node_color.size());
            color_type c = node_color[v->id()];
            if (c == color_white){
                node_color[v->id()] = color_grey;
                node_queue_.push(v);
                ++Nnodes;
            }
        }
        ++Nedges;
        return *edge_iter_;

    }
};
}
#endif
