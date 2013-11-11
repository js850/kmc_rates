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
    Edge * current_edge_;
    std::queue<Node *> node_queue_;
    bool start_;

public:
    BreadthFirstSearch(Graph * graph, node_id n):
        graph_(graph),
        node_color(graph_->number_of_nodes(), color_white),
        current_edge_(NULL),
        start_(true)
    {
        assert(n < graph_->number_of_nodes());
        node_color[n] = color_grey;
        Node * node = graph_->get_node(n);
        node_queue_.push(node);
    }

    Edge * get_next_edge()
    {
        Edge * next_edge = NULL;
        if ( node_queue_.empty() ){
            return NULL;
        } else if (start_){
            start_ = false;
            Node *u = node_queue_.front();
            next_edge = u->first_out_edge();
        } else if (current_edge_ == NULL){
            // done with this node.  time to go on to the next node.
            // pop the old node and color it black
            Node * uold = node_queue_.front();
            node_queue_.pop();
            assert(node_color[uold->id()] == color_grey);
            node_color[uold->id()] = color_black;

            // get the next node from the front of the queue
            if (node_queue_.empty()){
                next_edge = NULL;
            } else {
                Node * v = node_queue_.front();
                next_edge = v->first_out_edge();
            }
        } else {
            // got to next out edge from this node
            next_edge = current_edge_->next_edge();
        }

        current_edge_ = next_edge;
        if (next_edge == NULL){
            return get_next_edge();
        } else {
            // examine the head of the current edge if it's color is white then
            // color it grey and add it to the queue
            Node * v = next_edge->head();
            assert(v!=NULL);
            assert(v->id() < node_color.size());
            color_type c = node_color[v->id()];
            if (c == color_white){
                node_color[v->id()] = color_grey;
                node_queue_.push(v);
            }
        }
        return next_edge;

    }
};
}
#endif
