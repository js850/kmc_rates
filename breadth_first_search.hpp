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

//template <class Graph>
class BFSVisitor {
protected:
    typedef node_ptr node_t;
    typedef Edge * edge_t;
public:
    BFSVisitor(){}
    //void initialize_vertex(node_t u, Graph & g){}
    void discover_node(node_t u, Graph & g){}
    void examine_node(node_t u, Graph & g){}
    void examine_edge(edge_t e, Graph & g){}
    void tree_edge(edge_t e, Graph & g){}
    //void non_tree_edge(edge_t e, Graph & g){}
    //void gray_target(edge_t e, Graph & g){}
    //void black_target(edge_t e, Graph & g){}
    void finish_node(node_t u, Graph & g){}
};

template<class bfsvisitor>
class BreadthFirstSearch
{
    Graph & graph_;
    std::vector<color_type> node_color;
    std::queue<node_ptr> node_queue_;

    unsigned long Nedges = 0;
    unsigned long Nnodes = 0;

    bfsvisitor & vis;

public:
    BreadthFirstSearch(Graph & graph, node_id n, bfsvisitor & visitor):
        graph_(graph),
        node_color(graph_.number_of_nodes(), color_white),
        vis(visitor)
    {
        // add the first node to the queue
        assert(n < graph_.number_of_nodes());
        node_color[n] = color_grey;
        node_ptr node = graph_.get_node(n);
        node_queue_.push(node);
        vis.discover_node(node, graph_);
        ++Nnodes;
    }

    void run(){
        while (! node_queue_.empty() ){
            auto node = node_queue_.front();             
            vis.examine_node(node, graph_);            
            node_queue_.pop();
            for (auto edge: node->get_out_edges()){      
                vis.examine_edge(edge, graph_);
                auto u = edge->head();
                if (node_color[u->id()] == color_white){ 
                    vis.tree_edge(edge, graph_);
                    vis.discover_node(u, graph_);
                    node_color[u->id()] = color_grey;    
                    node_queue_.push(u);
                }
            }
            node_color[node->id()] = color_black;
            vis.finish_node(node, graph_);
        }
    }
};


class BreadthFirstSearchEdges
{
    Graph * graph_;
    std::vector<color_type> node_color;
    std::queue<node_ptr> node_queue_;

    // these hold the current state of the iteration
    Node::out_edge_list::iterator edge_iter_;
    Node::out_edge_list * edge_list = NULL;

    unsigned long Nedges = 0;
    unsigned long Nnodes = 0;

public:
    BreadthFirstSearchEdges(Graph * graph, node_id n):
        graph_(graph),
        node_color(graph_->number_of_nodes(), color_white)
    {
        // add the first node to the queue
        assert(n < graph_->number_of_nodes());
        node_color[n] = color_grey;
        node_ptr node = graph_->get_node(n);
        node_queue_.push(node);
        ++Nnodes;
    }

    Edge * get_next_edge()
    {
        if ( node_queue_.empty() ){
            return NULL;
        } else if (Nedges == 0){
            node_ptr u = node_queue_.front();
            edge_list = &(u->get_out_edges());
            edge_iter_ = edge_list->begin();
        } else if (edge_iter_ == edge_list->end()){
            // done with this node.  time to go on to the next node.
            // pop the old node and color it black
            node_ptr uold = node_queue_.front();
            node_queue_.pop();
            assert(node_color[uold->id()] == color_grey);
            node_color[uold->id()] = color_black;

            // get the next node from the front of the queue
            if (node_queue_.empty()){
                edge_iter_ = edge_list->end();
            } else {
                node_ptr u = node_queue_.front();
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
            node_ptr v = (*edge_iter_)->head();
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
