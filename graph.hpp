#ifndef _GRAPH_CPP_H_
#define _GRAPH_CPP_H_

#include <cstdlib>
#include <iostream>
#include <list>
#include <map>

namespace graph_ns
{
class Edge;
class Node;
typedef size_t node_id;
typedef Node * node_ptr;
typedef Edge * edge_ptr;
typedef int color_type;
color_type color_white = 0;
color_type color_grey = 1;
color_type color_black = 4;

/**
 * basic class for a node in the graph
 */
class Node{
public:
    typedef std::list<Edge *> out_edge_list;
private:
    node_id id_;
    out_edge_list out_edge_list_; // list of outgoing edges
public:

    //Edge * first_outgoing_; // first outgoing edge

    Node(node_id id):
        id_(id)
    {}

    void add_out_edge(Edge * edge);
    out_edge_list & get_out_edges(){ return out_edge_list_; }
    node_id id() const { return id_; }
};

/**
 * basic class for an edge (arc) in the graph
 */
class Edge{
public:
    node_ptr head_; // node the edge points to
    node_ptr tail_; // node the edge comes from

    Edge(node_ptr tail, node_ptr head):
        head_(head),
        tail_(tail)
    {}

    node_ptr head(){ return head_; }
    node_ptr tail(){ return tail_; }
};

void Node::add_out_edge(Edge * edge)
{
    out_edge_list_.push_back(edge);
}


class Graph
{
    std::map<node_id, node_ptr> node_map_;
    std::list<Edge *> edge_list_;

    node_id next_node_id_;

public:
    Graph():
        next_node_id_(0)
    {}

    ~Graph()
    {
        // delete all nodes
        for (auto & mapval : node_map_){
            node_ptr node = mapval.second;
            delete node;
        }
        node_map_.clear();

        // delete all edges
        for (auto & edge : edge_list_){
            delete edge;
        }
        edge_list_.clear();

    }

    node_id number_of_nodes() { return node_map_.size(); }
    node_id number_of_edges() { return edge_list_.size(); }

    /**
     * create a new node
     */
    node_id add_node(){
        node_ptr node = new Node(next_node_id_);
        next_node_id_++;
        node_map_.insert(std::pair<node_id, node_ptr> (node->id(), node));
        return node->id();
    }

    /**
     * create a n new nodes
     */
    void add_nodes(node_id n){
        for (node_id i = 0; i < n; ++i){
            add_node();
        }
    }

    /**
     * return a pointer to the node with given node id
     */
    node_ptr get_node(node_id nodeid)
    {
        typedef std::map<node_id, node_ptr> maptype;
        maptype::iterator iter = node_map_.find(nodeid);
        if (iter == node_map_.end()){
            return NULL;
        }
        return iter->second;
        return NULL;
    }

    /**
     * add an edge from tail to head
     */
    void add_edge(node_id tail, node_id head)
    {
        node_ptr node_tail = get_node(tail);
        node_ptr node_head = get_node(head);
        Edge * edge = new Edge(node_tail, node_head);
        edge_list_.push_back(edge);
        node_tail->add_out_edge(edge);
    }

};


}

#endif
