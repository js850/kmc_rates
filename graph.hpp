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
    node_id id_;
    Edge * first_outgoing_; // first outgoing edge

    Node(node_id id):
        id_(id),
        first_outgoing_(NULL)
    {}

    void add_out_edge(Edge * edge);
    Edge * first_out_edge(){ return first_outgoing_; }
    node_id id() const { return id_; }
};

/**
 * basic class for an edge (arc) in the graph
 */
class Edge{
public:
    Node * head_; // node the edge points to
    Node * tail_; // node the edge comes from
    Edge * next_; // next edge with same originating node

    Edge(Node * tail, Node * head):
        head_(head),
        tail_(tail),
        next_(NULL)
    {}

    Edge * next_edge(){ return next_; }
    Node * head(){ return head_; }
    Node * tail(){ return tail_; }
};

void Node::add_out_edge(Edge * edge)
{
    edge -> next_ = first_outgoing_;
    first_outgoing_ = edge;
}


class Graph
{
    std::list<Node *> node_list_;
    std::list<Edge *> edge_list_;

    node_id next_node_id_;

public:
    Graph():
        next_node_id_(0)
    {}

    ~Graph()
    {
        // delete all the nodes and edges
        for (std::list<Node *>::iterator iter = node_list_.begin();
                iter != node_list_.end(); ++iter){
            delete *iter;
        }
        for (std::list<Edge *>::iterator iter = edge_list_.begin();
                iter != edge_list_.end(); ++iter){
            delete *iter;
        }
    }

    node_id number_of_nodes() { return node_list_.size(); }
    node_id number_of_edges() { return edge_list_.size(); }

    /**
     * create a new node
     */
    node_id add_node(){
        Node *node = new Node(next_node_id_);
        next_node_id_++;
        node_list_.push_back(node);
        return node->id_;
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
    Node * get_node(node_id node_id)
    {
        for (std::list<Node *>::iterator iter = node_list_.begin();
                iter != node_list_.end(); ++iter){
            if ((*iter)->id_ == node_id){
                return *iter;
            }
        }
        return NULL;
    }

    /**
     * add an edge from tail to head
     */
    void add_edge(node_id tail, node_id head)
    {
        Node * node_tail = get_node(tail);
        Node * node_head = get_node(head);
        Edge * edge = new Edge(node_tail, node_head);
        edge_list_.push_back(edge);
        node_tail->add_out_edge(edge);
    }

};


}

#endif
