#ifndef _GRAPH_CPP_H_
#define _GRAPH_CPP_H_

#include <cstdlib>
#include <iostream>
#include <list>

namespace graph
{
class Edge;
class Node;

class Node{
public:
    size_t id_;
    Edge * first_outgoing_; // first outgoing edge

    Node(size_t id):
        id_(id),
        first_outgoing_(NULL)
    {}

    void add_out_edge(Edge * edge);
};

class Edge{
public:
    Node * head_; // node the edge points to
    Node * tail_; // node the edge comes from
    Edge * next_; // next edge with same originating node

    Edge(Node * head, Node * tail):
        head_(head),
        tail_(tail),
        next_(NULL)
    {}
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

    size_t next_node_id_;

public:
    Graph():
        next_node_id_(0)
    {}

    ~Graph()
    {
        for (std::list<Node *>::iterator iter = node_list_.begin();
                iter != node_list_.end(); ++iter){
            delete *iter;
        }
        for (std::list<Edge *>::iterator iter = edge_list_.begin();
                iter != edge_list_.end(); ++iter){
            delete *iter;
        }
    }

    /**
     * create a new node
     */
    size_t add_node(){
        Node *node = new Node(next_node_id_);
        next_node_id_++;
        node_list_.push_back(node);
        return node->id_;
    }

    /**
     * return the node with given node id
     */
    Node * get_node(size_t node_id)
    {
        for (std::list<Node *>::iterator iter = node_list_.begin();
                iter != node_list_.end(); ++iter){
            if ((*iter)->id_ == node_id){
                return *iter;
            }
        }
    }

    /**
     * add an edge from tail to head
     */
    void add_edge(size_t tail, size_t head)
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
