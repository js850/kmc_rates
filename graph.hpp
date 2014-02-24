#ifndef _GRAPH_CPP_H_
#define _GRAPH_CPP_H_

#include <cstdlib>
#include <iostream>
#include <list>
#include <map>
#include <set>

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
    typedef std::set<edge_ptr> edge_list;
    typedef typename edge_list::iterator edge_iterator;
private:
    node_id id_;
    edge_list out_edge_list_; // list of outgoing edges
    edge_list in_edge_list_; // list of outgoing edges
public:

    Node(node_id id):
        id_(id)
    {}

    void add_out_edge(edge_ptr edge){ out_edge_list_.insert(edge); }
    void remove_out_edge(edge_ptr edge){ out_edge_list_.erase(edge); }
    edge_list & get_out_edges(){ return out_edge_list_; }
    edge_iterator out_edge_begin(){ return out_edge_list_.begin(); }
    edge_iterator out_edge_end(){ return out_edge_list_.end(); }

    void add_in_edge(edge_ptr edge){ in_edge_list_.insert(edge); }
    void remove_in_edge(edge_ptr edge){ in_edge_list_.erase(edge); }
    edge_list & get_in_edges(){ return in_edge_list_; }
    edge_iterator in_edge_begin(){ return in_edge_list_.begin(); }
    edge_iterator in_edge_end(){ return in_edge_list_.end(); }

    node_id id() const { return id_; }
    size_t out_degree() const { return out_edge_list_.size(); }
    size_t in_degree() const { return in_edge_list_.size(); }
    size_t in_out_degree() const { return out_degree() + in_degree(); }
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



class Graph
{
    std::map<node_id, node_ptr> node_map_;
    std::set<edge_ptr> edge_list_;

    node_id next_node_id_;

public:
    Graph():
        next_node_id_(0)
    {}

    ~Graph()
    {
        // delete all nodes
        typedef std::map<node_id, Node *> maptype;
        for (maptype::iterator iter = node_map_.begin();
                iter != node_map_.end(); ++iter){
            node_ptr node = iter->second;
            delete node;
        }
        // delete all edges
        for (std::set<edge_ptr>::iterator iter = edge_list_.begin();
                iter != edge_list_.end(); ++iter){
            delete *iter;
        }
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
        edge_ptr edge = new Edge(node_tail, node_head);
        edge_list_.insert(edge);
        node_tail->add_out_edge(edge);
        node_head->add_in_edge(edge);
    }

    /**
     * remove a node and all edges connecting it
     */
    void remove_node(node_id nodeid)
    {
        node_ptr u = get_node(nodeid);
        Node::edge_iterator eiter;

        // remove the edges from the nodes connected to u
        for (eiter = u->out_edge_begin(); eiter != u->out_edge_end(); ++eiter){
            edge_ptr edge = *eiter;
            edge->head()->remove_in_edge(edge);
        }
        for (eiter = u->in_edge_begin(); eiter != u->in_edge_end(); ++eiter){
            edge_ptr edge = *eiter;
            edge->tail()->remove_out_edge(edge);
        }

        Node::edge_list to_delete;
        to_delete.insert(u->in_edge_begin(), u->in_edge_end());
        to_delete.insert(u->out_edge_begin(), u->out_edge_end());

        // remove the edges from the edge list
        for (eiter = to_delete.begin(); eiter != to_delete.end(); ++eiter){
            edge_ptr edge = *eiter;
            edge_list_.erase(edge);
        }

        // remove the node from the node list
        node_map_.erase(nodeid);

        // deallocate the memory
        for (eiter = to_delete.begin(); eiter != to_delete.end(); ++eiter){
            edge_ptr edge = *eiter;
            delete edge;
        }

        delete u;
    }

};


}

#endif
