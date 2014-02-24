#ifndef _NGT_HPP_
#define _NGT_HPP_

#include <cstdlib>
#include <iostream>
#include <list>
#include <queue>
#include <assert.h>
#include <stdexcept>

#include "graph.hpp"

namespace graph_ns
{

bool compare_degree(node_ptr u, node_ptr v){
    return u->in_out_degree() < v->in_out_degree();
}

class NGT {
    Graph & _graph;
    std::set<node_ptr> _A, _B;
    std::list<node_ptr> intermediates; //this will an up to date list of nodes keyed by the node degree

    /*
    NGT(Graph & graph, std::vector<node_id> A, std::vector<node_id> B) :
        _graph(graph)
    {
        for (std::vector<node_id>::iterator node_iter = A.begin(); node_iter != A.end(); ++node_iter){
            _A.insert(graph.get_node(*node_iter));
        }
        for (std::vector<node_id>::iterator node_iter = B.begin(); node_iter != B.end(); ++node_iter){
            _B.insert(graph.get_node(*node_iter));
        }

        for (Graph::node_map_t::iterator miter = graph.node_map_.begin(); miter != graph.node_map_.end(); ++miter){
            node_ptr u = miter->second;
        }
    }
    */

    void sort_intermediates(){
        node_ptr x = *intermediates.begin();
        if (x->in_out_degree() > 1) {
            intermediates.sort(compare_degree);
        }
    }
    
    inline double get_node_tau(node_ptr u){ return u->tau; }
    inline double get_edge_P(edge_ptr edge){ return edge->P; }
    edge_ptr get_node_self_edge(node_ptr u){ 
        Node::edge_iterator eiter;
        for (eiter = u->out_edge_begin(); eiter != u->out_edge_end(); eiter++){
            node_ptr v = (*eiter)->head();
            if (u == v) return *eiter;
        }
        throw std::runtime_error("no edge connecting itself"); 
    }
    double get_node_P(node_ptr u){ return get_edge_P(get_node_self_edge(u)); }
    double get_node_one_minus_P(node_ptr u){ return 1. - get_edge_P(get_node_self_edge(u)); }

    inline void set_node_tau(node_ptr u, double tau){ u->tau = tau; }
    inline void set_edge_P(edge_ptr edge, double P){ edge->P = P; }

    /*
     * node x is being deleted, so update P and tau for node u
     */
    void update_tau(edge_ptr ux, double omPxx, double taux){
        node_ptr u = ux->tail();
        double Pux = get_edge_P(ux);
        double tau_u = get_node_tau(u);
        double new_tau_u = tau_u + Pux * taux / omPxx; 
        set_node_tau(u, new_tau_u);
    }

    //void update_edge(edge_ptr uv, ){
        //double newPux = Pux * tau_x / omPxx;
    //}

    void remove_node(node_ptr x){
        double taux = get_node_tau(x);
        double Pxx = get_node_P(x);
        double omPxx = get_node_one_minus_P(x);

        // update the node data for all the neighbors
        Node::edge_iterator eiter;
        for (eiter = x->in_edge_begin(); eiter != x->in_edge_end(); eiter++){
            edge_ptr edge = *eiter;
            update_tau(edge, omPxx, taux);
        }

        std::set<node_ptr> neibs = x->in_out_neighbors();
        neibs.erase(x);
        std::set<std::pair<node_ptr, node_ptr> > edges;
        // make a set of all the edges that should be present
        for (std::set<node_ptr>::iterator uiter = neibs.begin(); uiter != neibs.end(); ++uiter){
            node_ptr u = *uiter;
            for (std::set<node_ptr>::iterator viter = neibs.begin(); viter != neibs.end(); ++viter){
                node_ptr v = *viter;
                edges.insert(std::pair<node_ptr, node_ptr>(u, v));
            }
        }

        /*
        // now loop through all the edges that are present
        for (std::set<node_ptr>::iterator uiter = neibs.begin(); uiter != neibs.end(); ++uiter){
            node_ptr u = *uiter;
            for (eiter = u->in_edge_begin(); eiter != u->in_edge_end(); eiter++){
                edge_ptr uv = *eiter;
                node_ptr v = uv->head();
                edges.erase(std::pair<node_ptr, node_ptr>(u, v));
                update_edge(uv);
            }
        }

        // for each edge remaining in edges, first add the edge then update it
        */

    }

    void remove_intermediates(){
        while (intermediates.size() > 0){
            sort_intermediates();

            node_ptr x = *intermediates.begin();

            remove_node(x);
        }
    }


};

}
#endif
