#ifndef _NGT_HPP_
#define _NGT_HPP_

#include <cstdlib>
#include <iostream>
#include <list>
#include <queue>
#include <assert.h>
#include <stdexcept>

#include "graph.hpp"

using std::cout;

namespace graph_ns
{

bool compare_degree(node_ptr u, node_ptr v){
    return u->in_out_degree() < v->in_out_degree();
}

class NGT {
public:
    typedef std::map<std::pair<node_id, node_id>, double> rate_map_t;

    Graph * _graph;
    std::set<node_ptr> _A, _B;
    std::list<node_ptr> intermediates; //this will an up to date list of nodes keyed by the node degree
    bool debug;
    bool own_graph;

    std::map<node_id, double> final_omPxx;
    std::map<node_id, double> final_tau;

    
    ~NGT()
    {
        if (own_graph && _graph != NULL){
            delete _graph;
            _graph = NULL;
        }
    }

    template<class Acontainer, class Bcontainer>
    NGT(Graph & graph, Acontainer &A, Bcontainer &B) :
        _graph(& graph),
        debug(false),
        own_graph(false)
    {
    	for (typename Acontainer::iterator iter = A.begin(); iter != A.end(); ++iter){
    		_A.insert(_graph->get_node(*iter));
    	}
    	for (typename Bcontainer::iterator iter = B.begin(); iter != B.end(); ++iter){
    		_B.insert(_graph->get_node(*iter));
    	}

    	// make intermediates
    	for (Graph::node_map_t::iterator miter = _graph->node_map_.begin(); miter != _graph->node_map_.end(); ++miter){
    		node_ptr u = miter->second;
    		if (_A.find(u) == _A.end() and _B.find(u) == _B.end()){
    			intermediates.push_back(u);
    		}
    	}

		cout << "number of nodes " << _graph->number_of_nodes() << "\n";
		cout << "A.size() " << _A.size() << "\n";
		cout << "B.size() " << _B.size() << "\n";
		cout << "intermediates.size() " << intermediates.size() << "\n";
		assert(intermediates.size() + _A.size() + _B.size() == _graph->number_of_nodes());

    }

    NGT(rate_map_t &rate_constants, std::vector<node_id> &A, std::vector<node_id> &B) :
        _graph(new Graph()),
        debug(false),
        own_graph(true)
    {
        std::set<node_ptr> nodes;

        // add nodes to the graph and sum the rate constants for all out edges for each node.
        std::map<node_ptr, double> sum_out_rates;
        typedef std::map<std::pair<node_id, node_id>, double> maptype;
        for (maptype::iterator iter = rate_constants.begin(); iter != rate_constants.end(); ++iter){
            node_ptr u = _graph->add_node(iter->first.first);
            node_ptr v = _graph->add_node(iter->first.second);
            double k = iter->second;
            nodes.insert(u);
            nodes.insert(v);

            try {
                sum_out_rates.at(u) += k;
            } catch (std::out_of_range & e) {
                sum_out_rates[u] = k;
            }
        }

        // set tau_x for each node
        // add edge Pxx for each node and initialize P to 0.
        for (std::set<node_ptr>::iterator uiter = nodes.begin(); uiter != nodes.end(); ++uiter){
            node_ptr x = *uiter;
            double tau_x = 1. / sum_out_rates[x];
            set_tau(x, tau_x);
            edge_ptr xx = _graph->_add_edge(x, x);
            set_P(xx, 0.);
        }

        // set Puv for each edge
        typedef std::map<std::pair<node_id, node_id>, double> maptype;
        for (maptype::iterator iter = rate_constants.begin(); iter != rate_constants.end(); ++iter){
            node_ptr u = _graph->get_node(iter->first.first);
            node_ptr v = _graph->get_node(iter->first.second);
            double k = iter->second;

            edge_ptr uv = _graph->_add_edge(u, v);
            double tau_u = get_tau(u);
            double Puv = k * tau_u;
            set_P(uv, Puv);

            try {
                sum_out_rates.at(u) += k;
            } catch (std::out_of_range & e) {
                sum_out_rates[u] = k;
            }
        }


        // make the set of A and B
        for (std::vector<node_id>::iterator node_iter = A.begin(); node_iter != A.end(); ++node_iter){
            _A.insert(_graph->get_node(*node_iter));
        }
        for (std::vector<node_id>::iterator node_iter = B.begin(); node_iter != B.end(); ++node_iter){
            _B.insert(_graph->get_node(*node_iter));
        }

        // make a list of intermediates
        for (std::set<node_ptr>::iterator uiter = _A.begin(); uiter != _A.end(); ++uiter){
            nodes.erase(*uiter);
        }
        for (std::set<node_ptr>::iterator uiter = _B.begin(); uiter != _B.end(); ++uiter){
            nodes.erase(*uiter);
        }
        intermediates.assign(nodes.begin(), nodes.end());


//        cout << _graph->number_of_nodes() << "\n";
//        cout << _A.size() << "\n";
//        cout << _B.size() << "\n";
//        cout << intermediates.size() << "\n";
//        cout << nodes.size() << "\n";
        assert(intermediates.size() + _A.size() + _B.size() == _graph->number_of_nodes());


    }

    void sort_intermediates(){
        node_ptr x = *intermediates.begin();
        if (debug){
            cout << "smallest node degree " << x->in_out_degree() << "\n";
        }
        if (x->in_out_degree() > 4) {
            intermediates.sort(compare_degree);
        }
    }
    
    inline double get_tau(node_ptr u){ return u->tau; }
    inline double get_P(edge_ptr edge){ return edge->P; }
    edge_ptr get_node_self_edge(node_ptr u){ 
        return u->get_successor_edge(u);
    }
    double get_node_P(node_ptr u){ return get_P(u->get_successor_edge(u)); }
    double get_node_one_minus_P(node_ptr u){ return 1. - get_P(u->get_successor_edge(u)); }

    inline void set_tau(node_ptr u, double tau){ u->tau = tau; }
    inline void set_P(edge_ptr edge, double P){ edge->P = P; }

    /*
     * node x is being deleted, so update P and tau for node u
     */
    void update_tau(edge_ptr ux, double omPxx, double tau_x){
        node_ptr u = ux->tail();
        double Pux = get_P(ux);
        double tau_u = get_tau(u);
        double new_tau_u = tau_u + Pux * tau_x / omPxx;
        if (debug){
            cout << "updating node " << u->id() << " tau " << tau_u << " -> " << new_tau_u << "\n";
        }
        set_tau(u, new_tau_u);
    }

    edge_ptr add_edge(node_ptr u, node_ptr v){
       edge_ptr edge = _graph->_add_edge(u, v);
       set_P(edge, 0.);
       return edge;
    }

    void update_edge(node_ptr u, node_ptr v, node_ptr x, edge_ptr ux, double omPxx){
        edge_ptr xv = x->get_successor_edge(v);
        if (ux == NULL || xv == NULL){
            // no need to do anything if either of these don't exist
            return;
        }
        edge_ptr uv = u->get_successor_edge(v);
        if (uv == NULL){
            uv = add_edge(u, v);
        }

        double Pux = get_P(ux);
        double Pxv = get_P(xv);
        double Puv = get_P(uv);

        double newPuv = Puv + Pux * Pxv / omPxx;
        if (debug) {
            cout << "updating edge " << u->id() << " -> " << v->id() << " Puv " << Puv << " -> " << newPuv
                    << " 1-Pxx " << omPxx
                    << " Pux " << Pux
                    << " Pxv " << Pxv
                    << "\n";
        }
        set_P(uv, newPuv);
    }

    void remove_node(node_ptr x){
        double taux = get_tau(x);
//        double Pxx = get_node_P(x);
        double omPxx = get_node_one_minus_P(x);

        // update the node data for all the neighbors
        Node::edge_iterator eiter;
        for (eiter = x->in_edge_begin(); eiter != x->in_edge_end(); eiter++){
            edge_ptr edge = *eiter;
            if (edge->tail() != edge->head()){
                update_tau(edge, omPxx, taux);
            }
        }

        std::set<node_ptr> neibs = x->in_out_neighbors();
        neibs.erase(x);

        //
        for (std::set<node_ptr>::iterator uiter = neibs.begin(); uiter != neibs.end(); ++uiter){
            node_ptr u = *uiter;
            edge_ptr ux = u->get_successor_edge(x);
            if (ux == NULL) continue;
            for (std::set<node_ptr>::iterator viter = neibs.begin(); viter != neibs.end(); ++viter){
                node_ptr v = *viter;
                update_edge(u, v, x, ux, omPxx);
            }
        }

        // remove the node from the graph
        _graph->_remove_node(x);

    }

    void remove_intermediates(){
        while (intermediates.size() > 0){
            sort_intermediates();

            node_ptr x = intermediates.front();
            intermediates.pop_front();
            std::cout << "removing node " << x->id() << "\n";

            remove_node(x);
        }
    }

    void phase_one(){
    	remove_intermediates();
    }

    void reduce_all_in_group(std::set<node_ptr> &to_remove, std::set<node_ptr> & to_keep){
    	std::set<node_id> Aids, Bids;
    	for (std::set<node_ptr>::iterator iter = to_remove.begin(); iter != to_remove.end(); ++iter){
    		Aids.insert((*iter)->id());
    	}
    	for (std::set<node_ptr>::iterator iter = to_keep.begin(); iter != to_keep.end(); ++iter){
    		Bids.insert((*iter)->id());
    	}
    	for (std::set<node_id>::iterator iter = Aids.begin(); iter != Aids.end(); ++iter){
    		node_id x = *iter;
    		std::set<node_id> newAids;
    		newAids.insert(x);

    		Graph new_graph(*_graph);
    		NGT new_ngt(new_graph, newAids, Bids);
    		new_ngt.remove_intermediates();
    		node_ptr xptr = new_graph.get_node(x);
    		final_omPxx[x] = new_ngt.get_node_one_minus_P(xptr);
    		final_tau[x] = new_ngt.get_tau(xptr);
    	}
    }

    void phase_two(){
    	reduce_all_in_group(_A, _B);
    	reduce_all_in_group(_B, _A);
    }

    double get_rate_AB(){
        assert(_A.size() == 1);
        node_ptr a = *_A.begin();
        double omPxx = final_omPxx.at(a->id());
        double tau_a = final_tau.at(a->id());
        return omPxx / tau_a;
    }

    double get_rate_BA(){
        assert(_B.size() == 1);
        node_ptr b = *_B.begin();
        double PbA = get_node_one_minus_P(b);
        double tau_b = get_tau(b);
        return PbA / tau_b;
    }


};

}
#endif
