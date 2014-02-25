#include "graph.hpp"
#include "breadth_first_search.hpp"
#include "grid_graph.hpp"
#include "ngt.hpp"
#include <iostream>
#include <fstream>
#include <sstream>

using std::cout;
using namespace graph_ns;
using graph_ns::Graph;
using graph_ns::node_id;
using graph_ns::edge_ptr;
using graph_ns::node_ptr;

NGT::rate_map_t make_rates(int nnodes){
    typedef std::pair<node_id, node_id> pair_t;
    NGT::rate_map_t rates;

	for (int i=0; i<nnodes; ++i){
		for (int j=i+1; j<nnodes; ++j){
			rates[pair_t(i,j)] = float(i+j) / (i+1);
			rates[pair_t(j,i)] = float(i+j) / (j+1);
		}
	}
	return rates;
}

template<class T>
size_t read1(const std::string file_name, T& table) {
    table.clear();
    size_t lines_read = 0;
    std::ifstream in(file_name.c_str());
    std::string line;
    while (std::getline(in, line)) {
        ++lines_read;
        std::istringstream iss(line);
        node_id u;
        iss >> u;
//        cout << u << "\n";
        table.push_back(u); // so that we don't loose trailing tabs
    }
    return lines_read;
}
template<class T>
size_t read2(const std::string file_name, T& table) {
    table.clear();
    size_t lines_read = 0;
    std::ifstream in(file_name.c_str());
    std::string line;
    while (std::getline(in, line)) {
        ++lines_read;
        std::istringstream iss(line);
        node_id u;
        double p;
        iss >> u;
        iss >> p;
//        cout << u << " " << p << "\n";
        table[u] = p; // so that we don't loose trailing tabs
    }
    return lines_read;
}
template<class T>
size_t read3(const std::string file_name, T& table) {
    table.clear();
    size_t lines_read = 0;
    std::ifstream in(file_name.c_str());
    std::string line;
    while (std::getline(in, line)) {
        ++lines_read;
        std::istringstream iss(line);
        node_id u, v;
        double k;
        iss >> u;
        iss >> v;
        iss >> k;
//        cout << u << " " << v << " " << k << "\n";
        table[std::pair<node_id, node_id>(u,v)] = k; // so that we don't loose trailing tabs
    }
    return lines_read;
}

void run_from_optim(){
    std::list<node_id> A, B;
    std::map<node_id, double> Peq;
    NGT::rate_map_t rates;
    read1("in.A", A);
    read1("in.B", B);
    read2("in.Peq", Peq);
    read3("in.rates", rates);

    NGT ngt(rates, A, B);
    ngt.set_node_occupation_probabilities(Peq);
    ngt.compute();
    cout << "rate A -> B " << ngt.get_rate_AB() << "\n";
    cout << "rate B -> A " << ngt.get_rate_BA() << "\n";

}

void run(){
    std::vector<node_id> A, B;
    A.push_back(0);
    A.push_back(1);
    A.push_back(2);
    B.push_back(3);
    B.push_back(4);

    std::map<node_id, double> Peq;
    Peq[0] = 1.;
    Peq[1] = 1.;
    Peq[2] = 5e3;
    Peq[3] = 1.;
    Peq[4] = 1.;
    Peq[5] = 1.;

    NGT::rate_map_t rates = make_rates(100);
    NGT ngt(rates, A, B);
    ngt.set_node_occupation_probabilities(Peq);
//    ngt.debug = true;

    ngt.phase_one();
    ngt.phase_two();
    cout << "rate A -> B " << ngt.get_rate_AB() << "\n";
    cout << "rate B -> A " << ngt.get_rate_BA() << "\n";

}


void run3()
{
    std::vector<node_id> A, B;
    A.push_back(0);
    B.push_back(1);

    typedef std::pair<node_id, node_id> pair_t;
    NGT::rate_map_t rates;
    rates[pair_t(0,1)] = 1.;
    rates[pair_t(1,0)] = 1.;
    rates[pair_t(0,2)] = 1.;
    rates[pair_t(2,0)] = 1.;
    rates[pair_t(1,2)] = 1.;
    rates[pair_t(2,1)] = 1.;

    NGT ngt(rates, A, B);
//    ngt.debug = true;
    node_ptr a = ngt._graph->get_node(0);
    node_ptr b = ngt._graph->get_node(1);
    edge_ptr ab = a->get_successor_edge(b);

    double tau_a = ngt.get_tau(a);
    double Paa = ngt.get_node_P(a);
    double Pab = ngt.get_P(ab);
    cout << "a tau " << tau_a << "\n";
    cout << "a P " << Paa << "\n";
    cout << "ab P " << Pab << "\n";

    ngt.phase_one();
    ngt.phase_two();


    tau_a = ngt.get_tau(a);
    Paa = ngt.get_node_P(a);
    Pab = ngt.get_P(ab);
    cout << "tau_a " << tau_a << "\n";
    cout << "Paa " << Paa << "\n";
    cout << "Pab " << Pab << "\n";
    cout << "rate A -> B " << ngt.get_rate_AB() << "\n";
    cout << "rate B -> A " << ngt.get_rate_BA() << "\n";

}
int main(){
    run_from_optim();
//    run();
}
