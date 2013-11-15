#ifndef _GRID_GRAPH_HPP_
#define _GRID_GRAPH_HPP_

#include <cstdlib>
#include <iostream>
#include <list>
#include <queue>
#include <assert.h>

#include "graph.hpp"

namespace graph_ns
{
class BuildGridGraph2d{

    Graph & graph_;
    size_t Lx;
    size_t Ly;
    bool periodic_;

public:
    BuildGridGraph2d(Graph & graph, size_t n, size_t m, bool periodic):
        graph_(graph), Lx(n), Ly(m), periodic_(periodic)
    {
        assert(graph.number_of_nodes() == 0);
    }

    size_t xy2i(size_t x, size_t y){
        return x*Lx + y;
    }

    void build_graph()
    {
        graph_.add_nodes(Lx*Ly);

        // add edge to site to the right
        for (size_t x=0; x<Lx-1; ++x){
            for (size_t y=0; y<Ly; ++y){
                graph_.add_edge(xy2i(x,y), xy2i(x+1,y));
            }
        }

        // add edge to site above
        for (size_t x=0; x<Lx; ++x){
            for (size_t y=0; y<Ly-1; ++y){
                graph_.add_edge(xy2i(x,y), xy2i(x,y+1));
            }
        }

        if (periodic_){
            add_periodic_edges();
        }
    }

    void add_periodic_edges()
    {
        // connect top and bottom
        for (size_t x=0; x<Lx; ++x){
            graph_.add_edge(xy2i(x,0), xy2i(x,Ly-1));
        }

        // connect left and right
        for (size_t y=0; y<Ly; ++y){
            graph_.add_edge(xy2i(0,y), xy2i(Lx-1,y));
        }
    }
};

}
#endif
