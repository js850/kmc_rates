"""
# distutils: language = C++
"""
import time

from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp.pair cimport pair
from libcpp.list cimport list as stdlist


cdef extern from "source/graph.hpp" namespace "graph_ns":
    ctypedef unsigned long node_id

ctypedef pair[node_id, node_id] pair_t
ctypedef map[pair_t, double] rate_map_t

cdef extern from "source/ngt.hpp" namespace "graph_ns":
    cdef cppclass cNGT "graph_ns::NGT":
        cNGT(rate_map_t &, stdlist[node_id] &, stdlist[node_id] &) except +
        void compute() except +
        double get_rate_AB() except +
        double get_rate_BA() except +
        void set_node_occupation_probabilities(map[node_id, double] &) except +
        void set_debug() except +




cdef class BaseNGT(object):
    """
    class to apply the graph reduction method for finding transition rates between two groups of nodes
    
    Parameters
    ----------
    rate_constants : dict
        a dictionary of rates.  the keys are tuples of nodes (u,v), the values
        are the rates.
        
            rate_uv = rate[(u,v)]
        
    A, B : iterables
        Groups of nodes specifying the reactant and product groups.  The rates
        returned will be the rate from A to B and vice versa.
    weights : dict
        Dictionary with nodes as keys and weights as values.  The weights are
        the equilibrium occupation probabilities of the nodes in A and B.  They
        are used to do the weighted mean for the final average over inverse
        mean first passage times.
    
    Notes
    -----
    This follows the new graph transformation procedure (NGT) described by 
    David Wales, J. Chem. Phys., 2009 http://dx.doi.org/10.1063/1.3133782
    
    The rate, rAB computed by this calculation (returned by
    self.get_final_rates) is the inverse mean first passage time averaged over
    the states in A
    
    """
    cdef cNGT* thisptr
    cdef node_list
    cdef node2id
    time_solve = 0.
    def __cinit__(self, rate_constants, A, B, debug=False, weights=None):
        # assign ids to all the nodes
        nodes = set()
        for u, v in rate_constants.iterkeys():
            nodes.add(u)
            nodes.add(v)
        self.node_list = list(nodes)
        self.node2id = dict(( (u, i) for i, u in enumerate(self.node_list) ))
        
        cdef rate_map_t rate_map
        cdef node_id uid, vid
        for (u, v), k in rate_constants.iteritems():
            uid = self.node2id[u]
            vid = self.node2id[v]
            rate_map[pair_t(uid, vid)] = k
        
        cdef stdlist[node_id] _A, _B
        for u in A:
            uid = self.node2id[u]
            _A.push_back(uid)
        for u in B:
            uid = self.node2id[u]
            _B.push_back(uid)
            
        
        self.thisptr = new cNGT(rate_map, _A, _B)
        cdef map[node_id, double] Peq 
        if weights is not None:
            for u, p in weights.iteritems():
                try:
                    uid = self.node2id[u]
                    Peq[uid] = p
                except KeyError:
                    pass
                        
            self.thisptr.set_node_occupation_probabilities(Peq)
        
        if debug:
            self.thisptr.set_debug()
    
    def __dealloc__(self):
        if self.thisptr != NULL:
            del self.thisptr
            self.thisptr = NULL
    
    def compute(self):
        t0 = time.clock()
        self.thisptr.compute()
        self.time_solve = time.clock() - t0 
    
    def get_rate_AB(self):
        return self.thisptr.get_rate_AB()
    def get_rate_BA(self):
        return self.thisptr.get_rate_BA()
        
    

class NGT(BaseNGT):
    pass
    
