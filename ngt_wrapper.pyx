from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp.pair cimport pair
from libcpp.list cimport list as stdlist


cdef extern from "graph.hpp" namespace "graph_ns":
    ctypedef unsigned long node_id

ctypedef pair[node_id, node_id] pair_t
ctypedef map[pair_t, double] rate_map_t

cdef extern from "ngt.hpp" namespace "graph_ns":
    cdef cppclass cNGT "graph_ns::NGT":
        cNGT(rate_map_t &, stdlist[node_id] &, stdlist[node_id] &) except +
        void compute()
        double get_rate_AB()
        double get_rate_BA()



cdef make(rate_constants):
    cdef map[node_id, double] Peq
    Peq[1] = 0.
    
    cdef pair[node_id, node_id] p
    cdef map[pair_t, double] rates
    rates[pair_t(1, 2)] = 1.

cdef class BaseNGT(object):
    cdef cNGT* thisptr
    cdef node_list
    cdef node2id
    
    def __cinit__(self, rates, A, B):
        # assign ids to all the nodes
        nodes = set()
        for u, v in rates.iterkeys():
            nodes.add(u)
            nodes.add(v)
        self.node_list = list(nodes)
        self.node2id = dict(( (u, i) for i, u in enumerate(self.node_list) ))
        
        cdef rate_map_t rate_map
        cdef node_id uid, vid
        for (u, v), k in rates.iteritems():
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
    
    def compute(self):
        self.thisptr.compute()
    
    def get_rate_AB(self):
        return self.thisptr.get_rate_AB()
    def get_rate_BA(self):
        return self.thisptr.get_rate_BA()
        
    

class NGT(BaseNGT):
    pass
    
