all: test test_ngt


test: test.cpp graph.hpp breadth_first_search.hpp grid_graph.hpp
	g++ -o test -g -Wall $<

test_ngt: test_ngt.cpp graph.hpp breadth_first_search.hpp grid_graph.hpp ngt.hpp
	g++ -o test_ngt -g -Wall $<
