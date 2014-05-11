#include <gtest/gtest.h>
//#include <pair>
#include "kmc_rates/graph.hpp"
#include "kmc_rates/ngt.hpp"

using graph_ns::NGT;
using graph_ns::node_id;

class NGT3 :  public ::testing::Test
{
public:
	NGT::rate_map_t rate_map;
	std::set<node_id> A, B;
    virtual void SetUp(){
        for (int i=0; i<3; ++i){
            for (int j=0; j<3; ++j){
                rate_map[std::pair<node_id, node_id>(i,j)] = 1;
            }
        }
        A.insert(0);
        B.insert(1);


    }
};

TEST_F(NGT3, Rate_Correct){
	NGT ngt(rate_map, A, B);
	ngt.compute_rates();
	ASSERT_NEAR(ngt.get_rate_AB(), 1, 1e-9);
	ASSERT_NEAR(ngt.get_rate_BA(), 1, 1e-9);
	ASSERT_NEAR(ngt.get_rate_AB_SS(), 1.5, 1e-9);
	ASSERT_NEAR(ngt.get_rate_BA_SS(), 1.5, 1e-9);
}

TEST_F(NGT3, RateCommittors_Correct){
	NGT ngt(rate_map, A, B);
	ngt.compute_rates_and_committors();
	ASSERT_NEAR(ngt.get_rate_AB(), 1, 1e-9);
	ASSERT_NEAR(ngt.get_rate_BA(), 1, 1e-9);
	ASSERT_NEAR(ngt.get_rate_AB_SS(), 1.5, 1e-9);
	ASSERT_NEAR(ngt.get_rate_BA_SS(), 1.5, 1e-9);

	auto committors = ngt.get_committors();
	ASSERT_NEAR(committors[2], 0.5, 1e-9);

}

