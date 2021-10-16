//Link to Boost
#define BOOST_TEST_DYN_LINK

//VERY IMPORTANT - include this last
//#include <boost/test/included/unit_test.hpp>
#include <boost/test/unit_test.hpp>

#include "test.h"
#include "../BetheRootCache.h"

// test suite
BOOST_FIXTURE_TEST_SUITE(BetheRootCache_suite, SimpleTestFixture, * utf::label("BetheRootCache"))

BOOST_AUTO_TEST_CASE(Demo1)
{
    BetheRootCache cache(6, 2);
    for (int i = 0; i < cache.numberOfRoots(); i++) {
        std::cout << i << ": " << cache.getRoot(i) << ", momentum=" << cache.momentum(i) << std::endl;
    }
}

BOOST_AUTO_TEST_SUITE_END()
