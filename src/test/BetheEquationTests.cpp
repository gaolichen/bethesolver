//Link to Boost
#define BOOST_TEST_DYN_LINK

//VERY IMPORTANT - include this last
//#include <boost/test/included/unit_test.hpp>
#include <boost/test/unit_test.hpp>

#include "test.h"

// test suite
BOOST_FIXTURE_TEST_SUITE(BetheEquation_suite, SimpleTestFixture, * utf::label("BetheEquation"))
/*
BOOST_DATA_TEST_CASE(eval_test, bdata::random(4, 10) ^ bdata::xrange(20), L, index)
{
    int M = random64(L/2);
    if (M <= 1) M = 2;
    SolvedBetheEquation be(L, M);
    BetheRootCache cache(L, M);
    for (int i = 0; i < cache.numberOfRoots(); i++) {
        if (cache.isSingular(i)) continue;
        Vector root = cache.getRoot(i);
        Vector res = be.eval(root);
        for (int j = 0; j < res.size(); j++) {
            BOOST_TEST_INFO("M=" << M << ", root=" << root);
            BOOST_TEST(std::abs(res[j] - var_t(1.0, 0)) == 0, tt::tolerance(1e-5));
        }
    }
}
*/
/*
BOOST_DATA_TEST_CASE(diff_test, bdata::random(4, 10) ^ bdata::xrange(20), L, index)
{
    int M = random64(L/2);
    if (M <= 1) M = 2;
    SolvedBetheEquation be(L, M);
    BetheRootCache cache(L, M);
    elem_t r = 1.0e-4;
    elem_t phi = random(0, 2 * Pi);
    var_t epsilon(r * cos(phi), r * sin(phi));
    Vector varies(M);
    for (int i = 0; i < M; i++) {
        varies[i] = epsilon;
    }
    for (int i = 0; i < cache.numberOfRoots(); i++) {
        Vector root = cache.getRoot(i);
        Vector oldRes = be.eval(root);
        root += varies;
        Vector newRes = be.eval(root);
        Vector expected = be.diff(root) * varies;
        BOOST_TEST_INFO("oldRes=" << oldRes << ", newRes=" << newRes << ", expected dif=" << expected);
        //std::cout << "oldRes=" << oldRes << ", newRes=" << newRes << ", expected dif=" << expected << std::endl;
        BOOST_TEST((newRes - oldRes -expected).norm() < 1e-3);
    }
}*/

/*
BOOST_DATA_TEST_CASE(impliciteDiff_test, bdata::random(4, 10) ^ bdata::xrange(5), L, index)
{
    int M = random64(L/2);
    if (M <= 1) M = 2;
    SolvedBetheEquation be(L, M);
    SimpleBetheHomotopy bh(0.9, &be);
    
    BetheRootCache cache(L, M);
    elem_t deltaBeta = 1.0e-3;
    for (int i = 0; i < cache.numberOfRoots(); i++) {
        Vector root = cache.getRoot(i);
        Vector uDiff = bh.impliciteDiff(0, root);
        Vector res1 = be.eval(root);
        Vector res2 = be.eval(root + uDiff * deltaBeta);
        Vector diff(M);
        for (int j = 0; j < res1.size(); j++) {
            diff[j] = log(res2[j]/res1[j])/var_t(0, (2.0 * L) * deltaBeta) - 1.0;
        }
        BOOST_TEST_INFO("M=" << M << ", root=" << root <<", diff=" << uDiff * deltaBeta << ", res1=" << res1 << ", res2=" << res2);
        BOOST_TEST(diff.norm() < 1e-1);
    }
}*/

BOOST_AUTO_TEST_SUITE_END()
