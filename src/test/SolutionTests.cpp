//Link to Boost
#define BOOST_TEST_DYN_LINK

//VERY IMPORTANT - include this last
//#include <boost/test/included/unit_test.hpp>
#include <boost/test/unit_test.hpp>

#include "test.h"
#include "../Solution.h"
using namespace std::literals::complex_literals;

// test suite
BOOST_FIXTURE_TEST_SUITE(BetheSolution_suite, SimpleTestFixture, * utf::label("BetheSolution"))

BOOST_DATA_TEST_CASE(regular_roots_test, bdata::random(4, 16) ^ bdata::xrange(20), L, index)
{
    int M = random64(L/2);
    if (M <= 1) M = 2;
    Vector roots(M);
    for (int i = 0; i < M; i++) {
        roots(i) = randomComplex(0.1, 0.9);
    }
    
    BetheSolution sol(roots, L, 0.0);
    
    BOOST_TEST(sol.size() == M);
    BOOST_TEST(sol.indexOfSingularRoot1() == -1);
    BOOST_TEST(sol.indexOfSingularRoot2() == -1);
    BOOST_TEST(sol.indexOfBlowupRoot() == -1);
    BOOST_TEST(sol.omega() == 0.0il);
    BOOST_TEST(sol.c() == 0.0il);
    for (int i = 0; i < sol.size(); i++) {
        BOOST_TEST(sol.root(i) == roots(i));
        BOOST_TEST(sol.get(i) == roots(i));
        BOOST_TEST(sol.get2(i) == roots(i));
        BOOST_TEST(sol.isSingular(i) == false);
        BOOST_TEST(sol.isLargestRoot(i) == false);
        BOOST_TEST(sol.type(i) == regular);
    }
}

BOOST_DATA_TEST_CASE(singular_roots_test_evenL, bdata::random(2, 8) ^ bdata::xrange(20), halfL, index)
{
    int L = 2 * halfL;
    int M = 2;
    Vector roots(M);
    roots(0) = -0.5il;
    roots(1) = 0.5il;
    BetheSolution sol(roots, L, 0.0);
    BOOST_TEST(sol.size() == M);
    
    BOOST_TEST(chop(sol.c()) == (halfL % 2 == 0 ? 2.0il : -2.0il));
    BOOST_TEST(chop(sol.epsilon()) == 0.0il);
    BOOST_TEST(sol.indexOfSingularRoot1() == 1);
    BOOST_TEST(sol.indexOfSingularRoot2() == 0);

    BOOST_TEST(sol.indexOfBlowupRoot() == -1);
    BOOST_TEST(chop(sol.omega()) == 0.0il);
    
    for (int i = 0; i < sol.size(); i++) {
        BOOST_TEST(chop(sol.get(i)) == roots(i));
        BOOST_TEST(chop(sol.get2(i)) == roots(i));
        BOOST_TEST(sol.isSingular(i));
        BOOST_TEST(sol.isLargestRoot(i) == false);
        BOOST_TEST(sol.type(i) == regular);
    }
}
/*
BOOST_DATA_TEST_CASE(singular_roots_test_oddL, bdata::random(2, 8) ^ bdata::xrange(20), halfL, index)
{
    int L = 2 * halfL + 1;
    int M = 2;
    Vector roots(M);
    roots(0) = -0.5il;
    roots(1) = 0.5il;
    BetheSolution sol(roots, L, 0.0);
    BOOST_TEST(sol.size() == M);
    
    BOOST_TEST(chop(sol.c()) == 0.0L);
    BOOST_TEST(chop(sol.epsilon()) == 0.0L);
    BOOST_TEST(sol.indexOfSingularRoot1() == -1);
    BOOST_TEST(sol.indexOfSingularRoot2() == -1);

    BOOST_TEST(sol.indexOfBlowupRoot() == -1);
    BOOST_TEST(chop(sol.omega()) == 0.0il);
    
    for (int i = 0; i < sol.size(); i++) {
        BOOST_TEST(chop(sol.get(i)) == roots(i));
        BOOST_TEST(chop(sol.get2(i)) == roots(i));
        BOOST_TEST(sol.isSingular(i) == false);
        BOOST_TEST(sol.isLargestRoot(i) == false);
        BOOST_TEST(sol.type(i) == regular);
    }
}*/

BOOST_DATA_TEST_CASE(blowup_roots_test, bdata::random(4, 16) ^ bdata::xrange(20), L, index)
{
    int M = 2 + random64(L/2 - 2);
    int indexOfLargest = random64(M - 1);
    
    Vector roots(M);
    roots(indexOfLargest) = randomComplex(blowup_root_cutoff/0.3, blowup_root_cutoff/0.2);
    
    std::vector<RootType> types(M, regular);
    types[indexOfLargest] = largest;
    for (int i = 0; i < M; i++) {
        if (i == indexOfLargest) continue;
        if (i % 2 == 0) {
            types[i] = large;
            roots(i) = randomComplex(blowup_root_cutoff/0.7, blowup_root_cutoff/0.6);
        } else {
            roots(i) = randomComplex(blowup_root_cutoff * 0.6, blowup_root_cutoff);
        }
    }
    
    BetheSolution sol(roots, L, 0.0);
    BOOST_TEST(sol.indexOfBlowupRoot() == indexOfLargest);
    BOOST_TEST(sol.isLargestRoot(indexOfLargest));
    BOOST_TEST(std::abs(sol.omega() - 1.0L/roots(indexOfLargest)) == 0, tt::tolerance(eps));
    
    BOOST_TEST(chop(sol.c()) == 0.0L);
    BOOST_TEST(chop(sol.epsilon()) == 0.0L);
    BOOST_TEST(sol.indexOfSingularRoot1() == -1);
    BOOST_TEST(sol.indexOfSingularRoot2() == -1);
    
    for (int i = 0; i < sol.size(); i++) {
        BOOST_TEST(sol.type(i) == types[i]);
        BOOST_TEST_INFO("sol.get(i)=" << std::setprecision(15) << sol.get(i) << ", roots(i)=" << roots(i));
        BOOST_TEST(std::abs(sol.get(i) - roots(i)) == 0, tt::tolerance(eps));
        BOOST_TEST(sol.isSingular(i) == false);
    }
}
/*
TODO:
BOOST_DATA_TEST_CASE(update_test, bdata::random(4, 16) ^ bdata::xrange(20), halfL, index)
{
    int L = 2 * halfL;
    int M = halfL;
    Vector roots(M);
    int indexOfSingular1 = 0;
    int indexOfSingular2 = 1;
    int indexOfLargest = 3;
    roots(indexOfSingular1) = -0.5il;
    roots(indexOfSingular2) = 0.5il;
    roots(indexOfLargest) = randomComplex(blowup_root_cutoff/0.3, blowup_root_cutoff/0.2);
    for (int i = 0; i < roots.size(); i++) {
        if (i == indexOfSingular1 || i == indexOfSingular2 || i == indexOfLargest) continue;
        roots(i) = randomComplex(blowup_root_cutoff/0.3, blowup_root_cutoff/0.2)
    }
    BetheSolution sol(roots, L, 0.0);

    int M = 2 + random64(L/2 - 2);
    int indexOfLargest = random64(M - 1);
    
    Vector roots(M);
    roots(indexOfLargest) = randomComplex(blowup_root_cutoff/0.3, blowup_root_cutoff/0.2);
    
    std::vector<RootType> types(M, regular);
    types[indexOfLargest] = largest;
    for (int i = 0; i < M; i++) {
        if (i == indexOfLargest) continue;
        if (i % 2 == 0) {
            types[i] = large;
            roots(i) = randomComplex(blowup_root_cutoff/0.7, blowup_root_cutoff/0.6);
        } else {
            roots(i) = randomComplex(blowup_root_cutoff * 0.6, blowup_root_cutoff);
        }
    }
    
    BetheSolution sol(roots, L, 0.0);
    BOOST_TEST(sol.indexOfBlowupRoot() == indexOfLargest);
    BOOST_TEST(sol.isLargestRoot(indexOfLargest));
    BOOST_TEST(std::abs(sol.omega() - 1.0L/roots(indexOfLargest)) == 0, tt::tolerance(eps));
    
    BOOST_TEST(chop(sol.c()) == 0.0L);
    BOOST_TEST(chop(sol.epsilon()) == 0.0L);
    BOOST_TEST(sol.indexOfSingularRoot1() == -1);
    BOOST_TEST(sol.indexOfSingularRoot2() == -1);
    
    for (int i = 0; i < sol.size(); i++) {
        BOOST_TEST(sol.type(i) == types[i]);
        BOOST_TEST_INFO("sol.get(i)=" << std::setprecision(15) << sol.get(i) << ", roots(i)=" << roots(i));
        BOOST_TEST(std::abs(sol.get(i) - roots(i)) == 0, tt::tolerance(eps));
        BOOST_TEST(sol.isSingular(i) == false);
    }
}*/


BOOST_AUTO_TEST_SUITE_END()

