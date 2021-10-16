//Link to Boost
#define BOOST_TEST_DYN_LINK

//VERY IMPORTANT - include this last
//#include <boost/test/included/unit_test.hpp>
#include <boost/test/unit_test.hpp>

#include "test.h"
#include "../InfiniteRoots.h"
using namespace std::literals::complex_literals;

// test suite
BOOST_FIXTURE_TEST_SUITE(InfiniteRoots_suite, SimpleTestFixture, * utf::label("InfiniteRoots"))

BOOST_DATA_TEST_CASE(twoRoots_test, bdata::random(0, 8) ^ bdata::random(2, 8) ^ bdata::xrange(20), diffLM, M, index)
{
    int L =  2 * M + diffLM;
    int m = 2;
    InfiniteRoots ir(L, M);
    std::vector<var_t>& roots = ir.getRoot(m);
    var_t r1 = 2.0L * L / var_t(3.0L + L - 2.0L * M, sqrt(3.0L + L - 2.0L * M));
    var_t r2 = 2.0L * L / var_t(3.0L + L - 2.0L * M, -sqrt(3.0L + L - 2.0L * M));
    BOOST_TEST(std::abs(roots[0] - std::conj(roots[1])) < EPS);
    BOOST_TEST(std::abs(roots[0].real() - r1.real()) < EPS);
    BOOST_TEST(abs(abs(roots[0].imag()) - abs(r1.imag())) < EPS);
}

BOOST_DATA_TEST_CASE(moreThanTwoRoots_test, bdata::random(0, 8) ^ bdata::random(3, 16) ^ bdata::xrange(20), diffLM, M, index)
{
    int L = 2 * M + diffLM;
    int m = random64(M - 2) + 2;
    InfiniteRoots ir(L, M);
    std::vector<var_t>& roots = ir.getRoot(m);
    int R = L - 2 * M + 2 * m;
    std::cout <<"L=" << L << ", M=" << M << ", roots=" << roots << std::endl;
    for (int i = 0; i < m; i++) {
        var_t res = roots[i] * (elem_t)R - 2.0L * L;
        for (int j = 0; j < m; j++) {
            if (j == i) continue;
            res += 2.0L * roots[i] * roots[j] / (roots[i] - roots[j]);
        }
        BOOST_TEST_INFO("m=" << m << ", i=" << i);
        BOOST_TEST(std::abs(res) < 1e-8);
    }
}

BOOST_AUTO_TEST_SUITE_END()
