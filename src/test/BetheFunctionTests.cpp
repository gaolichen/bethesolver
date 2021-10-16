//Link to Boost
#define BOOST_TEST_DYN_LINK

//VERY IMPORTANT - include this last
//#include <boost/test/included/unit_test.hpp>
#include <boost/test/unit_test.hpp>

#include "test.h"
#include "../BetheFunctions.h"
#include "../BethePolynomials.h"

using namespace std::literals::complex_literals;

// test suite
BOOST_FIXTURE_TEST_SUITE(BetheFactor_suite, SimpleTestFixture, * utf::label("BetheFactor"))

BOOST_AUTO_TEST_CASE(coef_test)
{
    int L = 6;
    int M = 2;
    BetheFactorLeft left(L, M);
    BetheFactorRight right(L, M);
    for (int i = 0; i < 3; i++) {
        var_t *coef1Left = left.coef1((RootType)i);
        var_t *coef1Right = right.coef1((RootType)i);
        for (int j = 0; j < 3; j++) {
            BOOST_TEST_INFO("i=" << i << ", j=" << j);
            BOOST_TEST(isZero(std::conj(coef1Left[j]) - coef1Right[j]));
        }
    }
    
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            var_t *coef2Left = left.coef2((RootType)i, (RootType)j);
            var_t *coef2Right = right.coef2((RootType)i, (RootType)j);
            for (int k = 0; k < 6; k++) {
                BOOST_TEST_INFO("i=" << i << ", j=" << j << ", k=" << k);
                BOOST_TEST(isZero(std::conj(coef2Left[k]) - coef2Right[k]));
            }
        }
    }
}


BOOST_DATA_TEST_CASE(eval_test_regular, bdata::random(4, 16) ^ bdata::xrange(20), L, index)
{
    int M = random64(L/2);
    if (M <= 1) M = 2;
    Vector roots(M);
    for (int i = 0; i < M; i++) {
        roots(i) = randomComplex(0.1, 0.9);
    }
        
    BetheSolution sol1(roots, L, 0.0, 0.0);
    InvertableSolution sol2(roots, L);
    
    BetheFactorLeft left1(L, M);
    BetheFactorRight right1(L, M);
    
    //BethePolynomial left2(L, M, true);
    BethePolynomialLeft left2(L, M);
    //BethePolynomial right2(L, M, false);
    BethePolynomialRight right2(L, M);
    
    
    BOOST_TEST((left1.eval(sol1) - left2.eval(sol2)).norm() == 0, tt::tolerance(eps));
    Vector v1 = right1.eval(sol1);
    Vector v2 = right2.eval(sol2);
    BOOST_TEST((v1 - v2).norm() == 0, tt::tolerance(eps));    
}

BOOST_DATA_TEST_CASE(eval_test_largest, bdata::random(4, 16) ^ bdata::xrange(20), L, index)
{
    int M = random64(L/2);
    if (M <= 1) M = 2;
    int indexOfLargest = random64(M - 1);    
    Vector roots(M);
    roots(indexOfLargest) = randomComplex(blowup_root_cutoff/0.3, blowup_root_cutoff/0.2);

    for (int i = 0; i < M; i++) {
        if (i == indexOfLargest) continue;
        roots(i) = randomComplex(0.1, 0.9);
    }
        
    BetheSolution sol1(roots, L, 0.0, 0.01);
    InvertableSolution sol2(roots, L);
    
    BetheFactorLeft left1(L, M);
    BetheFactorRight right1(L, M);
    
    BethePolynomialLeft left2(L, M);
    BethePolynomialRight right2(L, M);
    
    
    BOOST_TEST((left1.eval(sol1) - left2.eval(sol2)).norm() == 0, tt::tolerance(eps));
    Vector v1 = right1.eval(sol1);
    Vector v2 = right2.eval(sol2);
//    std::cout << "indexOfLargest=" << indexOfLargest << std::endl;
//    std::cout << "v1=" << v1 << std::endl;
//    std::cout << "v2=" << v2 << std::endl;
    BOOST_TEST((v1 - v2).norm() == 0, tt::tolerance(eps));    
}

BOOST_DATA_TEST_CASE(eval_test_largest_large, bdata::random(4, 16) ^ bdata::xrange(1), L, index)
{
    int M = 2;
    int indexOfLargest = 0;
    int indexOfLarge = 1;
    
    Vector roots(M);
    roots(indexOfLargest) = randomComplex(blowup_root_cutoff/0.3, blowup_root_cutoff/0.2);
    roots(indexOfLarge) = randomComplex(blowup_root_cutoff/0.5, blowup_root_cutoff/0.4);
        
    BetheSolution sol1(roots, L, 0.0, 0.01);
    InvertableSolution sol2(roots, L);
    
    BetheFactorLeft left1(L, M);
    BetheFactorRight right1(L, M);
    
    BethePolynomialLeft left2(L, M);
    BethePolynomialRight right2(L, M);
    
    Vector v1 = left1.eval(sol1);
    Vector v2 = left2.eval(sol2);
    
    var_t ratio = roots(indexOfLarge) / roots(indexOfLargest);
    v2(indexOfLargest) *= roots(indexOfLarge);
    v2(indexOfLarge) *= roots(indexOfLarge) * std::pow(ratio, L);
    
//    std::cout << "v1=" << v1 << std::endl;
//    std::cout << "v2=" << v2 << std::endl;

    BOOST_TEST((v1 - v2).norm() == 0, tt::tolerance(eps));
    v1 = right1.eval(sol1);
    v2 = right2.eval(sol2);
    v2(indexOfLargest) *= roots(indexOfLarge);
    v2(indexOfLarge) *= roots(indexOfLarge) * std::pow(ratio, L); 
//    std::cout << "indexOfLargest=" << indexOfLargest << std::endl;
//    std::cout << "v1=" << v1 << std::endl;
//    std::cout << "v2=" << v2 << std::endl;
    BOOST_TEST((v1 - v2).norm() == 0, tt::tolerance(eps));    
}

BOOST_DATA_TEST_CASE(eval_test_singular, bdata::random(2, 16) ^ bdata::xrange(20), halfL, index)
{
    int L = 2 * halfL;
    int M = 2;
    int index1 = 0;
    int index2 = 1;
    
    Vector roots(M);
    roots(index2) = -0.5il;
    roots(index1) = 0.5il;
    BetheSolution sol(roots, L, 0.0, 0.01);
        
    BetheFactorLeft left(L, M);
    BetheFactorRight right(L, M);
    
//    std::cout << "left.eval() = " << left.eval(sol) << std::endl;
//    std::cout << "right.eval() = " << right.eval(sol) << std::endl;
    
    BOOST_TEST((left.eval(sol) - right.eval(sol)).norm() == 0, tt::tolerance(eps));
}


BOOST_DATA_TEST_CASE(diff_test_regular, bdata::random(4, 16) ^ bdata::xrange(20), L, index)
{
    int M = random64(L/2);
    if (M <= 1) M = 2;
    Vector roots(M);
    for (int i = 0; i < M; i++) {
        roots(i) = randomComplex(0.1, 0.9);
    }
        
    BetheSolution sol1(roots, L, 0.0, 0.01);
    InvertableSolution sol2(roots, L);
    
    BetheFactorLeft left1(L, M);
    BetheFactorRight right1(L, M);
    
    BethePolynomialLeft left2(L, M);
    BethePolynomialRight right2(L, M);
    
    BOOST_TEST((left1.diff(sol1) - left2.diff(sol2)).norm() == 0, tt::tolerance(eps));
    Matrix v1 = right1.diff(sol1);
    Matrix v2 = right2.diff(sol2);
//    std::cout << "v1-v2=" << v1 - v2 << std::endl;
    BOOST_TEST((v1 - v2).norm() == 0, tt::tolerance(eps));    
}

BOOST_DATA_TEST_CASE(diff_test_largest, bdata::random(4, 16) ^ bdata::xrange(20), L, index)
{
    int M = random64(L/2);
    if (M <= 1) M = 2;
    int indexOfLargest = random64(M - 1);    
    Vector roots(M);
    roots(indexOfLargest) = randomComplex(blowup_root_cutoff/0.3, blowup_root_cutoff/0.2);

    for (int i = 0; i < M; i++) {
        if (i == indexOfLargest) continue;
        roots(i) = randomComplex(0.1, 0.9);
    }
        
    BetheSolution sol1(roots, L, 0.0, 0.01);
    InvertableSolution sol2(roots, L);
    
    BetheFactorLeft left1(L, M);
    BetheFactorRight right1(L, M);
    
    BethePolynomialLeft left2(L, M);
    BethePolynomialRight right2(L, M);
    
    Matrix w1 = left1.diff(sol1);
    Matrix w2 = left2.diff(sol2);
//    std::cout << "M=" << M << ", indexOfLargest=" << indexOfLargest << std::endl;
//    std::cout << "w1=" << w1 << std::endl;
//    std::cout << "w2=" << w2 << std::endl;
//    std::cout << "diff=" << w1-w2 << std::endl;

    
    BOOST_TEST((w1 - w2).norm() == 0, tt::tolerance(eps));
    Matrix v1 = right1.diff(sol1);
    Matrix v2 = right2.diff(sol2);
//    std::cout << "M=" << M << ", indexOfLargest=" << indexOfLargest << std::endl;
//    std::cout << "v1=" << v1 << std::endl;
//    std::cout << "v2=" << v2 << std::endl;
//    std::cout << "diff=" << v1-v2 << std::endl;
    BOOST_TEST((v1 - v2).norm() == 0, tt::tolerance(eps));    
}

BOOST_DATA_TEST_CASE(diff_test_LeftFactor, bdata::random(4, 16) ^ bdata::xrange(20), L, index)
{
    int M = random64(L/2);
    if (M <= 1) M = 2;
    int indexOfLargest = random64(M - 1);    
    Vector roots(M);
    Vector diff1(M); // infinetesimal differential of roots for regular roots
    Vector diff2(M); // infinetesimal differential of roots for BetheSolution.root();
    elem_t delta = 1e-4;
    
    roots(indexOfLargest) = randomComplex(blowup_root_cutoff/0.3, blowup_root_cutoff/0.2);
    diff1(indexOfLargest) = delta * (random64(3) + 1) * randomComplex(1.0, 1.1);
    diff2(indexOfLargest) = -diff1(indexOfLargest) / (roots(indexOfLargest) * roots(indexOfLargest));

    for (int i = 0; i < M; i++) {
        if (i == indexOfLargest) continue;
        diff1(i) = delta * (random64(3) + 1) * randomComplex(1.0, 1.1);
        if (i % 2 == 0) {
            // large root
            roots(i) = randomComplex(blowup_root_cutoff/0.5, blowup_root_cutoff/0.4);
            diff2(i) = diff1(i) / roots(indexOfLargest) - diff1(indexOfLargest) * roots(i) / (roots(indexOfLargest) * roots(indexOfLargest));
        } else {
            // regular root
            roots(i) = randomComplex(blowup_root_cutoff/1.5, blowup_root_cutoff/1.2);
            diff2(i) = diff1(i);
        }
    }
    
    Vector roots2 = roots + diff1;
    Vector roots3 = roots + diff1 * 0.5L;
            
    BetheSolution sol1(roots, L, 0.0, 0.01);
    BetheSolution sol2(roots2, L, 0.0, 0.01);
    BetheSolution sol3(roots3, L, 0.0, 0.01);
    
    BetheFactorLeft left(L, M);
    Vector res1 = left.eval(sol1);
    Vector res2 = left.eval(sol2);
    Vector resDiff = res2 - res1;
    
    Matrix jacobian = left.diff(sol3);
    Vector resDiff2 = jacobian * diff2;
    
//    std::cout << "resDiff=" << resDiff << std::endl;
//    std::cout << "resDiff2=" << resDiff2 << std::endl;
    for (int i = 0; i < M; i++) {
        var_t v1 = resDiff(i);
        var_t v2 = resDiff2(i);
        BOOST_TEST_INFO("M=" << M << ",i=" << i << ", v1=" << v1 << ", v2=" << v2);
        BOOST_TEST(std::abs(v1-v2)/std::max(std::abs(v1), std::abs(v2)) == 0, tt::tolerance(delta * 10));
    }
}

BOOST_DATA_TEST_CASE(diff_test_RightFactor, bdata::random(4, 16) ^ bdata::xrange(20), L, index)
{
    int M = random64(L/2);
    if (M <= 1) M = 2;
    int indexOfLargest = random64(M - 1);    
    Vector roots(M);
    Vector diff1(M); // infinetesimal differential of roots for regular roots
    Vector diff2(M); // infinetesimal differential of roots for BetheSolution.root();
    elem_t delta = 1e-4;
    
    roots(indexOfLargest) = randomComplex(blowup_root_cutoff/0.3, blowup_root_cutoff/0.2);
    diff1(indexOfLargest) = delta * (random64(3) + 1) * randomComplex(1.0, 1.1);
    diff2(indexOfLargest) = -diff1(indexOfLargest) / (roots(indexOfLargest) * roots(indexOfLargest));

    for (int i = 0; i < M; i++) {
        if (i == indexOfLargest) continue;
        diff1(i) = delta * (random64(3) + 1) * randomComplex(1.0, 1.1);
        if (i % 2 == 0) {
            // large root
            roots(i) = randomComplex(blowup_root_cutoff/0.5, blowup_root_cutoff/0.4);
            diff2(i) = diff1(i) / roots(indexOfLargest) - diff1(indexOfLargest) * roots(i) / (roots(indexOfLargest) * roots(indexOfLargest));
        } else {
            // regular root
            roots(i) = randomComplex(blowup_root_cutoff/1.5, blowup_root_cutoff/1.2);
            diff2(i) = diff1(i);
        }
    }
    
    Vector roots2 = roots + diff1;
    Vector roots3 = roots + diff1 * 0.5L;
            
    BetheSolution sol1(roots, L, 0.0, 0.01);
    BetheSolution sol2(roots2, L, 0.0, 0.01);
    BetheSolution sol3(roots3, L, 0.0, 0.01);
    
    BetheFactorRight right(L, M);
    Vector res1 = right.eval(sol1);
    Vector res2 = right.eval(sol2);
    Vector resDiff = res2 - res1;
    
    Matrix jacobian = right.diff(sol3);
    Vector resDiff2 = jacobian * diff2;
    
//    std::cout << "resDiff=" << resDiff << std::endl;
//    std::cout << "resDiff2=" << resDiff2 << std::endl;
    for (int i = 0; i < M; i++) {
        var_t v1 = resDiff(i);
        var_t v2 = resDiff2(i);
        BOOST_TEST_INFO("M=" << M << ",i=" << i << ", v1=" << v1 << ", v2=" << v2);
        BOOST_TEST(std::abs(v1-v2)/std::max(std::abs(v1), std::abs(v2)) == 0, tt::tolerance(delta * 20));
    }
}

BOOST_AUTO_TEST_CASE(diff_test_singular)
{
    elem_t tmp = eps_may_be_singular;
    eps_may_be_singular = 0.1L;
    int L = 6;
    int M = 3;
    var_t r1(0.54236679461751705L,2.90957067279989408e-11L);
    var_t r2(0.0968423125274162171L,0.500001874059868433L);
    var_t r3(0.0968423125273995912L,-0.50000187418704665L);
    elem_t beta = 0.914203L;

    Vector roots(M);
    roots(0) = r1;
    roots(1) = r2;
    roots(2) = r3;
    BetheSolution sol(roots, L, beta, Pi/1000);
    BOOST_TEST_REQUIRE(sol.indexOfSingularRoot1() >= 1);
    BOOST_TEST_REQUIRE(sol.indexOfSingularRoot2() >= 1);

    BetheFactorLeft left(L, M);
    Vector res1 = left.eval(sol);
    Matrix jacobian = left.diff(sol);
    
    Vector inc(M);
    inc(0) = 0.0L;
    inc(sol.indexOfSingularRoot1()) = 0.0L;
    inc(sol.indexOfSingularRoot2()) = 0.0L;
    inc(sol.indexOfSingularRoot2()) = 0.0L;
    inc(0) = randomComplex(1e-3, 2e-3);
    inc(sol.indexOfSingularRoot1()) = randomComplex(1e-1, 2e-1);
    inc(sol.indexOfSingularRoot2()) = randomComplex(1e-3, 2e-3);
    
    sol.root() += inc;
    Vector res2 = left.eval(sol);
    Vector diff = res2 - res1;
        
    Vector expected = jacobian * inc;
//    std::cout << "res1=" << res1 << std::endl;
//    std::cout << "res2=" << res2 << std::endl;
//    std::cout << "diff=" << diff << std::endl;
//    std::cout << "expected=" << expected << std::endl;
    
    BOOST_TEST_INFO("inc=" << inc);
    BOOST_TEST((diff-expected).norm()/expected.norm() < 2e-2);

    eps_may_be_singular = tmp;    
}

BOOST_AUTO_TEST_CASE(diff_test_singular_RightFactor)
{
    elem_t tmp = eps_may_be_singular;
    eps_may_be_singular = 0.1L;
    int L = 6;
    int M = 3;
    var_t r1(0.54236679461751705L,2.90957067279989408e-11L);
    var_t r2(0.0968423125274162171L,0.500001874059868433L);
    var_t r3(0.0968423125273995912L,-0.50000187418704665L);
    elem_t beta = 0.914203L;

    Vector roots(M);
    roots(0) = r1;
    roots(1) = r2;
    roots(2) = r3;
    BetheSolution sol(roots, L, beta, Pi/1000);
    BOOST_TEST_REQUIRE(sol.indexOfSingularRoot1() >= 1);
    BOOST_TEST_REQUIRE(sol.indexOfSingularRoot2() >= 1);

    BetheFactorRight right(L, M);
    Vector res1 = right.eval(sol);
    Matrix jacobian = right.diff(sol);
    
    Vector inc(M);
    inc(0) = 0.0L;
    inc(sol.indexOfSingularRoot1()) = 0.0L;
    inc(sol.indexOfSingularRoot2()) = 0.0L;
    inc(sol.indexOfSingularRoot2()) = 0.0L;
    inc(0) = randomComplex(1e-3, 2e-3);
    inc(sol.indexOfSingularRoot1()) = randomComplex(1e-1, 2e-1);
    inc(sol.indexOfSingularRoot2()) = randomComplex(1e-3, 2e-3);
    
    sol.root() += inc;
    Vector res2 = right.eval(sol);
    Vector diff = res2 - res1;
        
    Vector expected = jacobian * inc;
//    std::cout << "res1=" << res1 << std::endl;
//    std::cout << "res2=" << res2 << std::endl;
//    std::cout << "diff=" << diff << std::endl;
//    std::cout << "expected=" << expected << std::endl;
    
    BOOST_TEST_INFO("inc=" << inc);
    BOOST_TEST((diff-expected).norm()/expected.norm() < 2e-2);

    eps_may_be_singular = tmp;    
}

BOOST_AUTO_TEST_CASE(diff_test_largest_LeftFactor)
{
    int L = 6;
    int M = 3;
    
    var_t r1(0.288675134570472836L, -1.98390501120640128e-11L);
    var_t r2(4833488737.1251517L, 8727942389.24525939L);
    var_t r3(-5137763027.0882383L, 8552377948.96713318L);
    elem_t beta = 1.5708L;

    Vector roots(M);
    roots(0) = r1;
    roots(1) = r2;
    roots(2) = r3;
    BetheSolution sol(roots, L, beta, Pi/1000);
    BOOST_TEST_REQUIRE(sol.indexOfBlowupRoot() >= 1);

    BetheFactorRight left(L, M);
    Vector res1 = left.eval(sol);
    Matrix jacobian = left.diff(sol);
    
    Vector inc = Vector::Zero(M);
    for (int i = 0; i < inc.size(); i++) {
        if (sol.type(i) == regular) {
            inc(i) = randomComplex(1e-4, 2e-4);
        } else if (sol.type(i) == largest) {
            elem_t phi = std::arg(sol.root(i));
            inc(i) = -std::abs(sol.root(i)) * random(1.0, 1.5) * std::polar(1.0L, phi);;
        } else {
            inc(i) = randomComplex(1e-4, 2e-4);
        }
    }
    
    sol.root() += inc;
    Vector res2 = left.eval(sol);
    Vector diff = res2 - res1;
        
    Vector expected = jacobian * inc;
//    std::cout << "res1=" << res1 << std::endl;
//    std::cout << "res2=" << res2 << std::endl;
//    std::cout << "diff=" << diff << std::endl;
//    std::cout << "expected=" << expected << std::endl;
    
    BOOST_TEST_INFO("inc=" << inc);
    BOOST_TEST((diff-expected).norm()/expected.norm() < 1e-3);
}


BOOST_AUTO_TEST_CASE(diff_test_largest_RightFactor)
{
    int L = 6;
    int M = 3;
    
    var_t r1(0.288675134570472836L, -1.98390501120640128e-11L);
    var_t r2(4833488737.1251517L, 8727942389.24525939L);
    var_t r3(-5137763027.0882383L, 8552377948.96713318L);
    elem_t beta = 1.5708L;

    Vector roots(M);
    roots(0) = r1;
    roots(1) = r2;
    roots(2) = r3;
    BetheSolution sol(roots, L, beta, Pi/1000);
    BOOST_TEST_REQUIRE(sol.indexOfBlowupRoot() >= 1);

    BetheFactorRight right(L, M);
    Vector res1 = right.eval(sol);
    Matrix jacobian = right.diff(sol);
    
    Vector inc = Vector::Zero(M);
    for (int i = 0; i < inc.size(); i++) {
        if (sol.type(i) == regular) {
            inc(i) = randomComplex(1e-4, 2e-4);
        } else if (sol.type(i) == largest) {
            elem_t phi = std::arg(sol.root(i));
            inc(i) = std::abs(sol.root(i)) * random(1.0, 1.5) * std::polar(1.0L, phi);;
        } else {
            inc(i) = randomComplex(1e-4, 2e-4);
        }
    }
    
    sol.root() += inc;
    Vector res2 = right.eval(sol);
    Vector diff = res2 - res1;
        
    Vector expected = jacobian * inc;
//    std::cout << "res1=" << res1 << std::endl;
//    std::cout << "res2=" << res2 << std::endl;
//    std::cout << "diff=" << diff << std::endl;
//    std::cout << "expected=" << expected << std::endl;
    
    BOOST_TEST_INFO("inc=" << inc);
    BOOST_TEST((diff-expected).norm()/expected.norm() < 1e-3);
}
/*
BOOST_AUTO_TEST_CASE(diff_test_largest2)
{
    int L = 6;
    int M = 2;
    var_t r1(-2002453610.04371309,18872186.6385834396);
    var_t r2(-2002453610.04371309,18872186.6385834396);
    elem_t beta = 0.785398;
    elem_t delta = Pi/1000;

    Vector roots(M);
    roots(0) = r1;
    roots(1) = r2;
    BetheSolution sol(roots, L, beta, delta);
    BOOST_TEST_REQUIRE(sol.indexOfBlowupRoot() >= 0);

    BetheFactorLeft left(L, M);
    BetheFactorRight right(L, M);
    Vector res1 = right.eval(sol);
    Matrix jacobian1 = left.diff(sol);
    Matrix jacobian2 = right.diff(sol);
    Matrix jacobian = jacobian1 - jacobian2 * std::polar(1.0L, 2.0L * L * beta);
    Matrix jacobianInv = jacobian.inverse();

    Vector rhs = left.eval(sol) * var_t((elem_t)0.0, 2.0 * L);
    
//    std::cout << "jacobian = " << jacobian << std::endl;
//    std::cout << "jacobianInv=" << jacobianInv << std::endl;
//    std::cout << "rhs=" << std::setprecision(18) << rhs << std::endl;
    
    Vector inc1 = jacobianInv * rhs * delta/2.0;
//    std::cout << "jacobianInv.rhs=" << std::setprecision(18) << inc1 << std::endl;
//    std::cout << "root=" << sol.root() << std::endl;
    sol.root() += inc1;
//    std::cout << "after root=" << sol.normalRoot() << std::endl;
}*/

BOOST_AUTO_TEST_SUITE_END()

