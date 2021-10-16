//Link to Boost
#define BOOST_TEST_DYN_LINK

//VERY IMPORTANT - include this last
//#include <boost/test/included/unit_test.hpp>
#include <boost/test/unit_test.hpp>

#include "test.h"
#include "../BetheHomotopy.h"

using namespace std::literals::complex_literals;

// test suite
BOOST_FIXTURE_TEST_SUITE(BetheHomotopy_suite, SimpleTestFixture, * utf::label("BetheHomotopy"))

BOOST_AUTO_TEST_CASE(diff_test_largest)
{
    int L = 6;
    int M = 2;
    var_t r1(-2002453610.04371309,18872186.6385834396);
    var_t r2(-2002453610.04371309,18872186.6385834396);
//    elem_t beta = 0.785398163397448277L;
    elem_t beta = 0.7853981634L;
    elem_t delta = Pi/1000;

    Vector roots(M);
    roots(0) = r1;
    roots(1) = r2;
    BetheSolution sol(roots, L, 0.0, delta);
    BOOST_TEST_REQUIRE(sol.indexOfBlowupRoot() >= 0);
//    std::cout << "sol.root()=" << std::setprecision(18) << sol.root() << std::endl;
//    std::cout << "sol.normalRoot()=" << std::setprecision(18) << sol.normalRoot() << std::endl;
    
    PoleFreeBetheHomotopy homotopy(L, M, Pi);
    Vector old = sol.root();
    Vector k1 = homotopy.impliciteDiff(beta, sol);
    std::cout << "k1=" << k1 << std::endl;
    sol.root() += k1 * (delta/2);
    Vector k2 = homotopy.impliciteDiff(beta + delta/2, sol);
    std::cout << "k2=" << k2 << std::endl;
    sol.setRoot(old + k2 * (delta/2));
    Vector k3 = homotopy.impliciteDiff(beta + delta/2, sol);
    std::cout << "k3=" << k3 << std::endl;
    sol.setRoot(old + k3 * (delta/2));
    Vector k4 = homotopy.impliciteDiff(beta + delta, sol);
    std::cout << "k4=" << k4 << std::endl;
    
    Vector change = k1 * (delta/6) + k2 * (delta/3) + k3 * (delta/3) + k4 * (delta/6);
    std::cout << "change=" << change << std::endl;
    sol.root() = old + change;
    std::cout << "sol.root()=" << sol.root() << std::endl;
    std::cout << "sol.normalRoot()=" << std::setprecision(18) << sol.normalRoot() << std::endl;
}

var_t getA(BetheSolution& sol) {
    var_t ret = 1.0L/(sol.root(0) + 0.5il) - 1.0L/(sol.root(0) - 0.5il);
    return ret * (elem_t)sol.chainLength();
}

var_t getB(BetheSolution& sol) {
    return 1.0L/(sol.root(0) - sol.root(1) - 1.0il) - 1.0L/(sol.root(0) - sol.root(1) + 1.0il);
}

var_t det(BetheSolution& sol) {
    var_t A = getA(sol);
    var_t B = getB(sol);
    return B * std::conj(B) - (A+B) * std::conj(A+B);
}

var_t eps2(BetheSolution& sol) {
    var_t A = getA(sol);
    var_t B = getB(sol);
    return (elem_t)sol.chainLength() * 2.0il * (A + B - std::conj(B)) / det(sol);
}

var_t eps1(BetheSolution& sol) {
    var_t A = getA(sol);
    var_t B = getB(sol);
    return -(elem_t)sol.chainLength() * 2.0il * (std::conj(A + B) - B) / det(sol);
}


BOOST_AUTO_TEST_CASE(diff_test_conjugate_roots)
{
    int L = 7;
    int M = 2;
    var_t r1(-1.21319620301732725,0.0943445611962252633);
    var_t r2(-1.21319620301732725,-0.0943445611962252632);
    elem_t beta = 1.630486587213102690762;
    elem_t delta = Pi/1000;

    Vector roots(M);
    roots(0) = r1;
    roots(1) = r2;
    BetheSolution sol(roots, L, 0.0, delta);
//    std::cout << "sol.root()=" << std::setprecision(18) << sol.root() << std::endl;
//    std::cout << "sol.normalRoot()=" << std::setprecision(18) << sol.normalRoot() << std::endl;
    
    var_t e1, e2;
    PoleFreeBetheHomotopy homotopy(L, M, Pi);
    Vector old = sol.root();
    e1 = eps1(sol);
    e2 = eps2(sol);
    std::cout << "det=" << det(sol) << std::endl;
    std::cout << "expected jacobian= " << det(sol) * homotopy.Left()->eval(sol) * homotopy.Right()->eval(sol) * std::polar(1.0L, 2.0 * L * beta) << std::endl;
    Vector k1 = homotopy.impliciteDiff(beta, sol);
//    Matrix jacobian = homotopy.diff(beta, sol);
//    std::cout << "jacobian=" << jacobian << std::endl;
    std::cout << "k1=" << k1 << std::endl;
    std::cout << "(e1,e2)=" << e1 << ' ' << e2 << std::endl;
    sol.root() += k1 * (delta/2);
    e1 = eps1(sol);
    e2 = eps2(sol);

    std::cout << "det=" << det(sol) << std::endl;
    std::cout << "expected jacobian= " << det(sol) * homotopy.Left()->eval(sol) * homotopy.Right()->eval(sol) * std::polar(1.0L, 2.0 * L * (beta + delta/2)) << std::endl;
//    jacobian = homotopy.diff(beta + delta / 2, sol);
//    std::cout << "jacobian=" << jacobian << std::endl;
    Vector k2 = homotopy.impliciteDiff(beta + delta/2, sol);
    std::cout << "k2=" << k2 << std::endl;
    std::cout << "(e1,e2)=" << e1 << ' ' << e2 << std::endl;
    sol.setRoot(old + k2 * (delta/2));
    e1 = eps1(sol);
    e2 = eps2(sol);
//    std::cout << "det=" << det(sol) << std::endl;
//    jacobian = homotopy.diff(beta + delta / 2, sol);
//    std::cout << "jacobian=" << jacobian << std::endl;
    Vector k3 = homotopy.impliciteDiff(beta + delta/2, sol);
    std::cout << "k3=" << k3 << std::endl;
    std::cout << "(e1,e2)=" << e1 << ' ' << e2 << std::endl;
    sol.setRoot(old + k3 * (delta/2));
    e1 = eps1(sol);
    e2 = eps2(sol);
    //std::cout << "det=" << det(sol) << std::endl;
//    jacobian = homotopy.diff(beta + delta, sol);
//    std::cout << "jacobian=" << jacobian << std::endl;
    Vector k4 = homotopy.impliciteDiff(beta + delta, sol);
//    std::cout << "k4=" << k4 << std::endl;
//    std::cout << "(e1,e2)=" << e1 << ' ' << e2 << std::endl;
    
    Vector change = k1 * (delta/6) + k2 * (delta/3) + k3 * (delta/3) + k4 * (delta/6);
//    std::cout << "old=" << old << std::endl;
//    std::cout << "change=" << change << std::endl;
    sol.root() = old + change;
//    std::cout << "sol.root()=" << sol.root() << std::endl;
//    std::cout << "sol.normalRoot()=" << std::setprecision(18) << sol.normalRoot() << std::endl;
}

BOOST_AUTO_TEST_CASE(diff_test_largest_roots)
{
    int L = 8;
    int M = 3;
    //1.77710281350482	(0.162962788930455523,-7.73479176308443192e-13)	(-1.50005668874134187,-2.50992222081144617e-12)	(0.162962789055243649,1.05111478538827497e-12)	1.39947

    var_t r1(0.162962788930455523L,-7.73479176308443192e-13L);
    var_t r2(-1.50005668874134187L,-2.50992222081144617e-12L);
    var_t r3(0.162962789055243649L,1.05111478538827497e-12L);
    elem_t beta = 1.77710281350482;
    elem_t totalPhase = 3.66417L;
    elem_t delta = Pi/1000;

    Vector roots(M);
    roots(0) = r1;
    roots(1) = r2;
    roots(2) = r3;
    SET_LOG_LEVEL(Info);
    BetheSolution sol(roots, L, totalPhase, delta);
    std::cout << "sol.root()=" << std::setprecision(18) << sol.root() << std::endl;
    std::cout << "sol.normalRoot()=" << std::setprecision(18) << sol.normalRoot() << std::endl;
    Vector startRoots(M);
    startRoots(0) = -0.5il;
    startRoots(1) = 0.0L;
    startRoots(2) = 0.5il;
    
    PoleFreeBetheHomotopy homotopy(L, M, Pi);
    homotopy.setStartRoot(startRoots);
    std::cout << homotopy.error(beta, sol) << std::endl;
    sol.output();
}


BOOST_AUTO_TEST_SUITE_END()

