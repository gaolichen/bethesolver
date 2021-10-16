//Link to Boost
#define BOOST_TEST_DYN_LINK

//VERY IMPORTANT - include this last
//#include <boost/test/included/unit_test.hpp>
#include <boost/test/unit_test.hpp>

#include "test.h"
#include "../HomotopyContinuation.h"
#include "../BetheHomotopy.h"

// test suite
BOOST_FIXTURE_TEST_SUITE(Demo_suite, SimpleTestFixture, * utf::label("UnityStartSystem"))

/*
BOOST_AUTO_TEST_CASE(Demo2)
{
    SolvedBetheEquation start(6, 2);
    BetaHomotopy hom(0.5235, &start);
    //BetaHomotopy hom(0.2578, &start);
    //BetaHomotopy hom(0.25502, &start);
    HomotopyContinuation hc(1000);
    
    std::cout << "begin demo2" << std::endl;
    hc.solve(hom);
    for (int i = 0; i < hom.numberOfRoots(); i++) {
        Vector root = hom.getRoot(i);
        std::cout << "demo2 root=" << root << ", eval=" << hom.evalTarget(root) << std::endl;
    }
    std::cout << "end demo2" << std::endl;    
}*/

var_t normalPhase(var_t phase) {
    if (phase.real() < .0) {
        phase += 2 * Pi;
    } else {
        phase -= floor(phase.real() / (2.0 * Pi)) * 2 * Pi;
    }
    return phase;
}

BOOST_AUTO_TEST_CASE(Demo3, * utf::disabled())
{
//    ExpandedBetaHomotopy hom(0.523599 * 2, &numerator, &denominator);
    elem_t beta = Pi;
    elem_t expectedMomentumDiff = .0;
    ExpandedBetaHomotopy hom(6, 2, beta);
    hom.setSteps(1000);
    BetheRootCache cache(6, 2);
    BetheHomotopyContinuation hc;
    
    std::cout << "begin demo3" << std::endl;
    for (int i = 0; i < cache.numberOfRoots(); i++) {
        std::string file = "BS_L6M2r" + ToString(i) + ".txt";
        hc.setTraceFile(file);
        Vector startRoot = cache.getRoot(i);
        InvertableSolution sol(startRoot, 6);
        hc.solve(hom, sol);
        elem_t m1 = momentum(startRoot);
        elem_t m2 = momentum(sol.normalRoot());
        BOOST_TEST_INFO("m1/Pi=" << m1/(elem_t)Pi << ", m2/Pi=" << m2/(elem_t)Pi);
        BOOST_TEST(std::abs(floatMod(m2 - m1, 2.0 * Pi) - expectedMomentumDiff) < 1e-4);
        std::cout << "from " << startRoot << std::endl;
        std::cout << "to " << sol.normalRoot() << ", error=" << hom.evalTarget(sol) << std::endl;
    }
    std::cout << "end demo3" << std::endl;    
}


BOOST_AUTO_TEST_CASE(Demo4, * utf::disabled())
{
    ExpandedBetaHomotopy hom(6, 3, 0.5235);
    hom.setSteps(1000);
    BetheRootCache cache(6, 3);
    BetheHomotopyContinuation hc;
    
    std::cout << "begin demo4" << std::endl;
    for (int i = 0; i < cache.numberOfRoots(); i++) {
        std::string file = "BS_L6M3r" + ToString(i) + ".txt";
        hc.setTraceFile(file);
        Vector startRoot = cache.getRoot(i);
        InvertableSolution sol(startRoot, 6);
        hc.solve(hom, sol);
        std::cout << "startRoot=" << startRoot << std::endl;
//        BOOST_TEST_INFO("root=" << sol.normalRoot() << ", startRoot=" << startRoot);
//        BOOST_TEST((sol.normalRoot()-startRoot).norm() < 1e-3);
        std::cout << "root=" << sol.normalRoot() << ", eval=" << hom.evalTarget(sol).norm() << std::endl;
    }
    std::cout << "end demo4" << std::endl;    
}

BOOST_AUTO_TEST_CASE(Demo5, * utf::disabled())
{
    ExpandedBetaHomotopy hom(7, 2, 0.4488);
    hom.setSteps(1000);
    BetheRootCache cache(7, 2);
    BetheHomotopyContinuation hc;
    
    std::cout << "begin demo5" << std::endl;
    for (int i = 0; i < cache.numberOfRoots(); i++) {
        std::string file = "BS_L7M2r" + ToString(i) + ".txt";
        hc.setTraceFile(file);
        Vector startRoot = cache.getRoot(i);
        InvertableSolution sol(startRoot, 7);
        hc.solve(hom, sol);
//        BOOST_TEST_INFO("root=" << sol.normalRoot() << ", startRoot=" << startRoot);
//        BOOST_TEST((sol.normalRoot() - startRoot).norm() < 1e-3);
        std::cout << "demo5 root=" << sol.normalRoot() << ", eval=" << hom.evalTarget(sol);
        std::cout << ", momentum="<< momentum(sol.normalRoot()) <<  std::endl;
        std::cout << "start root=" << startRoot;
        std::cout << ", momentum="<< momentum(startRoot) << std::endl << std::endl;
    }
    std::cout << "end demo5" << std::endl;    
}

BOOST_AUTO_TEST_CASE(Demo6, * utf::disabled())
{
    int L = 6;
    int M = 2;
    int steps = 1000;
    elem_t totalPhase = Pi;
    elem_t expectedMomentumDiff = .0;
    PoleFreeBetheHomotopy hom(L, M, totalPhase);
    hom.setSteps(steps);
    BetheRootCache cache(L, M);
    BetheHomotopyContinuation hc;
    
    std::cout << "begin demo6" << std::endl;
    for (int i = 0; i < cache.numberOfRoots(); i++) {
        std::string file = "Demo6_BS_L" + ToString(L)+ "M" + ToString(M)+ "r" + ToString(i) + ".txt";
        hc.setTraceFile(file);
        Vector startRoot = cache.getRoot(i);
        BetheSolution sol(startRoot, L, 0.0, totalPhase / steps);
        hc.solve(hom, sol);
        elem_t m1 = momentum(startRoot);
        elem_t m2 = momentum(sol.normalRoot());
        BOOST_TEST_INFO("m1/Pi=" << m1/(elem_t)Pi << ", m2/Pi=" << m2/(elem_t)Pi);
        BOOST_TEST(std::abs(floatMod(m2 - m1, 2.0 * Pi) - expectedMomentumDiff) < 1e-4);
        std::cout << "from " << startRoot << std::endl;
        std::cout << "to " << sol.normalRoot() << ", error=" << hom.evalTarget(sol) << std::endl;
    }
    std::cout << "end demo6" << std::endl;    
}

BOOST_AUTO_TEST_CASE(Demo7, * utf::disabled())
{
    int L = 6;
    int M = 3;
    int steps = 1000;
    elem_t totalPhase = Pi;
    elem_t expectedMomentumDiff = .0;
    PoleFreeBetheHomotopy hom(L, M, totalPhase);
    hom.setSteps(steps);
    BetheRootCache cache(L, M);
    BetheHomotopyContinuation hc;
    
    std::cout << "begin demo7" << std::endl;
    for (int i = 0; i < cache.numberOfRoots(); i++) {
        std::string file = "Demo7_BS_L" + ToString(L)+ "M" + ToString(M)+ "r" + ToString(i) + ".txt";
        hc.setTraceFile(file);
        Vector startRoot = cache.getRoot(i);
        BetheSolution sol(startRoot, L, 0.0, totalPhase / steps);
        hc.solve(hom, sol);
        elem_t m1 = momentum(startRoot);
        elem_t m2 = momentum(sol.normalRoot());
        BOOST_TEST_INFO("m1/Pi=" << m1/(elem_t)Pi << ", m2/Pi=" << m2/(elem_t)Pi);
        BOOST_TEST(std::abs(floatMod(m2 - m1, 2.0 * Pi) - expectedMomentumDiff) < 1e-4);
        std::cout << "from " << startRoot << std::endl;
        std::cout << "to " << sol.normalRoot() << ", error=" << hom.evalTarget(sol) << std::endl;
    }
    std::cout << "end demo7" << std::endl;    
}

BOOST_AUTO_TEST_CASE(Demo8, * utf::disabled())
{
    int L = 7;
    int M = 2;
    int steps = 3000;
    elem_t totalPhase = Pi;
    elem_t expectedMomentumDiff = .0;
    PoleFreeBetheHomotopy hom(L, M, totalPhase);
    hom.setSteps(steps);
    BetheRootCache cache(L, M);
    BetheHomotopyContinuation hc;
    SET_LOG_LEVEL(Info);
    
    // i = 1, 3, 8 failed.
    // i = 1 fail at two identical infinity roots
    // i = 3 and 8 fail at two conjugate finite roots.
    std::cout << "begin demo8" << std::endl;
    for (int i = 1; i < cache.numberOfRoots() && i < 2; i++) {
        std::string file = "Demo8_BS_L" + ToString(L)+ "M" + ToString(M)+ "r" + ToString(i) + ".txt";
        hc.setTraceFile(file);
        Vector startRoot = cache.getRoot(i);
        BetheSolution sol(startRoot, L, 0.0, totalPhase / hom.getSteps());
        hc.solve(hom, sol);
        elem_t m1 = momentum(startRoot);
        elem_t m2 = momentum(sol.normalRoot());
        BOOST_TEST_INFO("i=" << i << ", m1/Pi=" << m1/(elem_t)Pi << ", m2/Pi=" << m2/(elem_t)Pi);
        BOOST_TEST(std::abs(floatMod(m2 - m1, 2.0 * Pi) - expectedMomentumDiff) < 1e-4);
        std::cout << "from " << startRoot << std::endl;
        std::cout << "to " << sol.normalRoot() << ", error=" << hom.evalTarget(sol) << std::endl;
    }
    std::cout << "end demo8" << std::endl;
    SET_LOG_LEVEL(Error);
}

BOOST_AUTO_TEST_CASE(Demo9, * utf::disabled())
{
    int L = 7;
    int M = 3;
    int steps = 6000;
    elem_t totalPhase = Pi;
    elem_t expectedMomentumDiff = .0;
    PoleFreeBetheHomotopy hom(L, M, totalPhase);
    hom.setSteps(steps);
    BetheRootCache cache(L, M);
    BetheHomotopyContinuation hc;
    
    // i = 0 fail at 2 identical infinite roots
    // i = 1 fail at 2 identical infinite roots
    // i = 4 fail at 2 conjugate finite roots
    // i=10 fail at 2 identical infinities
    // i = 11 fail at 2 identical infinities
    std::cout << "begin demo9" << std::endl;
    std::vector<bool> hasInfinity;
    for (int i = 6; i < cache.numberOfRoots() && i < 7; i++) {
        std::string file = "BS_L" + ToString(L)+ "M" + ToString(M)+ "r" + ToString(i) + ".txt";
        hc.setTraceFile(file);
        hc.hasInfinity = false;
        Vector startRoot = cache.getRoot(i);
        BetheSolution sol(startRoot, L, 0.0, totalPhase / steps);
        hc.solve(hom, sol);
        hasInfinity.push_back(hc.hasInfinity);
        elem_t m1 = momentum(startRoot);
        elem_t m2 = momentum(sol.normalRoot());
        BOOST_TEST_INFO("i=" << i << ", m1/Pi=" << m1/(elem_t)Pi << ", m2/Pi=" << m2/(elem_t)Pi);
        BOOST_TEST(std::abs(floatMod(m2 - m1, 2.0 * Pi) - expectedMomentumDiff) < 1e-4);
        std::cout << "from " << startRoot << std::endl;
        std::cout << "to " << sol.normalRoot() << ", error=" << hom.evalTarget(sol) << std::endl;
    }
    std::cout << hasInfinity << std::endl;
    std::cout << "end demo9" << std::endl;    
}

BOOST_AUTO_TEST_CASE(Demo10a, * utf::disabled())
{
    int L = 8;
    int M = 3;
    elem_t totalPhase = random(1.0, 2.0) * Pi;
    elem_t expectedMomentumDiff = floatMod(2 * M * totalPhase, 2 * Pi);
    PoleFreeBetheHomotopy hom(L, M, totalPhase);
    BetheRootCache cache(L, M);
    BetheHomotopyContinuation hc;
    
    SET_LOG_LEVEL(Info);
    std::cout << "begin demo10a" << std::endl;
    std::vector<bool> hasInfinity;
    for (int i = 17; i < cache.numberOfRoots() && i < 18; i++) {
        std::string file = "BS_L" + ToString(L)+ "M" + ToString(M)+ "r" + ToString(i) + ".txt";
        hc.setTraceFile(file);
        hc.hasInfinity = false;
        Vector startRoot = cache.getRoot(i);
        BetheSolution sol(startRoot, L, 0.0, totalPhase / hom.getSteps());
        hc.solve(hom, sol);
        hasInfinity.push_back(hc.hasInfinity);
        elem_t m1 = momentum(startRoot);
        elem_t m2 = momentum(sol.normalRoot());
        BOOST_TEST_INFO("i=" << i << ", m1/Pi=" << m1/(elem_t)Pi << ", m2/Pi=" << m2/(elem_t)Pi);
        BOOST_TEST(std::abs(floatMod(m2 - m1, 2.0 * Pi) - expectedMomentumDiff) < 1e-4);
        std::cout << "from " << startRoot << std::endl;
        std::cout << "to " << sol.normalRoot() << ", error=" << hom.evalTarget(sol) << std::endl;
    }
    std::cout << hasInfinity << std::endl;
    std::cout << "end demo10" << std::endl;    
}


BOOST_AUTO_TEST_CASE(Demo10, * utf::disabled())
{
    int L = 8;
    int M = 4;
    int steps = 3000;
    elem_t totalPhase = Pi;
    elem_t expectedMomentumDiff = .0;
    PoleFreeBetheHomotopy hom(L, M, totalPhase);
    hom.setSteps(steps);
    BetheRootCache cache(L, M);
    BetheHomotopyContinuation hc;
    
    // i = 5 four infinities, two pairs of conjugates
    // i = 6 three infinities, of which two are conjugates
    // i = 10 two identical infinities, two conjugates
    // i = 13 three infinities, of which two are conjugates
    std::cout << "begin demo10" << std::endl;
    std::vector<bool> hasInfinity;
    for (int i = 0; i < cache.numberOfRoots(); i++) {
        std::string file = "BS_L" + ToString(L)+ "M" + ToString(M)+ "r" + ToString(i) + ".txt";
        hc.setTraceFile(file);
        hc.hasInfinity = false;
        Vector startRoot = cache.getRoot(i);
        BetheSolution sol(startRoot, L, 0.0, totalPhase / steps);
        hc.solve(hom, sol);
        hasInfinity.push_back(hc.hasInfinity);
        elem_t m1 = momentum(startRoot);
        elem_t m2 = momentum(sol.normalRoot());
        BOOST_TEST_INFO("i=" << i << ", m1/Pi=" << m1/(elem_t)Pi << ", m2/Pi=" << m2/(elem_t)Pi);
        BOOST_TEST(std::abs(floatMod(m2 - m1, 2.0 * Pi) - expectedMomentumDiff) < 1e-4);
        std::cout << "from " << startRoot << std::endl;
        std::cout << "to " << sol.normalRoot() << ", error=" << hom.evalTarget(sol) << std::endl;
    }
    std::cout << hasInfinity << std::endl;
    std::cout << "end demo10" << std::endl;    
}

BOOST_AUTO_TEST_CASE(Demo11, * utf::disabled())
{
    int L = 12;
    int M = 6;
    int steps = 2000;
    elem_t totalPhase = Pi;
    elem_t expectedMomentumDiff = .0;
    PoleFreeBetheHomotopy hom(L, M, totalPhase);
    hom.setSteps(steps);
    BetheRootCache cache(L, M);
    BetheHomotopyContinuation hc;
    
    // i = 7 two finite conjugates
    // i = 8 two finite conjugates, two identical finite roots
    // i = 9 two identical finite, two singulars
    // i = 10 two pairs of finite conjugates
    std::vector<int> failed;
    std::cout << "begin demo11" << std::endl;
    for (int i = 0; i < cache.numberOfRoots(); i++) {
        std::string file = "BS_L" + ToString(L)+ "M" + ToString(M)+ "r" + ToString(i) + ".txt";
        hc.setTraceFile(file);
        Vector startRoot = cache.getRoot(i);
        BetheSolution sol(startRoot, L, 0.0, totalPhase / steps);
        hc.solve(hom, sol);
        elem_t m1 = momentum(startRoot);
        elem_t m2 = momentum(sol.normalRoot());
        BOOST_TEST_INFO("i=" << i << ", m1/Pi=" << m1/(elem_t)Pi << ", m2/Pi=" << m2/(elem_t)Pi);
        BOOST_TEST(std::abs(floatMod(m2 - m1, 2.0 * Pi) - expectedMomentumDiff) < 1e-4);
        if (!(std::abs(floatMod(m2 - m1, 2.0 * Pi) - expectedMomentumDiff) < 1e-4)) {
            failed.push_back(i);
        }
        std::cout << "from " << startRoot << std::endl;
        std::cout << "to " << sol.normalRoot() << ", error=" << hom.evalTarget(sol) << std::endl;
    }
    std::cout << "Failed cases: " << failed << std::endl;
    std::cout << "end demo11" << std::endl;    
}

BOOST_AUTO_TEST_SUITE_END()
