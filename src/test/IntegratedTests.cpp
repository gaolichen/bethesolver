//Link to Boost
#define BOOST_TEST_DYN_LINK

//VERY IMPORTANT - include this last
//#include <boost/test/included/unit_test.hpp>
#include <boost/test/unit_test.hpp>

#include "test.h"
#include "../HomotopyContinuation.h"
#include "../BetheHomotopy.h"
#include "../BethePolynomialHomotopy.h"


// test suite
BOOST_FIXTURE_TEST_SUITE(Integrated_suite, SimpleTestFixture, * utf::label("IntegratedTests"))

BOOST_AUTO_TEST_CASE(SimpleHomotopy_test1)
{
    UnityEquation start(1, 2);
    std::vector<var_t> coefs({-20.0,1.0,1.0});
    PolynomialFunction target(coefs);
    SimpleHomotopy sh(&start, &target);
    sh.setSteps(100);
    SimpleHomotopyContinuation hc;
    
    // one root is 5 and the other one is -4
    bool found1 = false;
    bool found2 = false;
    
    for (int i = 0; i < start.numberOfRoots(); i++) {
        Vector st = start.getRoot(i);
        Solution sol(st);
        hc.solve(sh, sol);
        BOOST_TEST(sol.root().size() == 1);
        if (std::abs(sol.get(0) - (elem_t)4.0)) {
            found1 = true;
        }
        if (std::abs(sol.get(0) + (elem_t)5.0)) {
            found2 = true;
        }
    }
    BOOST_TEST(found1);
    BOOST_TEST(found2);
}

BOOST_AUTO_TEST_CASE(SimpleHomotopy_test2)
{
    UnityEquation start(1, 4);
    std::vector<var_t> coefs({-20.0,0,1.0,0,1.0});
    PolynomialFunction target(coefs);
    SimpleHomotopy sh(&start, &target);
    sh.setSteps(100);
    SimpleHomotopyContinuation hc;
    
    // there are found roots.
    bool found1 = false;
    bool found2 = false;
    bool found3 = false;
    bool found4 = false;
    
    for (int i = 0; i < start.numberOfRoots(); i++) {
        Vector startRoot = start.getRoot(i);
        Solution sol(startRoot);
        hc.solve(sh, sol);
        if (std::abs(sol.get(0) - (elem_t)2.0) < EPS) {
            found1 = true;
        }
        if (std::abs(sol.get(0) + (elem_t)2.0) < EPS) {
            found2 = true;
        }
        if (std::abs(sol.get(0) - var_t((elem_t).0, (elem_t)sqrt(5.0))) < 1e-6) {
            found3 = true;
        }
        if (std::abs(sol.get(0) + var_t((elem_t).0, (elem_t)sqrt(5.0))) < 1e-6) {
            found4 = true;
        }
    }
    BOOST_TEST(found1);
    BOOST_TEST(found2);
    BOOST_TEST(found3);
    BOOST_TEST(found4);
}

BOOST_AUTO_TEST_CASE(BethePolynomialHomotopy_L6M2)
{
    std::cout << "Running integrated test for L=6, M=2 cases" << std::endl;
    int L = 6;
    int M = 2;
    elem_t totalPhase = Pi;
    elem_t expectedMomentumDiff = floatMod(2 * M * totalPhase, 2 * Pi);
    BethePolynomialHomotopy hom(L, M, totalPhase);
    BetheRootCache cache(L, M);
    BetheHomotopyContinuation hc;
    
    SET_LOG_LEVEL(Error);
    for (int i = 0; i < cache.numberOfRoots(); i++) {
        LOG("processing root #" << i << ":\t" << cache.getRoot(i), Warning);
        std::string file = "BS_L" + ToString(L)+ "M" + ToString(M)+ "r" + ToString(i) + ".txt";
        hc.setTraceFile(file);

        Vector startRoot = cache.getRoot(i);
        InvertableSolution sol(startRoot, L);
        std::cout << "eps=" << sol.epsilon() << ", c=" << sol.c() << std::endl;
        hc.solve(hom, sol);
        std::cout <<"i=" << i << ", errors=" << hc.errors() << ", gaps=" << hc.gaps() << std::endl;
        BOOST_TEST_INFO("i=" << ", errors=" << hc.errors() << ", gaps=" << hc.gaps());
        BOOST_TEST(hc.gaps().max() < .5L);

        elem_t m1 = momentum(startRoot);
        elem_t m2 = momentum(sol.normalRoot());
        BOOST_TEST_INFO("i=" << i << ", m1/Pi=" << m1/(elem_t)Pi << ", m2/Pi=" << m2/(elem_t)Pi << ", expectedDiff/Pi=" << expectedMomentumDiff/(elem_t)Pi);
        BOOST_TEST(std::abs(floatMod(m2 - m1 - expectedMomentumDiff, 2.0 * Pi)) < 1e-4);
    }
    SET_LOG_LEVEL(Error);
}

BOOST_AUTO_TEST_CASE(BethePolynomialHomotopy_L6M3)
{
    std::cout << "Running integrated test for L=6, M=2 cases" << std::endl;
    int L = 6;
    int M = 3;
    elem_t totalPhase = Pi;
    elem_t expectedMomentumDiff = floatMod(2 * M * totalPhase, 2 * Pi);
    BethePolynomialHomotopy hom(L, M, totalPhase);
    BetheRootCache cache(L, M);
    BetheHomotopyContinuation hc;
    
    SET_LOG_LEVEL(Error);
    for (int i = 0; i < cache.numberOfRoots(); i++) {
        LOG("processing root #" << i << ":\t" << cache.getRoot(i), Warning);
        Vector startRoot = cache.getRoot(i);
        InvertableSolution sol(startRoot, L);
        std::cout << "eps=" << sol.epsilon() << ", c=" << sol.c() << std::endl;
        hc.solve(hom, sol);
        std::cout <<"i=" << i << ", errors=" << hc.errors() << ", gaps=" << hc.gaps() << std::endl;
        BOOST_TEST_INFO("i=" << ", errors=" << hc.errors() << ", gaps=" << hc.gaps());
        BOOST_TEST(hc.gaps().max() < .5L);

        elem_t m1 = momentum(startRoot);
        elem_t m2 = momentum(sol.normalRoot());
        BOOST_TEST_INFO("i=" << i << ", m1/Pi=" << m1/(elem_t)Pi << ", m2/Pi=" << m2/(elem_t)Pi << ", expectedDiff/Pi=" << expectedMomentumDiff/(elem_t)Pi);
        BOOST_TEST(std::abs(floatMod(m2 - m1 - expectedMomentumDiff, 2.0 * Pi)) < 1e-4);
    }
    SET_LOG_LEVEL(Error);
}

BOOST_AUTO_TEST_CASE(BethePolynomialHomotopy_L7M2)
{
    int L = 7;
    int M = 2;
    std::cout << "Running integrated test for L=" << L << ", M=" << M << " cases" << std::endl;
    elem_t totalPhase = Pi;
    elem_t expectedMomentumDiff = floatMod(2 * M * totalPhase, 2 * Pi);
    BethePolynomialHomotopy hom(L, M, totalPhase);
    BetheRootCache cache(L, M);
    BetheHomotopyContinuation hc;
    
    SET_LOG_LEVEL(Error);
    for (int i = 0; i < cache.numberOfRoots(); i++) {
        LOG("processing root #" << i << ":\t" << cache.getRoot(i), Warning);
        Vector startRoot = cache.getRoot(i);
        InvertableSolution sol(startRoot, L);
        std::cout << "eps=" << sol.epsilon() << ", c=" << sol.c() << std::endl;
        try {
            hc.solve(hom, sol);
        } catch (HomotopyContinuationException &ex) {
            std::cout << ex.what() << std::endl;
//            throw ex;
        }
        std::cout <<"i=" << i << ", errors=" << hc.errors() << ", gaps=" << hc.gaps() << std::endl;
        BOOST_TEST_INFO("i=" << ", errors=" << hc.errors() << ", gaps=" << hc.gaps());
        BOOST_TEST(hc.gaps().max() < .5L);

        elem_t m1 = momentum(startRoot);
        elem_t m2 = momentum(sol.normalRoot());
        BOOST_TEST_INFO("i=" << i << ", m1/Pi=" << m1/(elem_t)Pi << ", m2/Pi=" << m2/(elem_t)Pi << ", expectedDiff/Pi=" << expectedMomentumDiff/(elem_t)Pi);
        BOOST_TEST(std::abs(floatMod(m2 - m1 - expectedMomentumDiff, 2.0 * Pi)) < 1e-4);
    }
    SET_LOG_LEVEL(Error);
}

BOOST_AUTO_TEST_CASE(BethePolynomialHomotopy_L7M3)
{
    int L = 7;
    int M = 3;
    std::cout << "Running integrated test for L=" << L << ", M=" << M << " cases" << std::endl;
    elem_t totalPhase = Pi;
    elem_t expectedMomentumDiff = floatMod(2 * M * totalPhase, 2 * Pi);
    BethePolynomialHomotopy hom(L, M, totalPhase);
    BetheRootCache cache(L, M);
    BetheHomotopyContinuation hc;
    
    SET_LOG_LEVEL(Error);
    for (int i = 0; i < cache.numberOfRoots(); i++) {
        LOG("processing root #" << i << ":\t" << cache.getRoot(i), Warning);
        std::string file = "BS_L" + ToString(L)+ "M" + ToString(M)+ "r" + ToString(i) + ".txt";
//        hc.setTraceFile(file);

        Vector startRoot = cache.getRoot(i);
        InvertableSolution sol(startRoot, L);
        if (sol.indexOfSingularRoot1() >= 0) {
            std::cout << "eps=" << sol.epsilon() << ", c=" << sol.c() << std::endl;
        }
        try {
            hc.solve(hom, sol);
        } catch (HomotopyContinuationException &ex) {
            std::cout << ex.what() << std::endl;
            throw ex;
        }
        
        std::cout <<"i=" << i << std::endl << "errors=" << hc.errors() << std::endl << "gaps=" << hc.gaps() << std::endl;
        BOOST_TEST_INFO("i=" << ", errors=" << hc.errors() << ", gaps=" << hc.gaps());
        BOOST_TEST(hc.gaps().max() < .5L);
        elem_t m1 = momentum(startRoot);
        elem_t m2 = momentum(sol.normalRoot());
        BOOST_TEST_INFO("i=" << i << ", m1/Pi=" << m1/(elem_t)Pi << ", m2/Pi=" << m2/(elem_t)Pi << ", expectedDiff/Pi=" << expectedMomentumDiff/(elem_t)Pi);
        BOOST_TEST(std::abs(floatMod(m2 - m1 - expectedMomentumDiff, 2.0 * Pi)) < 1e-4);
    }
    SET_LOG_LEVEL(Error);
}

BOOST_AUTO_TEST_CASE(BethePolynomialHomotopy_L8M2)
{
    int L = 8;
    int M = 2;
    std::cout << "Running integrated test for L=" << L << ", M=" << M << " cases" << std::endl;
    elem_t totalPhase = Pi;
    elem_t expectedMomentumDiff = floatMod(2 * M * totalPhase, 2 * Pi);
    BethePolynomialHomotopy hom(L, M, totalPhase);
    BetheRootCache cache(L, M);
    BetheHomotopyContinuation hc;
    
    SET_LOG_LEVEL(Error);
    for (int i = 0; i < cache.numberOfRoots(); i++) {
        LOG("processing root #" << i << ":\t" << cache.getRoot(i), Warning);
        std::string file = "BS_L" + ToString(L)+ "M" + ToString(M)+ "r" + ToString(i) + ".txt";
//        hc.setTraceFile(file);
        Vector startRoot = cache.getRoot(i);
        InvertableSolution sol(startRoot, L);
        std::cout << "eps=" << sol.epsilon() << ", c=" << sol.c() << std::endl;
        hc.solve(hom, sol);
        std::cout <<"i=" << i << ", errors=" << hc.errors() << ", gaps=" << hc.gaps() << std::endl;
        BOOST_TEST_INFO("i=" << ", errors=" << hc.errors() << ", gaps=" << hc.gaps());
        BOOST_TEST(hc.gaps().max() < .5L);

        elem_t m1 = momentum(startRoot);
        elem_t m2 = momentum(sol.normalRoot());
        BOOST_TEST_INFO("i=" << i << ", m1/Pi=" << m1/(elem_t)Pi << ", m2/Pi=" << m2/(elem_t)Pi << ", expectedDiff/Pi=" << expectedMomentumDiff/(elem_t)Pi);
        BOOST_TEST(std::abs(floatMod(m2 - m1 - expectedMomentumDiff, 2.0 * Pi)) < 1e-4);
    }
    SET_LOG_LEVEL(Error);
}

BOOST_AUTO_TEST_CASE(BethePolynomialHomotopy_L8M3)
{
    int L = 8;
    int M = 3;
    std::cout << "Running integrated test for L=" << L << ", M=" << M << " cases" << std::endl;
    elem_t totalPhase = Pi;
    elem_t expectedMomentumDiff = floatMod(2 * M * totalPhase, 2 * Pi);
    BethePolynomialHomotopy hom(L, M, totalPhase);
    BetheRootCache cache(L, M);
    BetheHomotopyContinuation hc;
    
    SET_LOG_LEVEL(Error);
    for (int i = 0; i < cache.numberOfRoots(); i++) {
        LOG("processing root #" << i << ":\t" << cache.getRoot(i), Warning);
        Vector startRoot = cache.getRoot(i);
        InvertableSolution sol(startRoot, L);
        std::cout << "eps=" << sol.epsilon() << ", c=" << sol.c() << std::endl;
        hc.solve(hom, sol);
        std::cout <<"i=" << i << ", errors=" << hc.errors() << ", gaps=" << hc.gaps() << std::endl;
        BOOST_TEST_INFO("i=" << ", errors=" << hc.errors() << ", gaps=" << hc.gaps());
        BOOST_TEST(hc.gaps().max() < .5L);

        elem_t m1 = momentum(startRoot);
        elem_t m2 = momentum(sol.normalRoot());
        BOOST_TEST_INFO("i=" << i << ", m1/Pi=" << m1/(elem_t)Pi << ", m2/Pi=" << m2/(elem_t)Pi << ", expectedDiff/Pi=" << expectedMomentumDiff/(elem_t)Pi);
        BOOST_TEST(std::abs(floatMod(m2 - m1 - expectedMomentumDiff, 2.0 * Pi)) < 1e-4);
    }
    SET_LOG_LEVEL(Error);
}

BOOST_AUTO_TEST_CASE(BethePolynomialHomotopy_L8M4)
{
    int L = 8;
    int M = 4;
    std::cout << "Running integrated test for L=" << L << ", M=" << M << " cases" << std::endl;
    elem_t totalPhase = 2 * Pi;
    elem_t expectedMomentumDiff = floatMod(2 * M * totalPhase, 2 * Pi);
    BethePolynomialHomotopy hom(L, M, totalPhase);
    BetheRootCache cache(L, M);
    BetheHomotopyContinuation hc;
    hc.setMaxNumberOfNewtonRaphonEachStep(50);
    hc.setMaxDepth(14);
    
    SET_LOG_LEVEL(Error);
    for (int i = 0; i < cache.numberOfRoots(); i++) {
        LOG("processing root #" << i << ":\t" << cache.getRoot(i), Warning);
        std::string file = "BS_L" + ToString(L)+ "M" + ToString(M)+ "r" + ToString(i) + ".txt";
//        hc.setTraceFile(file);

        Vector startRoot = cache.getRoot(i);
        InvertableSolution sol(startRoot, L);
        std::cout << "eps=" << sol.epsilon() << ", c=" << sol.c() << std::endl;
        hc.solve(hom, sol);
        std::cout <<"i=" << i << ", errors=" << hc.errors() << std::endl << "gaps=" << hc.gaps() << std::endl;
        BOOST_TEST_INFO("i=" << ", errors=" << hc.errors() << ", gaps=" << hc.gaps());
        BOOST_TEST(hc.gaps().max() < .5L);

        elem_t m1 = momentum(startRoot);
        elem_t m2 = momentum(sol.normalRoot());
        BOOST_TEST_INFO("i=" << i << ", m1/Pi=" << m1/(elem_t)Pi << ", m2/Pi=" << m2/(elem_t)Pi << ", expectedDiff/Pi=" << expectedMomentumDiff/(elem_t)Pi);
        BOOST_TEST(std::abs(floatMod(m2 - m1 - expectedMomentumDiff, 2.0 * Pi)) < 1e-4);
    }
    SET_LOG_LEVEL(Error);
}

BOOST_AUTO_TEST_CASE(BethePolynomialHomotopy_L10M4)
{
    int L = 10;
    int M = 4;
    std::cout << "Running integrated test for L=" << L << ", M=" << M << " cases" << std::endl;
    elem_t totalPhase = Pi;
    elem_t expectedMomentumDiff = floatMod(2 * M * totalPhase, 2 * Pi);
    BethePolynomialHomotopy hom(L, M, totalPhase);
    BetheRootCache cache(L, M);
    BetheHomotopyContinuation hc;
    hc.setMaxNumberOfNewtonRaphonEachStep(50);
    
    SET_LOG_LEVEL(Error);
    for (int i = 0; i < cache.numberOfRoots(); i++) {
        LOG("processing root #" << i << ":\t" << cache.getRoot(i), Warning);
        std::string file = "BS_L" + ToString(L)+ "M" + ToString(M)+ "r" + ToString(i) + ".txt";
//        hc.setTraceFile(file);
        Vector startRoot = cache.getRoot(i);
        InvertableSolution sol(startRoot, L);
        if (sol.indexOfSingularRoot1() >= 0) { 
            std::cout << "eps=" << sol.epsilon() << ", c=" << sol.c() << std::endl;
        }
        
        hc.solve(hom, sol);
        std::cout <<"i=" << i << std::endl << "errors=" << hc.errors() << std::endl << "gaps=" << hc.gaps() << std::endl;
        BOOST_TEST_INFO("i=" << ", errors=" << hc.errors() << ", gaps=" << hc.gaps());
        BOOST_TEST(hc.gaps().max() < .5L);
        elem_t m1 = momentum(startRoot);
        elem_t m2 = momentum(sol.normalRoot());
        BOOST_TEST_INFO("i=" << i << ", m1/Pi=" << m1/(elem_t)Pi << ", m2/Pi=" << m2/(elem_t)Pi << ", expectedDiff/Pi=" << expectedMomentumDiff/(elem_t)Pi << ",normalRoot=" << sol.root());
        BOOST_TEST(std::abs(floatMod(m2 - m1 - expectedMomentumDiff, 2.0 * Pi)) < 1e-4);
    }
    SET_LOG_LEVEL(Error);
}

BOOST_AUTO_TEST_CASE(BethePolynomialHomotopy_L10M5)
{
    int L = 10;
    int M = 5;
    std::cout << "Running integrated test for L=" << L << ", M=" << M << " cases" << std::endl;
    elem_t totalPhase = Pi;
    elem_t expectedMomentumDiff = floatMod(2 * M * totalPhase, 2 * Pi);
    BethePolynomialHomotopy hom(L, M, totalPhase);
    BetheRootCache cache(L, M);
    BetheHomotopyContinuation hc;
    hc.setMaxNumberOfNewtonRaphonEachStep(50);
    
    SET_LOG_LEVEL(Error);
    for (int i = 0; i < cache.numberOfRoots(); i++) {
//        if (i == 15) continue;
        LOG("processing root #" << i << ":\t" << cache.getRoot(i), Warning);
        std::string file = "BS_L" + ToString(L)+ "M" + ToString(M)+ "r" + ToString(i) + ".txt";
//        hc.setTraceFile(file);
        Vector startRoot = cache.getRoot(i);
        InvertableSolution sol(startRoot, L);
        if (sol.indexOfSingularRoot1() >= 0) { 
            std::cout << "eps=" << sol.epsilon() << ", c=" << sol.c() << std::endl;
        }
        
        hc.solve(hom, sol);
        std::cout <<"i=" << i << std::endl << "errors=" << hc.errors() << std::endl << "gaps=" << hc.gaps() << std::endl;
        BOOST_TEST_INFO("i=" << ", errors=" << hc.errors() << ", gaps=" << hc.gaps());
        BOOST_TEST(hc.gaps().max() < .5L);
        elem_t m1 = momentum(startRoot);
        elem_t m2 = momentum(sol.normalRoot());
        BOOST_TEST_INFO("i=" << i << ", m1/Pi=" << m1/(elem_t)Pi << ", m2/Pi=" << m2/(elem_t)Pi << ", expectedDiff/Pi=" << expectedMomentumDiff/(elem_t)Pi << ",normalRoot=" << sol.root());
        BOOST_TEST(std::abs(floatMod(m2 - m1 - expectedMomentumDiff, 2.0 * Pi)) < 1e-4);
    }
    SET_LOG_LEVEL(Error);
}


BOOST_AUTO_TEST_CASE(PoleFreeBetheHomotopy_L6M2)
{
    std::cout << "Running integrated test for L=6, M=2 cases" << std::endl;
    int L = 6;
    int M = 2;
    elem_t totalPhase = Pi;
    elem_t expectedMomentumDiff = floatMod(2 * M * totalPhase, 2 * Pi);
    PoleFreeBetheHomotopy hom(L, M, totalPhase);
    BetheRootCache cache(L, M);
    BetheHomotopyContinuation hc;
    
    SET_LOG_LEVEL(Error);
    for (int i = 0; i < cache.numberOfRoots(); i++) {
        LOG("processing root #" << i << ":\t" << cache.getRoot(i), Warning);
        Vector startRoot = cache.getRoot(i);
        BetheSolution sol(startRoot, L, 0.0, totalPhase / hom.getSteps());
        std::cout << "eps=" << sol.epsilon() << ", c=" << sol.c() << std::endl;
        hc.solve(hom, sol);
        elem_t m1 = momentum(startRoot);
        elem_t m2 = momentum(sol.normalRoot());
        BOOST_TEST_INFO("i=" << i << ", m1/Pi=" << m1/(elem_t)Pi << ", m2/Pi=" << m2/(elem_t)Pi << ", expectedDiff/Pi=" << expectedMomentumDiff/(elem_t)Pi);
        BOOST_TEST(std::abs(floatMod(m2 - m1 - expectedMomentumDiff, 2.0 * Pi)) < 1e-4);
    }
    SET_LOG_LEVEL(Error);
}

BOOST_AUTO_TEST_CASE(PoleFreeBetheHomotopy_L6M3)
{
    std::cout << "Running integrated test for L=6, M=3 cases" << std::endl;

    int L = 6;
    int M = 3;
    elem_t totalPhase = Pi;
    elem_t expectedMomentumDiff = floatMod(2 * M * totalPhase, 2 * Pi);
    PoleFreeBetheHomotopy hom(L, M, totalPhase);
    BetheRootCache cache(L, M);
    BetheHomotopyContinuation hc;
    SET_LOG_LEVEL(Error);
    
    for (int i = 0; i < cache.numberOfRoots(); i++) {
        Vector startRoot = cache.getRoot(i);
        BetheSolution sol(startRoot, L, 0.0, totalPhase / hom.getSteps());
        hc.solve(hom, sol);
        elem_t m1 = momentum(startRoot);
        elem_t m2 = momentum(sol.normalRoot());
        BOOST_TEST_INFO("i=" << i << ", m1/Pi=" << m1/(elem_t)Pi << ", m2/Pi=" << m2/(elem_t)Pi << ", expectedDiff/Pi=" << expectedMomentumDiff/(elem_t)Pi);
        BOOST_TEST(std::abs(floatMod(m2 - m1 - expectedMomentumDiff, 2.0 * Pi)) < 1e-4);
    }
}

BOOST_AUTO_TEST_CASE(PoleFreeBetheHomotopy_L7M2)
{
    std::cout << "Running integrated test for L=7, M=2 cases" << std::endl;
    int L = 7;
    int M = 2;
    elem_t totalPhase = random(1.0, 2.0) * Pi;
//    elem_t totalPhase = 1.27095 * Pi;
//    elem_t totalPhase = 1.81386 * Pi;
    elem_t expectedMomentumDiff = floatMod(2 * M * totalPhase, 2 * Pi);
    PoleFreeBetheHomotopy hom(L, M, totalPhase);
    BetheRootCache cache(L, M);
    BetheHomotopyContinuation hc;
    SET_LOG_LEVEL(Error);
    
    // i = 1: at two identical infinity roots
    // i = 3 and 8: two conjugate finite roots.
    for (int i = 0; i < cache.numberOfRoots(); i++) {
        Vector startRoot = cache.getRoot(i);
        BetheSolution sol(startRoot, L, 0.0, totalPhase / hom.getSteps());
        hc.solve(hom, sol);
        elem_t m1 = momentum(startRoot);
        elem_t m2 = momentum(sol.normalRoot());
        BOOST_TEST_INFO("totalPhass/Pi=" << totalPhase/Pi << ", i=" << i << ", m1/Pi=" << m1/(elem_t)Pi << ", m2/Pi=" << m2/(elem_t)Pi << ", expectedDiff/Pi=" << expectedMomentumDiff/(elem_t)Pi);
        BOOST_TEST(std::abs(floatMod(m2 - m1 - expectedMomentumDiff, 2.0 * Pi)) < 1e-4);
    }
}

BOOST_AUTO_TEST_CASE(PoleFreeBetheHomotopy_L7M3)
{
    std::cout << "Running integrated test for L=7, M=3 cases" << std::endl;
    int L = 7;
    int M = 3;
    elem_t totalPhase = 1.3453 * Pi;
    elem_t expectedMomentumDiff = floatMod(2 * M * totalPhase, 2 * Pi);
    PoleFreeBetheHomotopy hom(L, M, totalPhase);
    BetheRootCache cache(L, M);
    BetheHomotopyContinuation hc;
    SET_LOG_LEVEL(Error);
    
    // i = 1: at two identical infinity roots
    // i = 3 and 8: two conjugate finite roots.
    for (int i = 0; i < cache.numberOfRoots(); i++) {
        Vector startRoot = cache.getRoot(i);
        BetheSolution sol(startRoot, L, 0.0, totalPhase / hom.getSteps());
        hc.solve(hom, sol);
        elem_t m1 = momentum(startRoot);
        elem_t m2 = momentum(sol.normalRoot());
        BOOST_TEST_INFO("i=" << i << ", m1/Pi=" << m1/(elem_t)Pi << ", m2/Pi=" << m2/(elem_t)Pi << ", expectedDiff/Pi=" << expectedMomentumDiff/(elem_t)Pi);
        BOOST_TEST(std::abs(floatMod(m2 - m1 - expectedMomentumDiff, 2.0 * Pi)) < 1e-4);
    }
}

BOOST_AUTO_TEST_CASE(PoleFreeBetheHomotopy_L8M3)
{
    std::cout << "Running integrated test for L=8, M=3 cases" << std::endl;
    int L = 8;
    int M = 3;
//    elem_t totalPhase = random(1.0, 2.0) * Pi;
    elem_t totalPhase = 3.66417L;
    elem_t expectedMomentumDiff = floatMod(2 * M * totalPhase, 2 * Pi);
    PoleFreeBetheHomotopy hom(L, M, totalPhase);
    BetheRootCache cache(L, M);
    BetheHomotopyContinuation hc;
    SET_LOG_LEVEL(Error);
    
    // i = 1: at two identical infinity roots
    // i = 3 and 8: two conjugate finite roots.
    for (int i = 0; i < cache.numberOfRoots(); i++) {
//        std::string file = "BS_L" + ToString(L)+ "M" + ToString(M)+ "r" + ToString(i) + ".txt";
//        hc.setTraceFile(file);
        std::cout << "i=" << i << std::endl;

        Vector startRoot = cache.getRoot(i);
        BetheSolution sol(startRoot, L, 0.0, totalPhase / hom.getSteps());
        hc.solve(hom, sol);
        elem_t m1 = momentum(startRoot);
        elem_t m2 = momentum(sol.normalRoot());
        BOOST_TEST_INFO("i=" << i << ", m1/Pi=" << m1/(elem_t)Pi << ", m2/Pi=" << m2/(elem_t)Pi << ", expectedDiff/Pi=" << expectedMomentumDiff/(elem_t)Pi);
        BOOST_TEST(std::abs(floatMod(m2 - m1 - expectedMomentumDiff, 2.0 * Pi)) < 1e-4);
    }
}

BOOST_AUTO_TEST_SUITE_END()
