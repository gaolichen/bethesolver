//Link to Boost
#define BOOST_TEST_DYN_LINK

//Define our Module name (prints at testing)
#define BOOST_TEST_MODULE BetheSolver UnitTests

//VERY IMPORTANT - include this last
#include <boost/test/included/unit_test.hpp>
//#include <boost/test/unit_test.hpp>

#include "test.h"
#include "../HomotopyContinuation.h"

bool isInteger(elem_t n) {
    return std::abs(floor(n + eps) - n) < eps;
}

// test suite
BOOST_FIXTURE_TEST_SUITE(UnityEquation_suite, SimpleTestFixture, * utf::label("UnityEquation"))

BOOST_DATA_TEST_CASE(numberOfRoots_test, bdata::random(1, 4) ^ bdata::random(1, 7) ^ bdata::xrange(20), nEq, order, index)
{
    UnityEquation uss(nEq, order);
    int expected = 1;
    for (int i = 0; i < nEq; i++) {
        expected *= order;
    }
    
    BOOST_TEST(uss.numberOfRoots() == expected);
}

BOOST_DATA_TEST_CASE(getRoot_test, bdata::random(1, 4) ^ bdata::random(1, 7) ^ bdata::xrange(20), nEq, order, index)
{
    UnityEquation uss(nEq, order);
    std::vector<int> phases;
    for (int i = 0; i < uss.numberOfRoots(); i++) {
        Vector root = uss.getRoot(i);
        BOOST_TEST(root.size() == nEq);
        
        int index = 0;
        for (int j = 0; j < root.size(); j++) {
            BOOST_TEST(std::abs(std::abs(root[j]) - (elem_t)1.0) < eps);
            elem_t phi = std::arg(root[j]);
            BOOST_TEST(isInteger(phi/(2 * Pi/order)));
            int m = (int)floor(phi/(2 * Pi/order) + eps);
            if (m < 0) m += order;
            BOOST_TEST(m < order);
            index = index * order + m;
        }
        phases.push_back(index);
    }
    std::sort(phases.begin(), phases.end());
    std::vector<int>::iterator it = std::unique(phases.begin(), phases.end());
    phases.resize(std::distance(phases.begin(), it));
    BOOST_TEST(uss.numberOfRoots() == phases.size());
}

BOOST_DATA_TEST_CASE(eval_test, bdata::random(1, 4) ^ bdata::random(1, 7) ^ bdata::xrange(20), nEq, order, index)
{
    UnityEquation uss(nEq, order);
    Vector root = uss.getRoot(random64(uss.numberOfRoots()));
    Solution sol(root);
    Vector res = uss.eval(sol);
    for (int i = 0; i < res.size(); i++) {
        BOOST_TEST(abs(res[i]) < eps);
    }
}

BOOST_DATA_TEST_CASE(diff_test, bdata::random(1, 4) ^ bdata::random(1, 7) ^ bdata::xrange(20), nEq, order, index)
{
    UnityEquation uss(nEq, order);
    Vector root = uss.getRoot(random64(uss.numberOfRoots()));
    Solution sol(root);
    Matrix res = uss.diff(sol);
    for (int i = 0; i < res.rows(); i++) {
        for (int j = 0; j < res.cols(); j++) {
            if (i != j) {
                BOOST_TEST(std::abs(res(i, j)) < eps);
            } else {
                BOOST_TEST(abs(std::abs(res(i, j) * root[i]) - order) < eps);
            }
        }
    }
}


BOOST_AUTO_TEST_SUITE_END()

class QuardraticEquations : public Function
{
private:
    std::vector<var_t> _coefs1;
    std::vector<var_t> _coefs2;
public:
    QuardraticEquations() {
        _nEq = 2;
    }
    
    void setCoefs(std::vector<var_t> &coefs1, std::vector<var_t> &coefs2) {
        _coefs1 = coefs1;
        _coefs2 = coefs2;
    }
    
    static var_t eval(const Solution &sol, std::vector<var_t> coefs) {
        const Vector &vars = sol.root();
        return coefs[0] + coefs[1] * vars[0] + coefs[2] * vars[1] +
        coefs[3] * vars[0] * vars[0] + coefs[4] * vars[0] * vars[1] +
        coefs[5] * vars[1] * vars[1];
    }
    
    virtual Vector eval(const Solution &sol) {
        Vector ret(2);
        ret[0] = eval(sol, _coefs1);
        ret[1] = eval(sol, _coefs2);
        return ret;
    }
    
    virtual Matrix diff(const Solution &sol) {
        Matrix ret(2, 2);
        const Vector &vars = sol.root();
        ret(0, 0) = _coefs1[1] + _coefs1[3] * vars[0] * (elem_t)2.0 + _coefs1[4] * vars[1];
        ret(0, 1) = _coefs1[2] + _coefs1[4] * vars[0] + _coefs1[5] * vars[1] * (elem_t)2.0;
        ret(1, 0) = _coefs2[1] + _coefs2[3] * vars[0] * (elem_t)2.0 + _coefs2[4] * vars[1];
        ret(1, 1) = _coefs2[2] + _coefs2[4] * vars[0] + _coefs2[5] * vars[1] * (elem_t)2.0;
        return ret;
    }
};

// test suite
BOOST_FIXTURE_TEST_SUITE(MultiVariables_suite, SimpleTestFixture, * utf::label("MultiVariables"))


BOOST_AUTO_TEST_CASE(Demo1)
{
    UnityEquation start(2, 2);
    QuardraticEquations target;
    // construct equations with one solution (x=3,y=-2).
    // first equation is -12+3x-4y+x^2+xy-2y^2==0
    std::vector<var_t> coefs1({(elem_t)-12.0, (elem_t)3.0, (elem_t)-4.0, (elem_t)1.0, (elem_t)1.0, (elem_t)-2.0});
    // second equation is -32-2x+y+2x^2-3xy+y^2==0
    std::vector<var_t> coefs2({(elem_t)-32.0, (elem_t)-2.0, (elem_t)1.0, (elem_t)2.0, (elem_t)-3.0, (elem_t)1.0});
    target.setCoefs(coefs1, coefs2);
    SimpleHomotopy sh(&start, &target); 
    sh.setSteps(1000);   
    SimpleHomotopyContinuation hc;
    
    // there should be 3 solutions
    //.. todo
    
    Vector expected(2);
    expected << var_t(3.0), var_t(-2.0);
    bool found = false;
    for (int i = 0; i < start.numberOfRoots(); i++) {
        Vector startRoot = start.getRoot(i);
        Solution sol(startRoot);
        hc.solve(sh, sol);
        BOOST_TEST_INFO("root=" << sol.root() << " norm=" << sh.evalTarget(sol).norm());
        BOOST_TEST(sh.evalTarget(sol).norm() == 0, tt::tolerance(eps));
        if ((sol.root() - expected).norm() < eps) {
            found = true;
        }
//        std::cout << "solution " << i << ": " << target.getRoot(i) << std::endl;
//        std::cout << "verify: " << target.eval(target.getRoot(i)) << std::endl;
    }
    
    BOOST_TEST(found);
}

BOOST_AUTO_TEST_SUITE_END()

