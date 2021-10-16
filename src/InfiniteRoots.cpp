#include "InfiniteRoots.h"
#include "HomotopyContinuation.h"
#include <boost/math/tools/polynomial.hpp>

using namespace boost::math;
using namespace boost::math::tools; // for polynomial

InfiniteRoots::InfiniteRoots(int L, int M) {
    _L = L;
    _M = M;
    _roots.resize(M);
}

std::vector<var_t>& InfiniteRoots::getRoot(int m)
{
    if (_roots[m-1].size() == 0) {
        computeRoots(m);
    }
    
    return _roots[m-1];
}

std::vector<var_t> InfiniteRoots::equationSolver(std::vector<var_t> &coefs, std::vector<var_t> &out) {
    if (coefs.size() == 2) {
        return { -coefs[0]/coefs[1] };
    }
    if (coefs.size() <= 1) {
        return {};
    }
    
    PolynomialFunction target(coefs);
    UnityEquation start(1, coefs.size() - 1);
    
    SimpleHomotopy hom(&start, &target);
    hom.setSteps(coefs.size() * 300);
    hom.setRandSeed(coefs.size());
    SimpleHomotopyContinuation hc;
    
    std::vector<var_t> ret;
    
    for (int i = 0; i < start.numberOfRoots(); i++) {
        Vector root = start.getRoot(i);
        Solution sol(root);
        hc.solve(hom, sol);
        bool found = false;
        for (int j = 0; j < ret.size(); j++) {
            if (std::abs(ret[j] - sol.get(0)) < EPS) {
                found = true;
                break;
            }
        }
        if (!found) {
            ret.push_back(sol.get(0));
        }
    }
    
    // it produced some duplicated roots
    if (ret.size() < start.numberOfRoots()) {
        polynomial<var_t> eq(coefs.begin(), coefs.end());
        polynomial<var_t> eq2 {{var_t(1.0L, 0.0L)}};
        for (int i = 0; i < ret.size(); i++) {
            polynomial<var_t> solved {{-ret[i], var_t(1.0L, 0)}}; 
            eq2 *= solved;
        }
        polynomial<var_t> eq3 = eq / eq2;
        out = eq3.data();
    }
    return ret;
}

void InfiniteRoots::computeRoots(int m) {
    std::vector<var_t> coefs(m + 1);
    elem_t r = 1.0L;
    int R = _L - 2 * _M + 2 * m;
    for (int i = 0; i <= m; i++) {
        coefs[i] = var_t(r, 0.0L);
        r *= -(m-i) * (2 * _L) / (elem_t)((R - i) *  (i+1));
    }
    reverse(coefs.begin(), coefs.end());
    while (true) {
        std::vector<var_t> remain;
        std::vector<var_t> res = equationSolver(coefs, remain);
        _roots[m-1].insert(_roots[m-1].end(), res.begin(), res.end());
//        std::cout << "roots=" << _roots[m-1] << " remain=" << remain << std::endl;
        if (remain.size() <= 1) break;
        coefs = remain;
    }
}
