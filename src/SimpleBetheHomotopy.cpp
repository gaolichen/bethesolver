#include "SimpleBetheHomotopy.h"

Vector SolvedBetheEquation::eval(Solution &sol) {
    Vector res(_M);
    Vector &vars = sol.root();
    for (int j = 0; j < vars.size(); j++) {
        res[j] = std::pow(e2ip(vars[j]), _L);
        for (int k = 0; k < vars.size(); k++) {
            if (k == j) continue;
            res[j] *= sMatrix(vars[j], vars[k]);
        }
    }
    return res;
}

// return a matrix m where m(i, j) = d f_i/d x_j
Matrix SolvedBetheEquation::diff(Solution &sol) {
    Matrix ret(_M, _M);
    Vector &u = sol.root();
    Vector val = eval(sol);
    for (int j = 0; j < ret.rows(); j++) {
        ret(j, j) = var_t(0, -1.0 * _L)/(u[j] * u[j] + (elem_t)0.25) * val[j]; 
        for (int k = 0; k < ret.cols(); k++) {
            if (j == k) continue;
            ret(j, k) = val[j] * var_t(0, -2.0)/((u[j] - u[k])*(u[j] - u[k]) + (elem_t)1.0);
            ret(j, j) -= ret(j, k);
        }
    }
    return ret;
}
