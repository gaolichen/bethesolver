#pragma once
#include "common.h"
#include "Solution.h"

class Function
{
protected:
    int _nEq;
public:
    int numberOfEquations() {
        return _nEq;
    }
    
    virtual Vector eval(const Solution &sol) = 0;
    
    // return a matrix m where m(i, j) = d f_i/d x_j
    virtual Matrix diff(const Solution &sol) = 0;
};

class SolvedEquation : public Function
{
public:
    virtual int numberOfRoots() = 0;
    virtual Vector getRoot(int rIndex) = 0;
};

class UnityEquation : public SolvedEquation
{
private:
    int _order;
public:
    UnityEquation(int nEq, int order) {
        _nEq = nEq;
        _order = order;
    }
    
    virtual int numberOfRoots() {
        return pow(_order, _nEq);
    }
    
    virtual Vector getRoot(int index) {
        Vector root(_nEq);
        elem_t phi = 2 * Pi / _order;
        for (int n = 0; n < root.size(); n++) {
            int ord = index % _order;
            index /= _order;
            root(n) = var_t(cos(ord * phi), sin(ord * phi));
        }
        return root;        
    }
    
    virtual Vector eval(const Solution &sol) {
        Vector res(sol.size());
        if (sol.size() != _nEq) {
            std::cout << "argument error: vars.size() = " << sol.size() << ", _nEq=" << _nEq << std::endl;
            return res;
        }
        
        for (int n = 0; n < sol.size(); n++) {
            res(n) = std::pow(sol.root(n), _order) - var_t(1.0);
        }
        
        return res;
    }
    
    virtual Matrix diff(const Solution &sol) {
        Matrix ret = Matrix::Zero(sol.size(), sol.size());
        for (int n = 0; n < ret.rows(); n++) {
            ret(n, n) = std::pow(sol.root(n), _order - 1) * (elem_t)_order;
        }
        
        return ret;
    }
};

class PolynomialFunction : public Function
{
private:
    std::vector<var_t> _coefs;
public:
    PolynomialFunction(std::vector<var_t>& coefs) {
        _coefs = coefs;
        _nEq = 1;
    }
    
    virtual Vector eval(const Solution &sol) {
        Vector ret = Vector::Zero(1);
        for (int i = _coefs.size() - 1; i >= 0; i--) {
            ret(0) = ret(0) * sol.root(0) + _coefs[i];
        }
        return ret;
    }
    
    virtual Matrix diff(const Solution &sol) {
        Matrix ret = Matrix::Zero(1, 1);
        for (int i = _coefs.size() - 1; i > 0; i--) {
            ret(0, 0) = ret(0, 0) * sol.root(0) + _coefs[i] * (elem_t)i;
        }
        return ret;
    }
};

