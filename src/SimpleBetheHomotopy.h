#pragma once
#include "common.h"
#include "BetheRootCache.h"
#include "Homotopy.h"

class SolvedBetheEquation : public SolvedEquation
{
private:
    int _L;
    int _M;
public:
    SolvedBetheEquation(int L, int M) {
        _nEq = M;
        _L = L;
        _M = M;
    }
        
    int chainLength() { return _L; }
    
    int magnonNumber() { return _M; }
        
    virtual Vector eval(Solution &sol);
    
    // return a matrix m where m(i, j) = d f_i/d x_j
    virtual Matrix diff(Solution &sol);
};

class SimpleBetheHomotopy : public Homotopy
{
private:
    elem_t _beta;
    SolvedBetheEquation *_start;
public:
    SimpleBetheHomotopy(elem_t beta, SolvedBetheEquation *start) {
        _beta = beta;
        _start = start;
    }
    
    virtual elem_t tBegin() {
        return 0.0;
    }
    
    virtual elem_t tEnd() {
        // when beta * L = Pi, the deformed equation become the undeformed one.
        return _beta - floor(_beta * _start->chainLength()/Pi) * Pi/_start->chainLength();
    }
    
    virtual elem_t tolerance() {
        return 1e-3;
    }
        
    virtual Vector eval(elem_t t, Solution &sol) {
        Vector ret = _start->eval(sol);
        var_t val = var_t(cos(2.0 * _start->chainLength() * t), sin(2.0 * _start->chainLength() * t));
        for (int i = 0; i < ret.size(); i++) {
            ret[i] -= val;
        }
        return ret;
    }
    
    virtual Matrix diff(elem_t t, Solution &sol) {
        return _start->diff(sol);
    }
    
    virtual Vector impliciteDiff(elem_t t, Solution &sol) {
        Vector rhs = _start->eval(sol);
        for (int i = 0; i < rhs.size(); i++) {
            rhs[i] *= var_t(0, 2.0 * _start->chainLength());
        }
        return _start->diff(sol).inverse() * rhs;
    }
};

