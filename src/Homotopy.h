#pragma once
#include "common.h"
#include "Function.h"

class Homotopy
{
protected:
    Vector _startRoot;
    int _steps = 1000;
public:        
    virtual elem_t tBegin() = 0;
    virtual elem_t tEnd() = 0;
    
    virtual elem_t tolerance() {
        return 1e-8;
    }
    
    int getSteps() {
        return _steps;
    }
    
    void setSteps(int steps) {
        _steps = steps;
    }

    virtual void setStartRoot(const Vector& root) {
        _startRoot = root;
    }
    
    virtual elem_t error(elem_t t, const Solution &sol) {
        return eval(t, sol).norm();
    }

    Vector evalTarget(const Solution &sol) {
        return eval(tEnd(), sol);
    }
    
    virtual bool acceptCorrection(elem_t t, const Solution &sol, const Vector& correction) {
        return true;
    }

    virtual Vector eval(elem_t t, const Solution &sol) = 0;
    virtual Matrix diff(elem_t t, const Solution &sol) = 0;
    virtual Vector impliciteDiff(elem_t t, const Solution &sol) = 0;
};

class SimpleHomotopy : public Homotopy
{
private:
    var_t _gamma;
    SolvedEquation* _start;
    Function* _target;
public:
    SimpleHomotopy(SolvedEquation* start, Function* target) {
        _gamma = randomComplex(0.5, 1.0);
//        std::cout << "gamma=" << _gamma << std::endl;
        this->_start = start;
        this->_target = target;
    }
    
    virtual elem_t tBegin() {
        return 1.0;
    }
    
    virtual elem_t tEnd() {
        return 0.0;
    }
    
    void setRandSeed(elem_t r) {
        _gamma = randomComplex(0.75 * r, 1.25 * r);
    }

    virtual Vector eval(elem_t t, const Solution &sol);
    virtual Matrix diff(elem_t t, const Solution &sol);
    virtual Vector impliciteDiff(elem_t t, const Solution &sol);
};

