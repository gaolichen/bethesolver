#pragma once
#include "common.h"
#include "Statistics.h"
#include "Homotopy.h"
#include <fstream>
#include <exception>

/**
 * header file to define Homotopy Continuation algorithm
 */

class HomotopyContinuation
{
private:
    elem_t _accuracyLow = 1e-8L;
    elem_t _accuracyHigh = 1e-15L;
protected:
//    elem_t _maxGap = 0.0L;
//    elem_t _maxGapLocation = 0.0L;
    std::ofstream _out;
    Vector newtonRaphson(Solution &sol, Homotopy &homotopy, elem_t t);
public:
    HomotopyContinuation();
    
    ~HomotopyContinuation() {
        if (_out.is_open()) {
            _out.close();
        }
    }
    
    elem_t accuracyLow() {
        return _accuracyLow;
    }
        
    void setAccuracyLow(elem_t accuracy) {
        _accuracyLow = accuracy;
    }
    
    elem_t accuracyHigh() {
        return _accuracyHigh;
    }
    
/*    elem_t maxGap() {
        return _maxGap;
    }
    
    elem_t maxGapLocation() {
        return _maxGapLocation;
    }*/
    
    void setAccuracyHigh(elem_t accuracy) {
        _accuracyHigh = accuracy;
    }

    virtual void trace(elem_t t, Solution *sol, elem_t error) {
        if (_out.is_open()) {
            _out << std::setprecision(6) << t;
            for (int i = 0; i < sol->size(); i++) {
                _out << "\t";
                _out << std::setprecision(18) << chop(sol->get(i));
            }
            _out << "\t" << std::setprecision(6) << error << std::endl;
        }        
    }
    
    void closeTrace() {
        if (_out.is_open()) {
            _out.close();
        }
    }
    
    void setTraceFile(std::string file) {
        if (_out.is_open()) {
            _out.close();
        }
        _out.open(file.c_str());
    }
        
    virtual void solve(Homotopy &homotopy, Solution &sol) = 0;    
};

class SimpleHomotopyContinuation : public HomotopyContinuation
{
public:
    virtual void solve(Homotopy &homotopy, Solution &sol);
};

class BetheHomotopyContinuation : public HomotopyContinuation
{
private:
    Statistics<elem_t> _errors;
    Statistics<elem_t> _gaps;
    Vector _lastChanges[3];
    elem_t _deltas[3];
    int _lastChangeBegin;
    int _lastChangeEnd;
    RootFitter _fitter;
    int _maxNumberOfNewtonRaphonEachStep = 30;
    int _maxDepth = 10;
    
    void initLastChanges();
    void addLastChange(Vector &change, elem_t delta);
    Vector predictChange(elem_t delta);
    Vector nextChange(Homotopy &homotopy, const Solution &sol_in, elem_t t, elem_t delta);
    Vector moveOneStep(Homotopy &homotopy, Solution* sol, elem_t t, elem_t delta);
public:
    BetheHomotopyContinuation();
    bool hasInfinity;

    virtual void solve(Homotopy &homotopy, Solution &sol);
    
    virtual void trace(elem_t t, Solution *sol, elem_t error) {
        HomotopyContinuation::trace(t, sol, error);
        BetheSolution* bs = dynamic_cast<BetheSolution*>(sol);
        if (bs != NULL && bs->indexOfBlowupRoot()>= 0) {
            if (std::abs(bs->omega()) < 1e-2) {
                hasInfinity = true;
            }
        }
    }
    
    Statistics<elem_t>& errors() {
        return _errors;
    }
    
    Statistics<elem_t>& gaps() {
        return _gaps;
    }
    
    int maxNumberOfNewtonRaphonEachStep() {
        return _maxNumberOfNewtonRaphonEachStep;
    }
    
    void setMaxNumberOfNewtonRaphonEachStep(int n) {
        _maxNumberOfNewtonRaphonEachStep = n;
    }
    
    int maxDepth() {
        return _maxDepth;
    }
    
    void setMaxDepth(int val) {
        _maxDepth = val;
    }
};

class HomotopyContinuationException : public std::exception
{
private:
    elem_t _error;
    elem_t _t;
    elem_t _delta;
    std::string message;
public:
    HomotopyContinuationException(elem_t t, elem_t delta, elem_t error);
    virtual const char* what() const noexcept;
};
