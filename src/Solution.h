#pragma once
#include "common.h"
#include "RootFitter.h"
#include <complex>
using namespace std::literals::complex_literals;

class Solution
{
protected:
    Vector *_root;
    int _size;
public:
    Solution(const Solution &sol) {
        _size = sol._size;
        _root = new Vector(_size);
        for (int i = 0; i < _size; i++) {
            (*_root)[i] = (*sol._root)(i);
        }
    }
    
    Solution(const Vector &root) {
        _size = root.size();
        _root = new Vector(root.size());
        for (int i = 0; i < _size; i++) {
            (*_root)[i] = root[i];
        }
    }
    
    virtual Solution* copy() const {
        return new Solution(*this);
    }
        
    ~Solution() {
        delete _root;
    }
    
    void setRoot(const Vector & newRoot) {
        for (int i = 0; i < size(); i++) {
            (*_root)[i] = newRoot[i];
        }
    }
    
    virtual bool update(elem_t t, RootFitter* fitter = NULL) {
        return false;
    }
    
    int size() const {
        return _size;
    }
    
    Vector& root() const {
        return *_root;
    }
    
    var_t& root(int i) const {
        return (*_root)[i];
    }
    
    Vector normalRoot() const {
        Vector ret(size());
        for (int i = 0; i < ret.size(); i++) {
            ret[i] = this->get(i);
        }
        
        return ret;
    }
    
    virtual var_t get(int i) const {
        return (*_root)[i];
    }
    
    virtual void output() const {
    }
};

extern elem_t eps_may_be_singular;
extern elem_t blowup_root_cutoff;
bool mayBeSingular(var_t u);
bool cannotBeSingular(var_t &u);

class SingularRootSolution : public Solution
{
protected:
    int _indexOfSingularRoot1 = -1;
    int _indexOfSingularRoot2 = -1;
    int _L;
    var_t _zero = 0.0L;
    bool setupSingularRoots(elem_t beta);
public:
    SingularRootSolution(const SingularRootSolution &sol) : Solution(sol) {
        _indexOfSingularRoot1 = sol._indexOfSingularRoot1;
        _indexOfSingularRoot2 = sol._indexOfSingularRoot2;
        _L = sol._L;
    }
    
    SingularRootSolution(const Vector &root, int L) : Solution(root) {
        _L = L;
    }
    
    virtual Solution* copy() const {
        return new SingularRootSolution(*this);
    }
    
    int chainLength() const { return _L; }
    
    int indexOfSingularRoot1() const {
        return _indexOfSingularRoot1;
    }
    
    int indexOfSingularRoot2() const {
        return _indexOfSingularRoot2;
    }
    
    bool isSingular(int j) const {
        return j == _indexOfSingularRoot2 || j == _indexOfSingularRoot1;
    }
    
    virtual bool isRegular(int j) const {
        return j != _indexOfSingularRoot2 && j != _indexOfSingularRoot1;
    }
    
    const var_t& c() const {
        if (_indexOfSingularRoot1 < 0) {
            return _zero;
        }
        return root()[_indexOfSingularRoot1];
    }
    
    const var_t& epsilon() const {
        if (_indexOfSingularRoot2 < 0) {
            return _zero;
        }
        return root()[_indexOfSingularRoot2];
    }
    
    // return true value of roots.
    virtual var_t get(int index) const {
        if (index == _indexOfSingularRoot2) {
            return -0.5il + root(index);
        } else if (index == _indexOfSingularRoot1) {
            return 0.5il + root(_indexOfSingularRoot2) + root(index) * std::pow(root(_indexOfSingularRoot2), _L);
        }
        
        return Solution::get(index);
    }
    
    // return get(index) if it's singular, otherwise return root(index);
    var_t get2(int index) const {
        if (isSingular(index)) {
            return get(index);
        } else {
            return root(index);
        }
    }
    
    virtual void output() const {
        LOG("_indexOfSingularRoot1=" << _indexOfSingularRoot1 << ", _indexOfSingularRoot2=" << _indexOfSingularRoot2 << ", eps=" << epsilon() << ", c=" << c(), Info);
    }
};

class InvertableSolution : public SingularRootSolution
{
private:
    bool *_inverted;
public:
    InvertableSolution(const InvertableSolution& sol) : SingularRootSolution(sol) {
        _inverted = new bool[sol._size];
        for (int i = 0; i < _size; i++) {
            _inverted[i] = sol._inverted[i];
        }
    }
    
    InvertableSolution(Vector &root, int L) : SingularRootSolution(root, L) {
        _inverted = new bool[root.size()];
        for (int i = 0; i < _size; i++) {
            _inverted[i] = false;
        }
        update(0.0L);
    }
    
    virtual Solution* copy() const {
        return new InvertableSolution(*this);
    }
    
    ~InvertableSolution() {
        delete[] _inverted;
    }
    
    virtual bool isRegular(int j) const {
        return !_inverted[j];
    }
    
    virtual var_t get(int i) const {
        if (_inverted[i]) {
            return (elem_t)1.0 / (*_root)[i];
        }
        return SingularRootSolution::get(i);
    }    

    bool inverted(int i) const {
        return _inverted[i];
    }
    
//    bool isSingular(int i) {
//        elem_t eps = 1e-8;
//        return !_inverted[i] && std::abs((*_root)[i].real()) < 1e-2 &&
//        (std::abs((*_root)[i].imag() + 0.5) < EPS || std::abs((*_root)[i].imag() - 0.5) < EPS);
//    }
    
    virtual bool update(elem_t t, RootFitter* fitter = NULL) {
        bool ret = false;
        for (int i = 0; i < _size; i++) {
            if (!isSingular(i) && std::abs((*_root)[i]) > (elem_t)(1.0 + EPS)) {
                (*_root)[i] = (elem_t)1.0/(*_root)[i];
                _inverted[i] = !_inverted[i];
                ret = true;
                LOG("t=" << t << ", inverted roots changed", Info);
                if (fitter != NULL) {
                    fitter->invert(i);
                }
            }
        }
        
        if (setupSingularRoots(t)) {
            if (fitter != NULL) {
                fitter->clear();
                LOG("t=" << t << " RootFitter.clear()", Info);
            }
            ret = true;
        }
        return ret;
    }    
};


enum RootType {
    regular = 0, // for regular type root, we store it with its original value
    largest = 1, // for max type root, we store it as recipcal
    large = 2    // for large type root, we store the ratio between its value and max root.
};

class BetheSolution : public SingularRootSolution
{
private:
//    static const elem_t eps_may_be_singular = 1e-1;
    int _indexOfBlowupRoot = -1;
    elem_t _stepSize;
    std::vector<RootType>    _types;
    std::vector<int>    _identical;
    
    static elem_t magnitude(const var_t &u) {
        return abs(u.real()) > abs(u.imag()) ? abs(u.real()) : abs(u.imag());
    }
    
    bool setupBlowupRoots(elem_t beta);
//    bool setupSingularRoots(elem_t beta);
    
public:

    BetheSolution(const BetheSolution& sol) : SingularRootSolution(sol) {
        _indexOfBlowupRoot = sol._indexOfBlowupRoot;
        _stepSize = sol._stepSize;
        _types = sol._types;
        _identical = sol._identical;
    }
    
    BetheSolution(Vector &root, int L, elem_t beta, elem_t stepSize = 0.01L) : SingularRootSolution(root, L) {
        _stepSize = stepSize;
        _types.resize(root.size(), regular);
        _identical.resize(root.size(), 0);
        update(beta);
    }
    
    virtual Solution* copy() const {
        return new BetheSolution(*this);
    }
    
    virtual bool update(elem_t beta, RootFitter* fitter = NULL);
    
    virtual bool isRegular(int j) const {
        return _types[j] == regular;
    }
        
    bool isLargestRoot(int j) const {
        return j == _indexOfBlowupRoot;
    }
    
    bool isIdenticalRoot(int j) const {
        return _identical[j] == 1;
    }
    
    int indexOfBlowupRoot() const {
        return _indexOfBlowupRoot;
    }
    
    const var_t& omega() const {
        if (_indexOfBlowupRoot == -1) {
            return _zero;
        } else {
            return root(_indexOfBlowupRoot);
        }
    }
    
    RootType type(int index) const {
        return _types[index];
    }
        
    // return true value of roots.
    virtual var_t get(int index) const {        
        if (_types[index] == largest) {
            return (elem_t)1.0 / root(index);
        } else if (_types[index] == large) {
            return root(index) / root(_indexOfBlowupRoot);
        }
        
        return SingularRootSolution::get(index);
    }
        
    virtual void output() const {
        LOG("_types=" << _types << ",_identical=" << _identical << " _indexOfBlowupRoot=" << _indexOfBlowupRoot << ", _indexOfSingularRoot1=" << _indexOfSingularRoot1 << ", _indexOfSingularRoot2=" << _indexOfSingularRoot2, Info);
    }
};

