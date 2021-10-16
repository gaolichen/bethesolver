#pragma once
#include "common.h"
#include "BetheRootCache.h"
#include "BetheFunctions.h"
#include "Homotopy.h"

class ResultWrapper
{
private:
    Vector *_res;
    bool *_inverted;
    int _size;
    var_t _epsilon;
public:
    ResultWrapper(int size) {
        _size = size;
        _res = new Vector(_size);
        _inverted = new bool[_size];
        for (int i = 0; i < _size; i++) {
            (*_res)[i] = (elem_t)0.0;
            _inverted[i] = false;
        }
        _epsilon = randomComplex(EPS, 2.0 * EPS);
    }
    
    ResultWrapper(const Vector &res) {
        _size = res.size();
        _res = new Vector(_size);
        _inverted = new bool[_size];
        for (int i = 0; i < _size; i++) {
            (*_res)[i] = res[i];
            _inverted[i] = false;
        }
        _epsilon = randomComplex(EPS, 2.0 * EPS);
        update();
    }
        
    ~ResultWrapper() {
        delete _res;
        delete[] _inverted;
    }
    
    int size() {
        return _size;
    }
    
    Vector& res() {
        return *_res;
    }
    
    void apply(Vector& vec) {
        for (int i = 0; i < _size; i++) {
            if (_inverted[i]) {
                vec[i] = (elem_t)1.0 / (vec[i] + _epsilon);
            }
        }
    }
    
    elem_t norm() {
        return _res->norm();
    }
    
    var_t get(int i, bool normal = false) {
        if (!_inverted[i] || !normal) {
            return (*_res)[i];
        } else {
            return (elem_t)1.0 / (*_res)[i];
        }
    }
    
    void update() {
        for (int i = 0; i < _size; i++) {
            if (std::abs((*_res)[i]) < (elem_t)1.0 - EPS) {
                (*_res)[i] = (elem_t)1.0/((*_res)[i] + _epsilon);
                _inverted[i] = !_inverted[i];
            }
        }        
    }

};


class ExpandedBetaHomotopy : public Homotopy
{
private:
    BethePolynomial *_numerator;
    BethePolynomial *_denominator;
    elem_t _beta;
public:
    ExpandedBetaHomotopy(int L, int M, elem_t beta) {
        _beta = beta;
        _numerator = new BethePolynomial(L, M, true);
        _denominator = new BethePolynomial(L, M, false);
    }
    
    virtual elem_t tBegin() {
        return 0.0;
    }
    
    virtual elem_t tEnd() {
        return _beta;
        // when beta * L = Pi, the deformed equation become the undeformed one.
//        return _beta - floor(_beta * _numerator->chainLength()/Pi) * Pi/_numerator->chainLength();
    }
    
    virtual elem_t tolerance() {
        return 1e-8;
    }
        
    virtual elem_t error(elem_t t, const Solution &sol) {
        elem_t phi = 2.0 * _numerator->chainLength() * t;
        Vector num = _numerator->eval(sol);
        Vector deno = _denominator->eval(sol) * var_t(cos(phi), sin(phi));
        ResultWrapper wrap(num);
        wrap.apply(deno);
        return (wrap.res() - deno).norm();
    }
        
    virtual Vector eval(elem_t t, const Solution &sol) {
        elem_t phi = 2.0 * _numerator->chainLength() * t;
        return _numerator->eval(sol) - _denominator->eval(sol) * var_t(cos(phi), sin(phi));
    }
    
    virtual Matrix diff(elem_t t, const Solution &sol) {
        elem_t phi = 2.0 * _numerator->chainLength() * t;
/*        Vector num = _numerator->eval(sol);
        ResultWrapper wrap(num);
        Matrix mat1 = _numerator->diff(sol);
        Matrix mat2 = _denominator->diff(sol) * var_t(cos(phi), sin(phi));
        wrap.apply(mat1);
        wrap.apply(mat2);
        return mat1 - mat2;*/
//        std::cout << "_numerator->diff(sol)=" << _numerator->diff(sol) << std::endl;
//        std::cout << "_denominator->diff(sol)=" << _denominator->diff(sol) << std::endl;
        return _numerator->diff(sol) - _denominator->diff(sol) * var_t(cos(phi), sin(phi));
    }
    
    // calculate derivative of ln(u_j - u_k \pm i) w.r.t u_j
    static var_t diffJ_singular(const var_t &uj, const var_t &uk, bool positiveI, bool kInverted) {
        var_t iUnit = var_t((elem_t)0.0, positiveI ? (elem_t)1.0 : (elem_t)-1.0);
        if (kInverted) {
            return uk/(uj * uk - (elem_t)1.0 + uk * iUnit);
        } else {
            return (elem_t)1.0/(uj - uk + iUnit);
        }
    }
    
    // calculate derivative of ln(u_j - u_k \pm i) w.r.t u_k
    static var_t diffK_singular(const var_t &uj, const var_t &uk, bool positiveI, bool kInverted) {
        var_t iUnit = var_t((elem_t)0.0, positiveI ? (elem_t)1.0 : (elem_t)-1.0);
        if (kInverted) {
            return (uj + iUnit)/(uj * uk - (elem_t)1.0 + uk * iUnit);
        } else {
            return (elem_t)-1.0/(uj - uk + iUnit);
        }
    }
    
    static bool samePhase(var_t c1, var_t c2) {
        return floatMod(std::arg(c1) - std::arg(c2), Pi) < 1e-3;
    }
    
    virtual Vector impliciteDiff(elem_t t, const Solution &sol_in) {
        const InvertableSolution& sol = dynamic_cast<const InvertableSolution&>(sol_in);
        double L = (double)_numerator->chainLength();
        Vector rhs = _numerator->eval(sol) * var_t((elem_t)0.0, 2.0 * L);
        Matrix tmp = diff(t, sol);
        const Vector& u = sol.root();
        int infities = 0;
        elem_t delta = floatMod(2.0 * L * t, 2.0 * Pi);
        bool nearCritical = std::abs(delta) < 1.5 * 2.0 * L * std::abs(_beta/getSteps());
        var_t d1 = 2.0 * (elem_t)L /var_t(L-1.0, -sqrt(L-1.0));
        var_t d2 = 2.0 * (elem_t)L /var_t(L-1.0, sqrt(L-1.0));
        bool d1found = false;
        bool d2found = false;
        bool singularFound = false;

        if (nearCritical) {
            std::cout << std::setprecision(18) << "delta=" << delta << ", _beta/step=" << std::setprecision(6) << 1.5 * 2.0 * L * std::abs(_beta/getSteps()) << std::endl;
        }
        for (int j = 0; j < sol.size(); j++) {
            if (nearCritical && sol.inverted(j) && (samePhase(sol.root(j), d1) || samePhase(sol.root(j), d2))) {
                infities++;
                rhs[j] = (elem_t)1.0;
                std::cout << "found infinities t=" << t << std::endl;
                for (int k = 0; k < sol.size(); k++) {
                    if (j != k) {
                        tmp(j, k) = (elem_t)0.0;
                    }
                }
                if (samePhase(sol.root(j), d1) && !d1found) {
                    tmp(j, j) = (elem_t)1.0 / d1;
                    d1found = true;
                } else {
                    tmp(j, j) = (elem_t)1.0 / d2;
                    d2found = true;
                }
            } else if (sol.isSingular(j)) {
                singularFound = true;
                std::cout << "singular found!" << std::endl;
                rhs[j] = var_t(0.0, (elem_t)(2.0 * L));
                tmp(j, j) = var_t(0.0, -(elem_t)L);
                for (int k = 0; k < sol.size(); k++) {
                    if (j == k) continue;
                    if (sol.isSingular(k)) {
                        tmp(j, k) = (elem_t)0.0;
                        continue;
                    }
                    tmp(j, k) = diffK_singular(u[j], u[k], false, sol.inverted(k)) -
                        - diffK_singular(u[j], u[k], true, sol.inverted(k));
                    tmp(j, j) += diffJ_singular(u[j], u[k], false, sol.inverted(k)) -
                        diffJ_singular(u[j], u[k], true, sol.inverted(k));
                }
            }
        }
        if ((d1found && d2found) || singularFound) {
            std::cout << "sol=" << sol.root() << std::endl;
            std::cout << "rhs=" << rhs << ", tmp.inverse() * rhs=" << tmp.inverse() * rhs << std::endl;
        }
//        std::cout << "diff(t, sol).inverse()=" << tmp.inverse() << std::endl;
        return tmp.inverse() * rhs;
    }
    
    var_t _e2ipStart;
    
    virtual void setStartRoot(const Vector& root) {
        Homotopy::setStartRoot(root);
        _e2ipStart = e2ip(root);
    }
        
    static var_t e2ip2(Vector& root, const InvertableSolution &sol) {
        var_t ret((elem_t)1.0);
        bool hasSingular = false;
        for (int i = 0; i < root.size(); i++) {
            if (!sol.isSingular(i)) {
                if (!sol.inverted(i)) {
                    ret *= (root[i] + var_t(0.0, 0.5))/(root[i] - var_t(0.0, 0.5));
                } else {
                    ret *= ((elem_t)1.0 + var_t(0.0, 0.5) * root[i])/((elem_t)1.0 - var_t(0.0, 0.5) * root[i]);
                }
            } else {
                hasSingular = true;
            }
        }
        if (hasSingular) {
            ret = -ret;
        }
        return ret;
    }

    virtual bool acceptCorrection(elem_t t, const Solution &sol_in, Vector& correction) {
        const InvertableSolution &sol = dynamic_cast<const InvertableSolution&>(sol_in);
        Vector r = sol.root() + correction;
        int m = _numerator->magnonNumber();
        var_t expected = std::polar((elem_t)1.0, t* 2.0 * m) * _e2ipStart;
        var_t ip = e2ip2(r, sol);
        if (std::abs(ip - expected) > 0.1) {
            std::cout << "ip=" << std::log(ip)/(elem_t)Pi << ", expected=" << std::log(expected)/(elem_t)Pi << " diff=" << std::abs(ip - expected) << ", correction=" << correction << std::endl;
            return false;
        }
        return true;
    }
};

class PoleFreeBetheHomotopy : public Homotopy
{
private:
    int _L;
    int _M;
    elem_t _beta;
    BetheFactorBase *_left;
    BetheFactorBase *_right;
public:
    PoleFreeBetheHomotopy(int L, int M, elem_t beta) {
        _L = L;
        _M = M;
        _beta = beta;
        _left = new BetheFactorLeft(L, M);
        _right = new BetheFactorRight(L, M);
        setSteps((int)floor(L * 200 * beta / Pi + 0.5));
    }
    
    ~PoleFreeBetheHomotopy() {
        delete _left;
        delete _right;
    }
    
    BetheFactorBase* Left() { return _left; }
    BetheFactorBase* Right() { return _right; }
    
    virtual elem_t tBegin() {
        return 0.0L;
    }
    
    virtual elem_t tEnd() {
        return _beta;
    }
    
    virtual elem_t tolerance() {
        return 1e-15;
    }
    
    virtual Vector eval(elem_t t, const Solution &sol) {
        elem_t phi = 2.0L * _L * t;
        elem_t phi2 = 2.0L * _M * t;
        Vector ret = _left->eval(sol) - _right->eval(sol) * std::polar(1.0L, phi);
        ret(ret.size() - 1) = _left->evalMomentum(sol) - _right->evalMomentum(sol) * _e2ipStart * std::polar(1.0L, phi2);
        return ret;
    }
    
    virtual Matrix diff(elem_t t, const Solution &sol) {
        elem_t phi = 2.0L * _L * t;
        elem_t phi2 = 2.0L * _M * t;
        Matrix ret = _left->diff(sol) - _right->diff(sol) * std::polar(1.0L, phi);
        Vector lastRow = _left->diffMomentum(sol) - _right->diffMomentum(sol) * _e2ipStart * std::polar(1.0L, phi2);
        for (int i = 0; i < ret.cols(); i++) {
            ret(ret.rows() - 1, i) = lastRow(i);
        }
        
        return ret;
    }
    
    elem_t minJacobinNorm = 1000.0L;
    
    virtual Vector impliciteDiff(elem_t t, const Solution &sol) {
        elem_t phi = 2.0L * _L * t;
        elem_t phi2 = 2.0L * _M * t;
//        Vector rhs = _left->eval(sol) * var_t((elem_t)0.0, 2.0 * _L);
        Vector rhs = _right->eval(sol) * var_t((elem_t)0.0, 2.0 * _L) * std::polar(1.0L, phi);
        rhs(rhs.size() - 1) = _right->evalMomentum(sol) * var_t((elem_t)0.0, 2.0 * _M) * std::polar(1.0L, phi2);
        Matrix jacobian = diff(t, sol);
        if (std::abs(jacobian.determinant()) < minJacobinNorm) {
            minJacobinNorm = std::abs(jacobian.determinant());
//            std::cout << "jacobian=" << jacobian << std::endl;
            LOG("t=" << t << ", minJacobinNorm=" << std::setprecision(15) << minJacobinNorm, Info);
        }
//        std::cout << "jacobian=" << jacobian << std::endl;
//        std::cout << "jacobian.determinant=" << jacobian.determinant() << std::endl;
//        std::cout << "jacobian.inverse()=" << jacobian.inverse() << std::endl;
/*        if (std::abs(t-0.735132680940012) < 1e-3) {
            BetheSolution& bs = dynamic_cast<BetheSolution&>(sol);
            std::cout << "t=" << t << ", sol.root()=" << std::setprecision(18) << bs.root() << std::endl;
            std::cout << "sol.normalRoot()=" << bs.normalRoot() << std::endl;
            std::cout << "sol.indexOfBlowupRoot()=" << bs.indexOfBlowupRoot() << std::endl;
            std::cout << "sol.omega()=" << bs.omega() << std::endl;
            std::cout << "sol.indexSingularRoot1()=" << bs.indexOfSingularRoot1() << std::endl;
            std::cout << "jacobian=" << jacobian << std::endl;
            std::cout << "jacobian.inverse()=" << jacobian.inverse() << std::endl;
            std::cout << "jacobian.determinant=" << jacobian.determinant() << std::endl;
        }*/

        return jacobian.inverse() * rhs;
    }
    
    var_t _e2ipStart;
    
    virtual void setStartRoot(const Vector& root) {
        Homotopy::setStartRoot(root);
        _e2ipStart = e2ip(root);
    }
    
    static var_t e2ip2(const Vector& root, const BetheSolution &sol) {
        var_t ret((elem_t)1.0);
        bool hasSingular = false;
        if (sol.indexOfSingularRoot2() >= 0 && std::abs(sol.epsilon()) < EPS) {
            hasSingular = true;
            ret = -1.0L;
        }
        for (int i = 0; i < root.size(); i++) {
            if (sol.isSingular(i) && hasSingular) continue;
            var_t u = sol.get(i);
            ret *= (u + 0.5il)/(u - 0.5il);
        }
        return ret;
    }
    
    virtual bool acceptCorrection(elem_t t, const Solution &sol_in, const Vector& correction) {
        const BetheSolution &sol = dynamic_cast<const BetheSolution&>(sol_in);
        Vector r = sol.root() + correction;
        var_t expected = std::polar((elem_t)1.0, t * 2.0 * _M) * _e2ipStart;
        var_t ip = e2ip2(r, sol);
        if (std::abs(ip - expected) > 0.1) {
//            std::cout << "ip=" << std::log(ip)/(elem_t)Pi << ", expected=" << std::log(expected)/(elem_t)Pi << " diff=" << std::abs(ip - expected) << ", correction=" << correction << ", roots=" << sol.root() << std::endl;
            return false;
        }
        return true;
    }
    
    virtual elem_t error(elem_t t, const Solution &sol_in) {
        const BetheSolution &sol = dynamic_cast<const BetheSolution&>(sol_in);
        elem_t phi = 2.0 * _L * t;
        var_t expected = std::polar(1.0L, phi);
        Vector num = _left->eval(sol);
//        Vector deno = _right->eval(sol) * var_t(cos(phi), sin(phi));
        Vector deno = _right->eval(sol);
        Vector res(num.size());
        var_t largestError = .0L;
        for (int i = 0; i < res.size(); i++) {
            if (sol.isIdenticalRoot(i)) {
                res(i) = 0.0L;
            } else {
                res(i) = num(i)/deno(i) - expected;
                if (sol.isLargestRoot(i)) {
                    largestError = res(i);
                }
            }
        }
        
        for (int i = 0; i < res.size(); i++) {
            if (sol.isIdenticalRoot(i)) {
                res(i) = largestError;
            }
        }
        return std::sqrt(res.squaredNorm()/_M);
//        ResultWrapper wrap(num);
//        wrap.apply(deno);
//        return (wrap.res() - deno).norm();
    }
};
