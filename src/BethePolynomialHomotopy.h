#pragma once
#include "common.h"
#include "BethePolynomials.h"
#include "Homotopy.h"

class BethePolynomialHomotopy : public Homotopy
{
private:
    int _L;
    int _M;
    elem_t _beta;
    BethePolynomialBase *_left;
    BethePolynomialBase *_right;
public:
    BethePolynomialHomotopy(int L, int M, elem_t beta) {
        _L = L;
        _M = M;
        _beta = beta;
        _left = new BethePolynomialLeft(L, M);
        _right = new BethePolynomialRight(L, M);
        setSteps((int)floor(L * 200 * abs(beta) / Pi + 0.5));
    }
    
    ~BethePolynomialHomotopy() {
        delete _left;
        delete _right;
    }
    
    BethePolynomialBase* Left() { return _left; }
    BethePolynomialBase* Right() { return _right; }
    
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
    
    static var_t e2ip2(const Vector& root, const InvertableSolution &sol) {
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
        const InvertableSolution &sol = dynamic_cast<const InvertableSolution&>(sol_in);
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
        const InvertableSolution &sol = dynamic_cast<const InvertableSolution&>(sol_in);
        elem_t phi = 2.0 * _L * t;
        var_t expected = std::polar(1.0L, phi);
        Vector num = _left->eval(sol);
//        Vector deno = _right->eval(sol) * var_t(cos(phi), sin(phi));
        Vector deno = _right->eval(sol);
        Vector res(num.size());
        for (int i = 0; i < res.size(); i++) {
            res(i) = num(i)/deno(i) - expected;
        }
        
        elem_t ret = std::sqrt(res.squaredNorm()/_M);
        if (ret > 0.1L) {
//            LOG("ret=" << ret << std::endl << "res=" << res, Info);
//            LOG("num=" << num << ", deno=" << deno, Info);
        }
        return ret;
//        ResultWrapper wrap(num);
//        wrap.apply(deno);
//        return (wrap.res() - deno).norm();
    }
};
