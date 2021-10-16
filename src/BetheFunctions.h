#pragma once
#include "common.h"
#include "BetheRootCache.h"
#include "Function.h"

class BethePolynomial : public Function
{
private:
    var_t _iUnit;
    int _L;
    int _M;
public:
    BethePolynomial(int L, int M, bool isNumerator);
    
    int chainLength() { return _L; }
    
    int magnonNumber() { return _M; }
    
    bool isSingular(var_t u) {
        return abs(u.real()) < EPS && (abs(u.imag() + 0.5) < EPS || abs(u.imag() - 0.5) < EPS);
    }
    
    bool isSingular(Vector u) {
        for (int i = 0; i < u.size(); i++) {
            if (isSingular(u[i])) return true;
        }
        return false;
    }
    
    var_t factorJK(const var_t& uj, const var_t &uk, bool invertj, bool invertk) {
        if (!invertj && !invertk) {
            return uj - uk - _iUnit;
        } else if (!invertk) {
            return (elem_t)1.0 - uj * uk - _iUnit * uj;
        } else if (!invertj) {
            return uj * uk - (elem_t)1.0 - _iUnit * uk;
        } else {
            return uk - uj - _iUnit * uk * uj;
        }
    }
    
    
    virtual Vector eval(const Solution &sol_in) {
        const InvertableSolution& sol = dynamic_cast<const InvertableSolution&>(sol_in);
        Vector ret(sol.size());
//        Vector& u = sol.root();
        for (int j = 0; j < sol.size(); j++) {
            if (sol.isSingular(j)) {
                ret[j] = (elem_t)0.0;
            } else {
                if (sol.inverted(j)) {
                    ret[j] = std::pow((elem_t)1.0 + _iUnit * sol.root(j) * (elem_t).5, _L);
                } else {
                    ret[j] = std::pow(sol.root(j) + _iUnit * (elem_t).5, _L);
                }
                for (int k = 0; k < sol.size(); k++) {
                    if (k != j) {
                        ret[j] *= factorJK(sol.root(j), sol.root(k), sol.inverted(j), sol.inverted(k));
                    }
                }
            }
        }
        return ret;
    }
    
    var_t diffJ(const var_t &uj, const var_t &uk, bool invertj, bool invertk) {
        if (!invertj && !invertk) {
            return (elem_t)1.0;
        } else if (!invertk) {
            // if only uj is inverted
            return -uk - _iUnit;
        } else if (!invertj) {
            // if only uk is inverted
            return uk;
        } else {
            // if both are inverted
            return (elem_t)-1.0 - _iUnit * uk;
        }
    }
    
    var_t diffK(const var_t &uj, const var_t &uk, bool invertj, bool invertk) {
        if (!invertj && !invertk) {
            return (elem_t)-1.0;
        } else if (!invertk) {
            // if only uj is inverted
            return -uj;
        } else if (!invertj) {
            // if only uk is inverted
            return uj - _iUnit;
        } else {
            // if both are inverted
            return (elem_t)1.0 - _iUnit * uj;
        }
    }
    
    // return a matrix m where m(i, j) = d f_i/d x_j
    virtual Matrix diff(const Solution &sol_in) {
        const InvertableSolution& sol = dynamic_cast<const InvertableSolution&>(sol_in);
        Matrix ret(sol.size(), sol.size());
        const Vector &u = sol.root();
        var_t prefactor;
        for (int j = 0; j < ret.rows(); j++) {
            // for singluar solutions we need special treatment
            if (sol.isSingular(j)) {
                for (int k = 0; k < u.size(); k++) {
                    ret(j, k) = (elem_t)0.0;
                }
                continue;
            }
            // handle non-singular ones.
            if (sol.inverted(j)) {
                prefactor = std::pow((elem_t)1.0 + _iUnit * (elem_t)0.5 * u[j], _L);
                ret(j, j) = std::pow((elem_t)1.0 + _iUnit * (elem_t)0.5 * u[j], _L - 1) * _iUnit * (elem_t)(_L * 0.5);
            } else {
                prefactor = std::pow(u[j] + _iUnit * (elem_t)0.5, _L);
                ret(j, j) = std::pow(u[j] + _iUnit * (elem_t)0.5, _L - 1) * (elem_t)_L;
            }
            
            for (int k = 0; k < u.size(); k++) {
                if (k == j) continue;
                ret(j, j) *= factorJK(u[j], u[k], sol.inverted(j), sol.inverted(k));
            }
            
            for (int k = 0; k < u.size(); k++) {
                if (k == j) continue;
                
                ret(j, k) = prefactor;
                for (int h = 0; h < u.size(); h++) {
                    if (h == j || h == k) continue;
                    ret(j, k) *= factorJK(u[j], u[h], sol.inverted(j), sol.inverted(h));
                }
                ret(j, j) += ret(j, k) * diffJ(u[j], u[k], sol.inverted(j), sol.inverted(k));
                ret(j, k) *= diffK(u[j], u[k], sol.inverted(j), sol.inverted(k));
            }
        }
        
        return ret;
    }
};

extern var_t prefactorCoefLeft[3][3];

extern var_t prefactorCoefRight[3][3];

extern var_t coefLeft[3][3][6];

extern var_t coefRight[3][3][6];

class BetheFactorBase : public Function
{
protected:
    int _L;
    int _M;
public:
    BetheFactorBase(int L, int M);
    
    int chainLength() { return _L; }
    
    int magnonNumber() { return _M; }
    
    virtual var_t* coef1(RootType type) = 0;
    virtual var_t* coef2(RootType jType, RootType kType) = 0;
        
    virtual var_t halfIProductSingular1(const BetheSolution& sol) = 0;
    virtual var_t halfIProductSingular2(const BetheSolution& sol) = 0;
    virtual var_t momentumHalfIFactorSingular(const BetheSolution& sol) = 0;
    
    virtual var_t halfIProductDc1(const BetheSolution& sol) = 0;
    virtual var_t halfIProductDc2(const BetheSolution& sol) = 0;
    virtual var_t halfIProductMomentumDc(const BetheSolution& sol) = 0;

    virtual var_t halfIProductDeps1(const BetheSolution& sol) = 0;
    virtual var_t halfIProductDeps2(const BetheSolution& sol) = 0;
    virtual var_t halfIProductMomentumDeps(const BetheSolution& sol) = 0;

    
    // returns (uj +- i/2)^L
    var_t halfIProduct(const BetheSolution& sol, int j) {
        if (j == sol.indexOfSingularRoot1()) {
            return halfIProductSingular1(sol);
        } else if (j == sol.indexOfSingularRoot2()) {
            return halfIProductSingular2(sol);
        } else {
            return std::pow(halfIFactor(sol, j), _L);
        }
    }
    
    // returns (uj +- i/2)
    var_t halfIFactor(const BetheSolution& sol, int j) {
        var_t* c = coef1(sol.type(j));
        return c[0] + c[1] * sol.root(j) + c[2] * sol.omega();
    }
    
    // derivative of (uj +- i/2) with respect to omega
    var_t halfIFactorDiffOmega(const BetheSolution& sol, int j) {
        var_t* c = coef1(sol.type(j));
        return c[2];
    }
    
    // derivative of (uj +- i/2) with respect to uj
    var_t halfIFactorDiffJ(const BetheSolution& sol, int j) {
        var_t* c = coef1(sol.type(j));
        if (sol.isLargestRoot(j)) {
            return c[2];
        } else {
            return c[1];
        }
    }
    
    virtual var_t identicalRootFactor() = 0;

    // product of (uj - uk +- i) for all k != j and k != toSkip
    // if uj is singular, it only nonsigular uk is included. 
    var_t factorProduct(const BetheSolution& sol, int j, int toSkip = -1) {
        var_t ret = 1.0L;
        var_t uj = sol.get2(j);
        bool isSingular = sol.isSingular(j);
        
        for (int k = 0; k < sol.size(); k++) {
            if (k == j || k == toSkip) {
                continue;
            }
            if (isSingular && sol.isSingular(k)) {
                continue;
            }
            if (sol.isLargestRoot(j) && sol.isIdenticalRoot(k)) {
                ret *= identicalRootFactor();
                continue;
            }
            ret *= factorJK(uj, sol.get2(k), sol.omega(), sol.type(j), sol.type(k));
        }
        return ret;
    }
    
    // factor uj - uk +- i
    var_t factorJK(const BetheSolution& sol, int j, int k) {
        return factorJK(sol.get2(j), sol.get2(k), sol.omega(), sol.type(j), sol.type(k));
    }
    
    // factor uj - uk +- i
    var_t factorJK(var_t uj, var_t uk, const var_t &omega, RootType typeJ, RootType typeK) {
//        var_t vec[6] = {1.0L, uj, uk, sol.omega(), sol.omega() * uj, sol.omega() * uk};
        var_t *c = coef2(typeJ, typeK);
        return 1.0L * c[0] + uj * c[1] + uk * c[2] + c[3] * omega + c[4] * omega * uj + c[5] * omega * uk;
    }
    
    // derivative of the (uj-uk +- i) factor with respect to uk.
    var_t diffK(const BetheSolution& sol, int j, int k) {
        var_t *c = coef2(sol.type(j), sol.type(k));
        // if k-th root is the largest one, which equals 1/omega,
        // we need take directive with respect to omega
        if (sol.isLargestRoot(k)) {
            // here c[5] should be zero, so we do not need the c[5] term. 
            return c[3] + c[4] * sol.get2(j);
        } else {
            return c[2] + c[5] * sol.omega();
        }
    }
    
    // derivative of the (uj-uk +- i) factor with respect to uj.
    var_t diffJ(const BetheSolution& sol, int j, int k) {
        var_t *c = coef2(sol.type(j), sol.type(k));
        // if j-th root is the largest one, which equals 1/omega,
        // we need take directive with respect to omega
        if (sol.isLargestRoot(j)) {
            // c[4] should be zero, so we do not need c[4] term
            return c[3] + c[5] * sol.get2(k);
        } else {
            return c[1] + c[4] * sol.omega();
        }
    }
    
    // derivative of the (uj-uk +- i) factor with respect to omega.
    var_t diffOmega(const BetheSolution& sol, int j, int k) {
        var_t *c = coef2(sol.type(j), sol.type(k));
        return c[3] + c[4] * sol.get2(j) + c[5] * sol.get2(k);
    }
    
    virtual Vector eval(const Solution &sol_in) {
        const BetheSolution& sol = dynamic_cast<const BetheSolution &>(sol_in);
        Vector ret(sol.size());
        for (int j = 0; j < sol.size(); j++) {
            if (sol.isIdenticalRoot(j)) {
                ret[j] = 0.0L;
            } else {
                ret[j] = halfIProduct(sol, j) * factorProduct(sol, j);
            }
        }
        return ret;
    }
    
    var_t evalMomentum(const Solution &sol_in) {
        const BetheSolution& sol = dynamic_cast<const BetheSolution &>(sol_in);
        var_t ret = 1.0L;
        if (sol.indexOfSingularRoot1() >= 0 && sol.indexOfSingularRoot2() >= 0) {
            ret = momentumHalfIFactorSingular(sol);
        }

        for (int j = 0; j < sol.size(); j++) {
            if (!sol.isSingular(j)) {
                ret *= halfIFactor(sol, j);
            }
        }
        
        return ret;
    }
    
    Vector diffMomentum(const Solution &sol_in) {
        const BetheSolution& sol = dynamic_cast<const BetheSolution &>(sol_in);
//        Vector ret(sol.size());
        Vector ret = Vector::Zero(sol.size());
        var_t factor = 1.0L;
        for (int j = 0; j < sol.size(); j++) {
            if (!sol.isSingular(j)) {
                factor *= halfIFactor(sol, j);
            }
        }
        
        // handle singular roots
        if (sol.indexOfSingularRoot1() >= 0 && sol.indexOfSingularRoot2() >= 0) {
            int i1 = sol.indexOfSingularRoot1();
            int i2 = sol.indexOfSingularRoot2();
            ret(i1) = halfIProductMomentumDc(sol) * factor;
            ret(i2) = halfIProductMomentumDeps(sol) * factor;
            
            factor *= momentumHalfIFactorSingular(sol);
        }

        for (int j = 0; j < sol.size(); j++) {
            if (!sol.isSingular(j)) {
                ret(j) = factor * halfIFactorDiffJ(sol, j) / halfIFactor(sol, j);
            }
        }
        
        int k = sol.indexOfBlowupRoot();
        if (k >= 0) {
            // for the index of blow up root, there are additional terms.
            for (int h = 0; h < sol.size(); h++) {
                if (sol.type(h) != large) continue;
                ret(k) += factor * halfIFactorDiffOmega(sol, h) / halfIFactor(sol, h);
            }
        }
        
        return ret;
    }
    
    virtual var_t identicalRootDerivative() = 0;
    
    virtual Matrix diff(const Solution &sol_in) {
        const BetheSolution& sol = dynamic_cast<const BetheSolution &>(sol_in);
        Matrix ret(sol.size(), sol.size());
        
        // step 1: compute off-diagonal terms of entry (j,k) execept for:
        // 1. (j,k) are a pair of singular roots
        for (int j = 0; j < sol.size(); j++) {
            var_t prod = factorProduct(sol, j);
            var_t prefactor = halfIProduct(sol, j);
            for (int k = 0; k < sol.size(); k++) {
//                ret(j, k) = 0.0L;
                if (j == k) continue;
                if (sol.isIdenticalRoot(j) || sol.isIdenticalRoot(k)) {
                    ret(j, k) = 0.0L;
                    continue;
                }
                
                if (sol.isSingular(j) && sol.isSingular(k)) continue;
                
                if (!isZero(prod)) {
                    ret(j, k) = diffK(sol, j, k) / factorJK(sol, j, k);
                } else {
                    ret(j, k) = diffK(sol, j, k) * factorProduct(sol, j, k);
                }
                
                // if k is the index of blow up root, there are additional terms.
                if (k == sol.indexOfBlowupRoot()) {
                    for (int h = 0; h < sol.size(); h++) {
                        if (h == j || h == k || (sol.type(h) != large && sol.type(j) != large)) continue;
                        if (!isZero(prod)) {
                            ret(j, k) += diffOmega(sol, j, h) / factorJK(sol, j, h);
                        } else {
                            ret(j, k) += diffOmega(sol, j, h) * factorProduct(sol, j, h);
                        }
                    }
                }
                
                if (!isZero(prod)) {
                    ret(j, k) *= prod;
                }                
                ret(j, k) *= prefactor;
                
                // contribution from (uj +- i/2)^L
                if (k == sol.indexOfBlowupRoot() && sol.type(j) == large) {
                    ret(j, k) += (elem_t)_L * halfIFactorDiffOmega(sol, j) *
                        std::pow(halfIFactor(sol, j), _L-1) * prod;
                }
            }
        }
        
        // step 2: compute diagonal terms except for singular roots.
        for (int j = 0; j < sol.size(); j++) {
            //if (sol.isLargestRoot(j) || sol.isSingular(j)) continue;
            if (sol.isSingular(j)) continue;
            if (sol.isIdenticalRoot(j)) {
                ret(j, j) = identicalRootDerivative();
                continue;
            }
            var_t prefactor = halfIProduct(sol, j);
            var_t prod = factorProduct(sol, j);
            ret(j, j) = 0.0L;
            
            // contribute from uj-uk +- i terms
            for (int k = 0; k < sol.size(); k++) {
                if (k == j || (sol.isLargestRoot(j) && sol.isIdenticalRoot(k))) continue;
                if (!isZero(prod)) {
                    ret(j, j) += diffJ(sol, j, k) / factorJK(sol, j, k);
                } else {
                    ret(j, j) += diffJ(sol, j, k) * factorProduct(sol, j, k);
                }
            }
            if (!isZero(prod)) {
                ret(j, j) *= prod;
            }
            ret(j, j) *= prefactor;
            
            // contribute from (uj +- i/2)^L term
            ret(j, j) += prod * prefactor * halfIFactorDiffJ(sol, j) * (elem_t)_L / halfIFactor(sol, j);
        }

        // step 3: compute terms involving singular roots.
        if (sol.indexOfSingularRoot1() >= 0 && sol.indexOfSingularRoot2() >= 0) {
            int j1 = sol.indexOfSingularRoot1();
            int j2 = sol.indexOfSingularRoot2();
            var_t prod1 = factorProduct(sol, j1);
            var_t prod2 = factorProduct(sol, j2);
            var_t prefactor1 = halfIProductSingular1(sol);
            var_t prefactor2 = halfIProductSingular2(sol);
            
            // first find contribution of (uj - uk +- i) factors to diagonal terms.
            ret(j1, j1) = 0.0L;
            ret(j1, j2) = 0.0L;
            ret(j2, j2) = 0.0L;
            ret(j2, j1) = 0.0L;

            for (int k = 0; k < sol.size(); k++) {
                if (sol.isSingular(k)) continue;
                ret(j1, j1) += diffJ(sol, j1, k) / factorJK(sol, j1, k);
                ret(j2, j2) += diffJ(sol, j2, k) / factorJK(sol, j2, k);
            }
            ret(j1, j1) *= prod1 * prefactor1;
            ret(j2, j2) *= prod2 * prefactor2;
            
            // apply chain rule: d/(d eps) = d/du1 * (1 + c*L*eps^(L-1)) + d/du2
            // d/(d c) = d/du1 * eps^L
            var_t du1dEps = 1.0L + sol.c() * (elem_t)_L * std::pow(sol.epsilon(), _L - 1);
            var_t du1dc = std::pow(sol.epsilon(), _L);
            for (int j = 0; j < sol.size(); j++) {
                ret(j, j2) += ret(j, j1) * du1dEps;
                ret(j, j1) *= du1dc;
            }
            
            // compute contribution from (uj +- i/2)^L factors
            ret(j1, j1) += halfIProductDc1(sol) * prod1;
            ret(j1, j2) += halfIProductDeps1(sol) * prod1;
            ret(j2, j1) += halfIProductDc2(sol) * prod2;
            ret(j2, j2) += halfIProductDeps2(sol) * prod2;
        }
        
        return ret;        
    }
};

class BetheFactorLeft : public BetheFactorBase
{
private:
public:
    BetheFactorLeft(int L, int M) : BetheFactorBase(L, M) {
    }
    
    virtual var_t identicalRootFactor() {
        return -1.0L;
    }
    
    virtual var_t identicalRootDerivative() {
        return 1.0L;
    }
    
    virtual var_t* coef1(RootType type) {
        return prefactorCoefLeft[type];
    }
    
    virtual var_t* coef2(RootType jType, RootType kType) {
        return coefLeft[jType][kType];
    }
    
    virtual var_t halfIProductSingular1(const BetheSolution& sol) {
        // c * (i + eps + c * eps^L)^L
        return sol.c() * std::pow(1il + sol.epsilon() + sol.c() * std::pow(sol.epsilon(), _L), _L);
    }

    virtual var_t halfIProductSingular2(const BetheSolution& sol) {
        // 2i + c * eps^L
        return 2il + sol.c() * std::pow(sol.epsilon(), _L);
    }
    
    virtual var_t halfIProductDc1(const BetheSolution& sol) {
        // (i + eps + c * eps^L)^(L-1) * (i + eps + (L+1) * c * eps^L)
        var_t cEps = sol.c() * std::pow(sol.epsilon(), _L);
        return std::pow(1il + sol.epsilon() + cEps, _L-1) * (1il + sol.epsilon() + cEps * (elem_t)(_L+1));
    }
    
    virtual var_t halfIProductDc2(const BetheSolution& sol) {
        // eps^L
        return std::pow(sol.epsilon(), _L);
    }
    
    virtual var_t halfIProductDeps1(const BetheSolution& sol) {
        // (i + eps + c*eps^L)^(L-1)*c*L*(1 + c*L*eps^(L-1))
        var_t cEps2 = sol.c() * std::pow(sol.epsilon(), _L - 1);
        return std::pow(1il + sol.epsilon() * (1.0L + cEps2), _L - 1) * sol.c() * (elem_t)_L * (1.0L + (elem_t)_L * cEps2);
    }
    
    virtual var_t halfIProductDeps2(const BetheSolution& sol) {
        // c*L*eps^(L-1)
        return sol.c() * (elem_t)_L * std::pow(sol.epsilon(), _L - 1);
    }
    
    virtual var_t momentumHalfIFactorSingular(const BetheSolution& sol) {
        // i + eps + c * eps^L
        return 1.0il + sol.epsilon() + sol.c() * std::pow(sol.epsilon(), _L);
    }
    
    virtual var_t halfIProductMomentumDeps(const BetheSolution& sol) {
        // 1 + c * L * epsilon^(L-1) 
        return 1.0L + sol.c() * (elem_t)_L * std::pow(sol.epsilon(), _L - 1);
    }
    
    virtual var_t halfIProductMomentumDc(const BetheSolution& sol) {
        // epsilon^L 
        return std::pow(sol.epsilon(), _L);
    }
};

class BetheFactorRight : public BetheFactorBase
{
private:
public:
    BetheFactorRight(int L, int M) : BetheFactorBase(L, M) {
    }
    
    virtual var_t identicalRootFactor() {
        return 1.0L;
    }
    
    virtual var_t identicalRootDerivative() {
        return 0.0L;
    }
    
    virtual var_t* coef1(RootType type) {
        return prefactorCoefRight[type];
    }
    
    virtual var_t* coef2(RootType jType, RootType kType) {
        return coefRight[jType][kType];
    }
    
    virtual var_t halfIProductSingular1(const BetheSolution& sol) {
        // (1 + c eps^(L-1))^L * (2i + c * eps^L)
        return std::pow(1.0L + sol.c() * std::pow(sol.epsilon(), _L - 1), _L) * (2il + sol.c() * std::pow(sol.epsilon(), _L));
    }
    
    virtual var_t halfIProductSingular2(const BetheSolution& sol) {
        // c*(eps - i)^L
        return sol.c() * std::pow(sol.epsilon() - 1il, _L);
    }
    
    virtual var_t halfIProductDc1(const BetheSolution& sol) {
        // (eps + c * eps^L)^(L-1) * (L + eps + c * eps^L)
        var_t val = sol.epsilon() + sol.c() * std::pow(sol.epsilon(), _L);
        return std::pow(val, _L - 1) * ((elem_t)_L + val);
    }
    
    virtual var_t halfIProductDc2(const BetheSolution& sol) {
        // (eps - i)^L
        return std::pow(sol.epsilon() - 1il, _L);
    }
    
    virtual var_t halfIProductDeps1(const BetheSolution& sol) {
        // c*L*(eps + c*eps^L)^(L-2) * (L-1+eps+c*eps^L) * (1 + c*eps^(L-1))
        var_t val = 1.0L + sol.c() * std::pow(sol.epsilon(), _L-1);
        return sol.c() * (elem_t)_L * std::pow(val * sol.epsilon(), (elem_t)_L - 2.0L) * ((elem_t)_L - 1.0L + val * sol.epsilon()) * val;
    }
    
    virtual var_t halfIProductDeps2(const BetheSolution& sol) {
        // c*L*(eps-i)^(L-1)
        return std::pow(sol.epsilon() - 1il, _L-1) * sol.c() * (elem_t)_L;
    }
    
    virtual var_t momentumHalfIFactorSingular(const BetheSolution& sol) {
        // (1 + c*eps^(L-1))(-i+eps)
        return (1.0L + sol.c() * std::pow(sol.epsilon(), _L - 1)) * (-1.0il + sol.epsilon());
    }
    
    virtual var_t halfIProductMomentumDeps(const BetheSolution& sol) {
        // (1 + c*eps^(L-1)) + c * (L-1) eps^(L-2) (-i + eps) 
        // = 1 - i * c * (L-1) * eps^(L-2) + c * L * eps^(L-1)
        // =1 + c * L * eps^(L-2) * (-i(1 - 1/L) + eps)
        return 1.0L + sol.c() * (elem_t)_L * std::pow(sol.epsilon(), _L - 2) * (var_t(0, -1.0L + 1.0L/_L) + sol.epsilon());
    }
    
    virtual var_t halfIProductMomentumDc(const BetheSolution& sol) {
        // eps^(L-1) * (-i + eps) 
        return std::pow(sol.epsilon(), _L - 1) * (-1il + sol.epsilon());
    }
};
