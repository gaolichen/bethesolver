#pragma once
#include "common.h"
#include "Function.h"


extern var_t prefactorCoefLeftA[2][2];

extern var_t prefactorCoefRightA[2][2];

extern var_t coefLeftA[2][2][4];

extern var_t coefRightA[2][2][4];


class BethePolynomialBase : public Function
{
protected:
    int _L;
    int _M;
public:
    BethePolynomialBase(int L, int M);
    
    int chainLength() { return _L; }
    
    int magnonNumber() { return _M; }
    
    virtual var_t* coef1(bool inverted) = 0;
    virtual var_t* coef2(bool jInverted, bool kInverted) = 0;
        
    virtual var_t halfIProductSingular1(const InvertableSolution& sol) = 0;
    virtual var_t halfIProductSingular2(const InvertableSolution& sol) = 0;
    virtual var_t momentumHalfIFactorSingular(const InvertableSolution& sol) = 0;
    
    virtual var_t halfIProductDc1(const InvertableSolution& sol) = 0;
    virtual var_t halfIProductDc2(const InvertableSolution& sol) = 0;
    virtual var_t halfIProductMomentumDc(const InvertableSolution& sol) = 0;

    virtual var_t halfIProductDeps1(const InvertableSolution& sol) = 0;
    virtual var_t halfIProductDeps2(const InvertableSolution& sol) = 0;
    virtual var_t halfIProductMomentumDeps(const InvertableSolution& sol) = 0;

    
    // returns (uj +- i/2)^L
    var_t halfIProduct(const InvertableSolution& sol, int j) {
        if (j == sol.indexOfSingularRoot1()) {
            return halfIProductSingular1(sol);
        } else if (j == sol.indexOfSingularRoot2()) {
            return halfIProductSingular2(sol);
        } else {
            return std::pow(halfIFactor(sol, j), _L);
        }
    }
    
    // returns (uj +- i/2)
    var_t halfIFactor(const InvertableSolution& sol, int j) {
        var_t* c = coef1(sol.inverted(j));
        return c[0] + c[1] * sol.root(j);
    }
    
    // derivative of (uj +- i/2) with respect to uj
    var_t halfIFactorDiffJ(const InvertableSolution& sol, int j) {
        var_t* c = coef1(sol.inverted(j));
        return c[1];
    }
    
    // product of (uj - uk +- i) for all k != j and k != toSkip
    // if uj is singular, only nonsigular uk is included. 
    var_t factorProduct(const InvertableSolution& sol, int j, int toSkip = -1) {
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
            ret *= factorJK(uj, sol.get2(k), sol.inverted(j), sol.inverted(k));
        }
        return ret;
    }
    
    // factor uj - uk +- i
    var_t factorJK(const InvertableSolution& sol, int j, int k) {
        return factorJK(sol.get2(j), sol.get2(k), sol.inverted(j), sol.inverted(k));
    }
    
    // factor uj - uk +- i
    var_t factorJK(var_t uj, var_t uk, bool invertedJ, bool invertedK) {
//        var_t vec[4] = {1.0L, uj, uk, uj * uk};
        var_t *c = coef2(invertedJ, invertedK);
        return 1.0L * c[0] + uj * c[1] + uk * c[2] + c[3] * uj * uk;
    }
    
    // derivative of the (uj-uk +- i) factor with respect to uk.
    var_t diffK(const InvertableSolution& sol, int j, int k) {
        var_t *c = coef2(sol.inverted(j), sol.inverted(k));
        return c[2] + c[3] * sol.get2(j);
    }
    
    // derivative of the (uj-uk +- i) factor with respect to uj.
    var_t diffJ(const InvertableSolution& sol, int j, int k) {
        var_t *c = coef2(sol.inverted(j), sol.inverted(k));
        return c[1] + c[3] * sol.get2(k);
    }

    virtual Vector eval(const Solution &sol_in) {
        const InvertableSolution& sol = dynamic_cast<const InvertableSolution &>(sol_in);
        Vector ret(sol.size());
        for (int j = 0; j < sol.size(); j++) {
            ret[j] = halfIProduct(sol, j) * factorProduct(sol, j);
        }
        return ret;
    }
    
    var_t evalMomentum(const Solution &sol_in) {
        const InvertableSolution& sol = dynamic_cast<const InvertableSolution &>(sol_in);
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
        const InvertableSolution& sol = dynamic_cast<const InvertableSolution &>(sol_in);
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
        
        return ret;
    }
    
    virtual Matrix diff(const Solution &sol_in) {
        const InvertableSolution& sol = dynamic_cast<const InvertableSolution &>(sol_in);
        Matrix ret(sol.size(), sol.size());
        
        // step 1: compute off-diagonal terms of entry (j,k) execept for:
        // 1. (j,k) are a pair of singular roots
        for (int j = 0; j < sol.size(); j++) {
            var_t prod = factorProduct(sol, j);
            var_t prefactor = halfIProduct(sol, j);
            for (int k = 0; k < sol.size(); k++) {
//                ret(j, k) = 0.0L;
                if (j == k) continue;
                
                if (sol.isSingular(j) && sol.isSingular(k)) continue;
                
                if (!isZero(prod)) {
                    ret(j, k) = diffK(sol, j, k) / factorJK(sol, j, k);
                } else {
                    ret(j, k) = diffK(sol, j, k) * factorProduct(sol, j, k);
                }
                
                if (!isZero(prod)) {
                    ret(j, k) *= prod;
                }                
                ret(j, k) *= prefactor;                
            }
        }
        
        // step 2: compute diagonal terms except for singular roots.
        for (int j = 0; j < sol.size(); j++) {
            if (sol.isSingular(j)) continue;
            var_t prefactor = halfIProduct(sol, j);
            var_t prod = factorProduct(sol, j);
            ret(j, j) = 0.0L;
            
            // contribute from uj-uk +- i terms
            for (int k = 0; k < sol.size(); k++) {
                if (k == j) continue;
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


class BethePolynomialLeft : public BethePolynomialBase
{
private:
public:
    BethePolynomialLeft(int L, int M) : BethePolynomialBase(L, M) {
    }
        
    virtual var_t* coef1(bool inverted) {
        return prefactorCoefLeftA[inverted];
    }
    
    virtual var_t* coef2(bool jInverted, bool kInverted) {
        return coefLeftA[jInverted][kInverted];
    }
    
    virtual var_t halfIProductSingular1(const InvertableSolution& sol) {
        // c * (i + eps + c * eps^L)^L
        return sol.c() * std::pow(1il + sol.epsilon() + sol.c() * std::pow(sol.epsilon(), _L), _L);
    }

    virtual var_t halfIProductSingular2(const InvertableSolution& sol) {
        // 2i + c * eps^L
        return 2il + sol.c() * std::pow(sol.epsilon(), _L);
    }
    
    virtual var_t halfIProductDc1(const InvertableSolution& sol) {
        // (i + eps + c * eps^L)^(L-1) * (i + eps + (L+1) * c * eps^L)
        var_t cEps = sol.c() * std::pow(sol.epsilon(), _L);
        return std::pow(1il + sol.epsilon() + cEps, _L-1) * (1il + sol.epsilon() + cEps * (elem_t)(_L+1));
    }
    
    virtual var_t halfIProductDc2(const InvertableSolution& sol) {
        // eps^L
        return std::pow(sol.epsilon(), _L);
    }
    
    virtual var_t halfIProductDeps1(const InvertableSolution& sol) {
        // (i + eps + c*eps^L)^(L-1)*c*L*(1 + c*L*eps^(L-1))
        var_t cEps2 = sol.c() * std::pow(sol.epsilon(), _L - 1);
        return std::pow(1il + sol.epsilon() * (1.0L + cEps2), _L - 1) * sol.c() * (elem_t)_L * (1.0L + (elem_t)_L * cEps2);
    }
    
    virtual var_t halfIProductDeps2(const InvertableSolution& sol) {
        // c*L*eps^(L-1)
        return sol.c() * (elem_t)_L * std::pow(sol.epsilon(), _L - 1);
    }
    
    virtual var_t momentumHalfIFactorSingular(const InvertableSolution& sol) {
        // i + eps + c * eps^L
        return 1.0il + sol.epsilon() + sol.c() * std::pow(sol.epsilon(), _L);
    }
    
    virtual var_t halfIProductMomentumDeps(const InvertableSolution& sol) {
        // 1 + c * L * epsilon^(L-1) 
        return 1.0L + sol.c() * (elem_t)_L * std::pow(sol.epsilon(), _L - 1);
    }
    
    virtual var_t halfIProductMomentumDc(const InvertableSolution& sol) {
        // epsilon^L 
        return std::pow(sol.epsilon(), _L);
    }
};

class BethePolynomialRight : public BethePolynomialBase
{
private:
public:
    BethePolynomialRight(int L, int M) : BethePolynomialBase(L, M) {
    }
        
    virtual var_t* coef1(bool inverted) {
        return prefactorCoefRightA[inverted];
    }
    
    virtual var_t* coef2(bool jInverted, bool kInverted) {
        return coefRightA[jInverted][kInverted];
    }
    
    virtual var_t halfIProductSingular1(const InvertableSolution& sol) {
        // (1 + c eps^(L-1))^L * (2i + c * eps^L)
        return std::pow(1.0L + sol.c() * std::pow(sol.epsilon(), _L - 1), _L) * (2il + sol.c() * std::pow(sol.epsilon(), _L));
    }
    
    virtual var_t halfIProductSingular2(const InvertableSolution& sol) {
        // c*(eps - i)^L
        return sol.c() * std::pow(sol.epsilon() - 1il, _L);
    }
    
    virtual var_t halfIProductDc1(const InvertableSolution& sol) {
        // (eps + c * eps^L)^(L-1) * (L + eps + c * eps^L)
        var_t val = sol.epsilon() + sol.c() * std::pow(sol.epsilon(), _L);
        return std::pow(val, _L - 1) * ((elem_t)_L + val);
    }
    
    virtual var_t halfIProductDc2(const InvertableSolution& sol) {
        // (eps - i)^L
        return std::pow(sol.epsilon() - 1il, _L);
    }
    
    virtual var_t halfIProductDeps1(const InvertableSolution& sol) {
        // c*L*(eps + c*eps^L)^(L-2) * (L-1+eps+c*eps^L) * (1 + c*eps^(L-1))
        var_t val = 1.0L + sol.c() * std::pow(sol.epsilon(), _L-1);
        return sol.c() * (elem_t)_L * std::pow(val * sol.epsilon(), (elem_t)_L - 2.0L) * ((elem_t)_L - 1.0L + val * sol.epsilon()) * val;
    }
    
    virtual var_t halfIProductDeps2(const InvertableSolution& sol) {
        // c*L*(eps-i)^(L-1)
        return std::pow(sol.epsilon() - 1il, _L-1) * sol.c() * (elem_t)_L;
    }
    
    virtual var_t momentumHalfIFactorSingular(const InvertableSolution& sol) {
        // (1 + c*eps^(L-1))(-i+eps)
        return (1.0L + sol.c() * std::pow(sol.epsilon(), _L - 1)) * (-1.0il + sol.epsilon());
    }
    
    virtual var_t halfIProductMomentumDeps(const InvertableSolution& sol) {
        // (1 + c*eps^(L-1)) + c * (L-1) eps^(L-2) (-i + eps) 
        // = 1 - i * c * (L-1) * eps^(L-2) + c * L * eps^(L-1)
        // =1 + c * L * eps^(L-2) * (-i(1 - 1/L) + eps)
        return 1.0L + sol.c() * (elem_t)_L * std::pow(sol.epsilon(), _L - 2) * (var_t(0, -1.0L + 1.0L/_L) + sol.epsilon());
    }
    
    virtual var_t halfIProductMomentumDc(const InvertableSolution& sol) {
        // eps^(L-1) * (-i + eps) 
        return std::pow(sol.epsilon(), _L - 1) * (-1il + sol.epsilon());
    }
};
