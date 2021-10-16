#pragma once
#include <Eigen/Dense>
#include "common.h"

class BetheRootCache
{
private:
    std::string rootFolder = "/home/gaolichen/gitroot/mywork/chaos/BetheRoots/";
    int _L;
    int _M;
    std::vector<Vector> _roots;
    std::vector<int> _isSingular;
    void load();
    static bool momentumConserved(Vector& sol);
    static bool isSingular(elem_t real, elem_t imag);
    bool isPhysical(Vector& roots);
public:
    BetheRootCache(int L, int M) {
        _L = L;
        _M = M;
        load();
    }
    
    int chainLength() { return _L; }
    
    int magnonNumber() { return _M; }
    
    int numberOfRoots() {
        return _roots.size();
    }
    
    Vector getRoot(int index) {
        return _roots[index];
    }
    
    static var_t e2ipTotal(Vector &root) {
        var_t res(1.0);
        for (int i = 0; i < root.size(); i++) {
            if (isSingular(root[i].real(), root[i].imag())) {
                return 0.0;
            }
            res *= e2ip(root[i]);
        }
        return res;
    }
    
    var_t momentum(int index) {
        var_t val = e2ipTotal(_roots[index]);
        return log(val) * var_t(0.0, -1.0 / Pi);
    }
    
    bool isSingular(int index) {
        return _isSingular[index] > 0;
    }
};
