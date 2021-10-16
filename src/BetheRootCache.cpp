#include "BetheRootCache.h"
#include <fstream>

void BetheRootCache::load() {
    std::string file = rootFolder + "BetheRootsL" + ToString(_L) + "M" + ToString(_M) + ".txt";
    std::ifstream in(file.c_str());
    int cnt = 0;
    std::vector<elem_t> real(_M);
    std::vector<elem_t> imag(_M);
    std::string tmp;
    while (in >> tmp) {
        for (int i = 0; i < _M; i++) {
            in >> real[i];
        }
        for (int i = 0; i < _M; i++) {
            in >> imag[i];
        }
        Vector root(_M);
        bool singular = false;
        int pos = 0;
        for (int i = 0; i < _M; i++) {
            if (isSingular(real[i], imag[i])) {
                singular = true;
//                continue;
            }
            root(pos++) = var_t(real[i], imag[i]);
        }
        if (singular && !isPhysical(root)) {
            continue;
        }
//        if (momentumConserved(root)) {
            _roots.push_back(root);
            _isSingular.push_back(singular ? 1 : 0);
//        }
    }

    in.close();
}

bool BetheRootCache::isPhysical(Vector& roots) {
    var_t ip = (_L % 2 == 0 ? 1.0L : -1.0L);
    for (int i = 0; i < roots.size(); i++) {
        if (isSingular(roots(i).real(), roots(i).imag())) continue;
        ip *= std::pow(e2ip(roots(i)), _L);
    }
    
    return isZero(ip - 1.0L);
}

bool BetheRootCache::momentumConserved(Vector& sol) {
    var_t res(1.0);
    for (int i = 0; i < sol.size(); i++) {
        if (abs(sol[i].real()) < 1e-6 && (abs(sol[i].imag()-0.5) < 1e-6 || abs(sol[i].imag()+0.5) < 1e-6)) {
            return false;
        }
        res *= e2ip(sol[i]);
    }
    return true;
//    return std::abs(res - var_t(1.0)) < 1e-6;
}

bool BetheRootCache::isSingular(elem_t real, elem_t imag) {
    return abs(real) < EPS && (abs(imag + 0.5) < EPS || abs(imag - 0.5) < EPS);
}
