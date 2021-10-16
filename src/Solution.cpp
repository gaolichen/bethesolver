#include "Solution.h"

elem_t eps_may_be_singular = 0.1L;
elem_t blowup_root_cutoff = 1.5L;

bool BetheSolution::update(elem_t beta, RootFitter* fitter) {
    bool ret1 = setupBlowupRoots(beta);
    bool ret2 = setupSingularRoots(beta);
    if (ret1 || ret2) {
        if (fitter != NULL) {
            fitter->clear();
            LOG("t=" << beta << " RootFitter.clear()", Info);
        }
    }
    return ret1 || ret2;
}

bool mayBeSingular(var_t u) {
        return std::abs(u.real()) < eps_may_be_singular &&
        (std::abs(u.imag() + 0.5) < eps_may_be_singular || std::abs(u.imag() - 0.5) < eps_may_be_singular);
}

bool cannotBeSingular(var_t& epsilon) {
    return std::abs(epsilon.real()) > 2 * eps_may_be_singular || std::abs(epsilon.imag()) > 2 * eps_may_be_singular;
}

bool SingularRootSolution::setupSingularRoots(elem_t beta) {
    if (_indexOfSingularRoot1 >= 0 && _indexOfSingularRoot2 >= 0) {
        if (cannotBeSingular(this->root(_indexOfSingularRoot2))) {
            // the the solution is not singular any more.
            var_t v1 = this->get(_indexOfSingularRoot1);
            var_t v2 = this->get(_indexOfSingularRoot2);
            this->root(_indexOfSingularRoot2) = v2;
            this->root(_indexOfSingularRoot1) = v1;
            _indexOfSingularRoot2 = _indexOfSingularRoot1 = -1;
            LOG("t=" << beta << ", removed singular solutions: " << root(), Info);
            return true;
        }
        return false;
    }

    std::vector<int> list1, list2;
    for (int i = 0; i < size(); i++) {
        if (isRegular(i) && mayBeSingular(this->root(i))) {
            if (root(i).imag() < 0.0) {
                // this is the root -i/2
                list2.push_back(i);
            } else {
                // it is the root i/2
                list1.push_back(i);
            }
        }
    }
        
    if (list1.size() <= 0 || list2.size() <= 0) {
        return false;
    }
    
    int index1 = -1;
    int index2 = -1;
    bool found = false;
    for (int i = 0; i < list1.size(); i++) {
        if (found) break;
        for (int j = 0; j < list2.size(); j++) {
//            if ((root(list1[i]).imag() - 0.5) * (root(list2[j]).imag() + 0.5) < 0.0) continue;
            var_t v1 = root(list1[i]) + root(list2[j]);
            var_t v2 = root(list1[i]) - root(list2[j]);
            if (std::abs(v1.imag()) < 1e-6 && std::abs(v2.real()) < 1e-6) {
                index1 = list1[i];
                index2 = list2[j];
                found = true;
                break;
            }
        }
    }
    if (!found) return false;
    
//    std::cout <<"t=" << beta << ", index1=" << index1 << ", index2=" << index2 << std::endl;

    // comment out the following part, because the condition 
    // fails some cases that should be treated as singular roots.
    // check singular consistence condition.
    // (-1)^L * \prod_{k} (uk + i/2)^L/(uk-i/2)^L = exp(2i*L*(M+2)beta)
/*    var_t prod = (_L % 2 == 0 ? 1.0L : -1.0L);
    for (int k = 0; k < size(); k++) {
        if (k == index1 || k == index2) continue;
        prod *= std::pow((get(k) + 0.5il)/(get(k) - 0.5il), _L);
    }
    
    elem_t theta = std::arg(prod / std::polar(1.0L, 2.0L * _L * beta * (size() + 2)));
    // here looks strange, but if any of the condition is removed, some unit tests fails.
    if (std::abs(theta) > 0.1 && std::abs(theta)/(2 * _L * (size() + 2)) > 0.1) {
        std::cout << "prod=" << prod / std::polar(1.0L, 2.0L * _L * beta * (size() + 2)) << ", theta=" << std::abs(theta)/(2 * _L * (size() + 2)) << std::endl;
        return;
    }*/
    
    // found a pair of singular solutions.
    this->root(index2) += 0.5il;
    this->root(index1) -= 0.5il;
    this->root(index1) -= this->root(index2);
    // find c, which is defined as
    // u_1 = i/2 + eps + c * eps^N
    // u_2 = -i/2 + eps
    if (!isZero(this->root(index2)) && !isZero(this->root(index1))) {
        this->root(index1) /= std::pow(this->root(index2), _L);
        //std::cout << "this->root(index1)=" << this->root(index1) << std::endl;
    } else {
        this->root(index1) = std::polar(2.0L, 2.0L * _L * beta - Pi * 0.5 * (_L - 1));
        for (int k = 0; k < size(); k++) {
            if (k == index1 || k == index2) continue;
            this->root(index1) *= (this->get(k) - 1.5il);
            this->root(index1) /= (this->get(k) + 0.5il);
        }
    }
    _indexOfSingularRoot1 = index1;
    _indexOfSingularRoot2 = index2;
    LOG("t=" << beta<< ", find singular solutions: eps=" << epsilon() << ", c=" << c() << " ,root=" << root(), Info);
//    std::cout << "normalRoot()=" << normalRoot() << std::endl;
    return true;
}

bool BetheSolution::setupBlowupRoots(elem_t beta) {
    // 1. check whether or not _indexOfBlowupRoot need to change
    //       if _indexOfBlowupRoot >= 0, only select new from large root
    //       else select from regular roots.
    // 2. if need change, make change 
    // 3.
    int indexOfMaxRoot = -1;
    elem_t maxMagnitude = 0.0L;
    if (_indexOfBlowupRoot == -1) {
        for (int i = 0; i < this->size(); i++) {
            if (isSingular(i)) continue;
            if (indexOfMaxRoot == -1 || (magnitude(root(i)) > maxMagnitude)) {
                indexOfMaxRoot = i;
                maxMagnitude = magnitude(root(i));
            }
        }
//        std::cout << "maxMagnitude=" << maxMagnitude << std::endl;
        if (maxMagnitude < blowup_root_cutoff) {
            return false;
        }
//        std::cout << "beta=" << beta << ", root before: " << root() << std::endl;
        for (int i = 0; i < this->size(); i++) {
            if (i == indexOfMaxRoot || isSingular(i)) continue;
            if (magnitude(root(i)) > blowup_root_cutoff) {
                root(i) /= root(indexOfMaxRoot);
                _types[i] = large;
            }
        }
        root(indexOfMaxRoot) = 1.0L/root(indexOfMaxRoot);
        _types[indexOfMaxRoot] = largest;
        _indexOfBlowupRoot = indexOfMaxRoot;
        LOG("beta=" << beta  << "found blowup roots: " << root(), Info);
        
        return true;
    } else if (magnitude(omega()) > 1.0L/blowup_root_cutoff) {
        // remove largest root
        for (int i = 0; i < size(); i++) {
            if (type(i) == large) {
                _types[i] = regular;
                root(i) /= omega();
                _identical[i] = 0;
            }
        }
        root(_indexOfBlowupRoot) = 1.0L / root(_indexOfBlowupRoot);
        _types[_indexOfBlowupRoot] = regular;
        _indexOfBlowupRoot = -1;
        LOG("beta=" << beta << ", removed largest roots, root()=" << root(), Info);
        return true;
    } else {
        int indexOfMaxRoot = -1;
        for (int i = 0; i < this->size(); i++) {
            if (type(i) != large) continue;
            if (indexOfMaxRoot == -1 || magnitude(root(i)) > magnitude(root(indexOfMaxRoot)) + EPS) {
                indexOfMaxRoot = i;
            }
        }
        if (indexOfMaxRoot >= 0 && std::abs(root(indexOfMaxRoot)) > blowup_root_cutoff) {
            // _indexOfBlowupRoot need to be changed.
//            std::cout << "beta=" << beta << ", indexOfMaxRoot=" << indexOfMaxRoot <<", before roots changed: " << root() << std::endl;
            for (int i = 0; i < this->size(); i++) {
                if (i == _indexOfBlowupRoot || i == indexOfMaxRoot || _types[i] != large) continue;
                root(i) /= root(indexOfMaxRoot);
            }
            
            _types[indexOfMaxRoot] = largest;
            _types[_indexOfBlowupRoot] = large;
            root(indexOfMaxRoot) = root(_indexOfBlowupRoot) / root(indexOfMaxRoot);
            root(_indexOfBlowupRoot) = root(indexOfMaxRoot) / root(_indexOfBlowupRoot);
            _indexOfBlowupRoot = indexOfMaxRoot;
            
            LOG("beta=" << beta << "blowup roots updated: " << root(), Info);
            // for simplicity, just return and do not do others
            return true;
        } else {
            // if _indexOfBlowupRoot >=0 and is not changed, we update large roots
            elem_t maxAbs = 1.0L/std::abs(root(_indexOfBlowupRoot));
            bool updated = false;
            for (int i = 0; i < size(); i++) {
                if (_types[i] == large) {
                    if (magnitude(root(i)) < 0.7L) {
                        // this root is too small, we need to change to regular root
                        root(i) /= root(_indexOfBlowupRoot);
                        _types[i] = regular;
                        _identical[i] = 0;
                        updated = true;
                    } else if (std::abs(root(i) - 1.0L) < EPS){
                        // this is the identical infinity roots.
                        // if omega() is very small
                        if (std::abs(omega()) < 1e-2) {
                            if (!_identical[i]) {
                                LOG("beta=" << beta << ", identical root:" << i, Info);
                                _identical[i] = 1;
//                                root(i) = 1.0L;
                                updated = true;
                            }
                        } else if (_identical[i]) {
                            // if omega() is large.
                            _identical[i] = 0;
                            updated = true;
                        }
                    }
                } else if (_types[i] == regular && !isSingular(i) && std::abs(root(i)) > std::max(blowup_root_cutoff, maxAbs * 0.9)) {
                    // the root is large enough to be a large root.
                    _types[i] = large;
                    root(i) *= root(_indexOfBlowupRoot);
                    updated = true;
                } 
            }
            return updated;
        }
    }
}

