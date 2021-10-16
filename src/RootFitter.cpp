#include "RootFitter.h"

void RootFitter::clear() {
    _begin = 0;
    _end = 0;
}

void RootFitter::addLastRoot(const Vector& root, elem_t beta) {
    _lastRoots[_end] = root;
    _beta[_end] = beta;
    _end = NEXT(_end);
    if (_end == _begin) {
        _begin = NEXT(_begin);
    }
}

void RootFitter::invert(int i) {
    for (int j = _begin; j != _end; j = NEXT(j)) {
        _lastRoots[j](i) = 1.0L / _lastRoots[j](i);
    }
}

Vector RootFitter::fit(elem_t beta) const {
    if (_end == _begin || NEXT(_begin) == _end) {
        // if there is no root or only one root, return empty.
        return Vector(0);
    }
    
    int i1 = _begin;
    int i2 = NEXT(i1);
    int i3 = NEXT(i2);
    int i4 = NEXT(i3);

    if (i3 == _end) {
        // only two roots cached.
        elem_t ratio = (beta - _beta[i2]) / (_beta[i2] - _beta[i1]);
        return _lastRoots[i2] + (_lastRoots[i2] - _lastRoots[i1]) * ratio;
    } else if (i4 == _end) {
        // three roots cached.
        elem_t d21 = _beta[i2] - _beta[i1];
        elem_t d31 = _beta[i3] - _beta[i1];
        elem_t d = beta - _beta[i1];
        Eigen::Matrix<elem_t, 2, 2> mat;
        mat <<  d21 * d21, d21,
                d31 * d31, d31;
        Eigen::Matrix<elem_t, 1, 2> vec1(d * d, d);
        Eigen::Matrix<elem_t, 1, 2> res = vec1 * mat.inverse();
        return _lastRoots[i1] + (_lastRoots[i2] - _lastRoots[i1]) * res(0) + (_lastRoots[i3] - _lastRoots[i1]) * res(1);
    } else {
        // four changes
        elem_t d21 = _beta[i2] - _beta[i1];
        elem_t d31 = _beta[i3] - _beta[i1];
        elem_t d41 = _beta[i3] - _beta[i1];
        elem_t d = beta - _beta[i1];
        Eigen::Matrix<elem_t, 3, 3> mat;
        mat <<  d21 * d21 * d21, d21 * d21, d21,
                d31 * d31 * d31, d31 * d31, d31,
                d41 * d41 * d41, d41 * d41, d31;
        Eigen::Matrix<elem_t, 1, 3> vec1(d * d * d, d * d, d);
        Eigen::Matrix<elem_t, 1, 3> res = vec1 * mat.inverse();
        return _lastRoots[i1] + (_lastRoots[i2] - _lastRoots[i1]) * res(0)
            + (_lastRoots[i3] - _lastRoots[i1]) * res(1) + (_lastRoots[i4] - _lastRoots[i1]) * res(2);
    }
}

