#include "Homotopy.h"

Vector SimpleHomotopy::eval(elem_t t, const Solution &sol) {
    return _start->eval(sol) * (_gamma * t) + _target->eval(sol) * (1-t);
}

Matrix SimpleHomotopy::diff(elem_t t, const Solution &sol) {
    return _start->diff(sol) * (_gamma * t) + _target->diff(sol) * (1-t);
}

Vector SimpleHomotopy::impliciteDiff(elem_t t, const Solution &sol) {
    Matrix m = _target->diff(sol) * (1.0 - t) + _start->diff(sol) * (_gamma * t);
    return  m.inverse() * (_target->eval(sol) - _start->eval(sol) * _gamma);
}

