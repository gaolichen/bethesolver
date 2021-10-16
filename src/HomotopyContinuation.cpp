#include "HomotopyContinuation.h"

HomotopyContinuation::HomotopyContinuation() {
}

Vector HomotopyContinuation::newtonRaphson(Solution &sol, Homotopy& homotopy, elem_t t) {
    Matrix mDiff = homotopy.diff(t, sol);
    Vector vEval = homotopy.eval(t, sol);
    var_t det = mDiff.determinant();
    if (std::abs(det) < std::pow(EPS, mDiff.rows())) {
        return Vector(0);
    }
    return -mDiff.inverse() * vEval;
}

void SimpleHomotopyContinuation::solve(Homotopy &homotopy, Solution& sol) {
    Vector nr = sol.normalRoot();
    homotopy.setStartRoot(nr);
    elem_t delta = (homotopy.tEnd() - homotopy.tBegin()) / homotopy.getSteps();
    
    for (int i = 0; i < homotopy.getSteps(); i++) {
        elem_t t = homotopy.tBegin() + delta * i;
        Vector old = sol.root();
        trace(t, &sol, homotopy.error(t, sol));
        elem_t error0 = homotopy.error(t + delta, sol);
        Vector k1 = homotopy.impliciteDiff(t, sol);
        sol.root() += k1 * (delta/2);
        Vector k2 = homotopy.impliciteDiff(t + delta/2, sol);
        sol.setRoot(old + k2 * (delta/2));
        Vector k3 = homotopy.impliciteDiff(t + delta/2, sol);
        sol.setRoot(old + k3 * (delta/2));
        Vector k4 = homotopy.impliciteDiff(t + delta, sol);
        Vector change = k1 * (delta/6) + k2 * (delta/3) + k3 * (delta/3) + k4 * (delta/6);
        sol.setRoot(old + change);
//        elem_t error = homotopy.error(t + delta, sol);
        
        t += delta;
        for (int j = 0; j < 20; j++) {
            Vector correction = newtonRaphson(sol, homotopy, t);
            sol.root() += correction;
            if (correction.norm() < accuracyLow()) {
                break;
            }
//            error = homotopy.error(t, sol);
        }
        sol.update(t);
    }
    
    for (int j = 0; j < 50; j++) {
            Vector correction = newtonRaphson(sol, homotopy, homotopy.tEnd());
            sol.root() += correction;
            if (correction.norm() < accuracyHigh()) {
                break;
            }
    }
}

/*
void BetheHomotopyContinuation::solve(Homotopy &homotopy, Solution& sol) {
    elem_t maxDelta = (homotopy.tEnd() - homotopy.tBegin()) / homotopy.getSteps();
    elem_t minDelta = maxDelta / 2.0;
    elem_t delta = maxDelta;
    Vector nr = sol.normalRoot();
    homotopy.setStartRoot(nr);
    
    elem_t error, totalError, lastError = 1.0e20;
    elem_t t = homotopy.tBegin();
    int steps = 0;
    int foundSingular = 0;
    bool hasChange;
    
    elem_t maxChange = 0.0L;
    while (std::abs(t - homotopy.tEnd()) > EPS) {
        hasChange = false;
        Vector old = sol.root();
        trace(t, &sol, homotopy.error(t, sol));
        elem_t error0 = homotopy.error(t + delta, sol);
        Vector k1 = homotopy.impliciteDiff(t, sol);
        sol.root() += k1 * (delta/2);
        Vector k2 = homotopy.impliciteDiff(t + delta/2, sol);
        sol.setRoot(old + k2 * (delta/2));
        Vector k3 = homotopy.impliciteDiff(t + delta/2, sol);
        sol.setRoot(old + k3 * (delta/2));
        Vector k4 = homotopy.impliciteDiff(t + delta, sol);
        Vector change = k1 * (delta/6) + k2 * (delta/3) + k3 * (delta/3) + k4 * (delta/6);
        
        elem_t rootNorm = old.norm();
        if (rootNorm > 0.1 && change.norm() > rootNorm * 0.1) {
            LOG("t=" << t << ", old.norm()=" << rootNorm << " change.norm()=" << change.norm(), Info);
            LOG("old=" << std::setprecision(18) << old << ", change=" << change, Info);
            error = error0;
            sol.setRoot(old);
        } else {
            sol.setRoot(old + change);
            error = homotopy.error(t + delta, sol);
            hasChange = true;
        }

        t += delta;
        steps++;
                
        int times = 0;
        while (times < 20 && error > homotopy.tolerance()) {
            Vector correction = newtonRaphson(sol, homotopy, t);
            if (!homotopy.acceptCorrection(t, sol, correction)) {
                LOG("correction rejected: t=" << t << ", times=" << times, Warning);
                std::cout << "correction.norm()=" << correction.norm() << ", rootNorm=" << rootNorm << std::endl;
                break;
            }
            if (correction.norm() > rootNorm * 0.1) {
                break;
            }
            hasChange = true;
            sol.root() += correction;
            error = homotopy.error(t, sol);
            times++;
        }
        if (times == 20) {
            LOG("t=" << t << ", run newtonRaphson more than 20 times, error=" << error, Warning);
        }
        if (!hasChange) {
            LOG("process terminated at t=" << t, Error);
            break;
        }
        sol.update(t);
    }
    if (std::abs(t - homotopy.tEnd()) < EPS) {
        trace(homotopy.tEnd(), &sol, error);
    }
}
*/

BetheHomotopyContinuation::BetheHomotopyContinuation() {
    initLastChanges();
}

void BetheHomotopyContinuation::initLastChanges() {
    _lastChangeBegin = 0;
    _lastChangeEnd = 0;
}

void BetheHomotopyContinuation::addLastChange(Vector& change, elem_t delta) {
    _lastChanges[_lastChangeEnd] = change;
    _deltas[_lastChangeEnd] = delta;
    _lastChangeEnd = (_lastChangeEnd + 1) % 3;
    if (_lastChangeEnd == _lastChangeBegin) {
        _lastChangeBegin = (_lastChangeBegin + 1) % 3;
    }
}

/*Vector BetheHomotopyContinuation::predictChange(Homotopy &homotopy, Solution& sol, elem_t t) {
    if (_lastChangeEnd == _lastChangeBegin) {
        return Vector(0);
    }
    
    int i1 = _lastChangeBegin;
    int i2 = (_lastChangeBegin + 1) % 3;
    int i3 = (_lastChangeBegin + 2) % 3;

    if (i2 == _lastChangeEnd) {
        // only one change available.
        return _lastChanges[i1];
    } else {
        Vector ret = _lastChanges[i1];
        Vector old = sol.root(); 
        
        sol.root() += _lastChanges[i1];
        elem_t minError = homotopy.error(t, sol);
        
        sol.root() += (_lastChanges[i2] - _lastChanges[i1]) * 2.0L;
        elem_t error2 = homotopy.error(t, sol);
        if (error2 < minError) {
            minError = error2;
            ret = _lastChanges[i2] + _lastChanges[i2] - _lastChanges[i1];
        }
        if (i3 == _lastChangeEnd) {
            sol.root() = old;
            return ret;
            // two changes available.
//            return _lastChanges[i2] + _lastChanges[i2] - _lastChanges[i1];
        } else {
            sol.root() += _lastChanges[i3] * 3.0 - _lastChanges[i2] * 5.0L + _lastChanges[i1] * 2.0L;
            elem_t error3 = homotopy.error(t, sol);
            sol.root() = old;
            if (error3 < minError) {
                return (_lastChanges[i3] - _lastChanges[i2]) * 3.0L + _lastChanges[i1];
            } else {
                return ret;
            }
            // three changes
//            return (_lastChanges[i3] - _lastChanges[i2]) * 3.0L + _lastChanges[i1];
        }
    }
}*/

Vector BetheHomotopyContinuation::predictChange(elem_t delta) {
    if (_lastChangeEnd == _lastChangeBegin) {
        return Vector(0);
    }
    
    int i1 = _lastChangeBegin;
    int i2 = (_lastChangeBegin + 1) % 3;
    int i3 = (_lastChangeBegin + 2) % 3;

    if (i2 == _lastChangeEnd) {
        // only one change available.
        return _lastChanges[i1];
    } else if (i3 == _lastChangeEnd) {
        // two changes available.
        elem_t ratio = (delta + _deltas[i2]) / (_deltas[i2] + _deltas[i1]);
        return _lastChanges[i2] + (_lastChanges[i2] - _lastChanges[i1]) * ratio;
    } else {
        // three changes
        elem_t d1 = (_deltas[i2] + _deltas[i1]) * 0.5L;
        elem_t d12 = _deltas[i2] + (_deltas[i3] + _deltas[i1]) * 0.5L;
        elem_t d123 = _deltas[i2] + _deltas[i3] + (_deltas[i1] + delta) * 0.5L;
        Eigen::Matrix<elem_t, 2, 2> mat;
        mat <<  d1 * d1, d1,
                d12 * d12, d12;
        Eigen::Matrix<elem_t, 1, 2> vec1(d123 * d123, d123);
        Eigen::Matrix<elem_t, 1, 2> res = vec1 * mat.inverse();
        return _lastChanges[i1] + (_lastChanges[i2] - _lastChanges[i1]) * res(0) + (_lastChanges[i3] - _lastChanges[i1]) * res(1);
    }
}

Vector BetheHomotopyContinuation::nextChange(Homotopy &homotopy, const Solution &sol_in, elem_t t, elem_t delta) {
        Solution* sol = sol_in.copy();
        
        Vector k1 = homotopy.impliciteDiff(t, *sol);
        sol->root() += k1 * (delta/2);
        Vector k2 = homotopy.impliciteDiff(t + delta/2, *sol);
        Vector tmp1 = sol_in.root() + k2 * (delta /2);
        sol->setRoot(tmp1);
        Vector k3 = homotopy.impliciteDiff(t + delta/2, *sol);
        Vector tmp2 = sol_in.root() + k3 * (delta/2);
        sol->setRoot(tmp2);
        Vector k4 = homotopy.impliciteDiff(t + delta, *sol);
        delete sol;
        
        return k1 * (delta/6) + k2 * (delta/3) + k3 * (delta/3) + k4 * (delta/6);
}

Vector BetheHomotopyContinuation::moveOneStep(Homotopy &homotopy, Solution* sol, elem_t t, elem_t nextT) {
    elem_t tEnd = nextT;
    int depth = 0;
    elem_t delta = nextT - t;
    elem_t error;
    Vector r0 = sol->root();
    Vector last = sol->root();
    while(abs(t - tEnd) > EPS) {
//        Vector predict = predictChange(delta) * delta;
        Vector predict = _fitter.fit(t + delta);
        if (predict.size() > 0) {
//            sol->root() += predict;
            sol->setRoot(predict);
        } else {
            sol->root() += nextChange(homotopy, *sol, t, delta);
        }
        
        // if depth less than maxDepth, we check error.
        // if error is too big, we reduce delta to one-half of itself
        if (depth < maxDepth()) {
            error = homotopy.error(t + delta, *sol);
            if (error > 0.1L) {
                sol->setRoot(last);
                depth++;
                delta *= 0.5L;
                continue;
            }
        }
        
        // when delta is small enough, we can Newton-Raphson
        int j = 0;
        bool valid = true;
        for (j = 0; j < maxNumberOfNewtonRaphonEachStep(); j++) {
            Vector correction = newtonRaphson(*sol, homotopy, t + delta);
            if (correction.size() != sol->root().size()) {
                break;
            }
            sol->root() += correction;
            if (predict.size() > 0) {
                // if the root is too far away from predict,
                // it is very likely a bad root
//                if ((sol->root() - last - predict).norm()  > 0.5L) {
                if ((sol->root() - predict).norm()  > 0.01L) {
                    valid = false;
                    break;
                }
            }
            if (correction.norm() < accuracyLow()) {
                break;
            }
        }
        
        if (j == maxNumberOfNewtonRaphonEachStep()) {
            LOG("t=" << t + delta << ", run newtonRaphson more than " << maxNumberOfNewtonRaphonEachStep() << " times", Warning);
        }

        if (!valid || homotopy.error(t + delta, *sol) > 0.1) {
            if (depth >= maxDepth()) {
                throw HomotopyContinuationException(t + delta, delta, homotopy.error(t + delta, *sol));
            } else {
                sol->setRoot(last);
                depth++;
                delta *= 0.5L;
                continue;
            }
        }
        
        t += delta;
//        Vector slope = (sol->root() - last)/delta;
//        addLastChange(slope, delta);
        chop(sol->root(), 1e-18L);
        _fitter.addLastRoot(sol->root(), t);
//        elem_t gapNorm = slope.norm() * delta;
        elem_t gapNorm = (sol->root() - last).norm();
//        if (gapNorm > _maxGap) {
//            _maxGap = gapNorm;
//            _maxGapLocation = t;
//            LOG("t=" << t << ", maxGap=" << _maxGap, Info);
//        }
        
        // after update
        last = sol->root();
        _errors.addSample(homotopy.error(t, *sol), t);
        _gaps.addSample(gapNorm, t);
    }
    
    return sol->root() - r0;
}

void BetheHomotopyContinuation::solve(Homotopy &homotopy, Solution& sol) {
    Vector nr = sol.normalRoot();
    homotopy.setStartRoot(nr);
    elem_t gap = (homotopy.tEnd() - homotopy.tBegin()) / homotopy.getSteps();
    elem_t gapMin = gap / 32.0;
    _fitter.clear();
    _errors.clear();
    _gaps.clear();
    
    elem_t t = homotopy.tBegin();
    elem_t maxDiff = 0.0;
    elem_t maxDiff2 = 0.0;
    elem_t lastT = homotopy.tBegin();
//    _maxGap = 0.0;
    int largeErrorTimes = 0;
//    _maxGapLocation = homotopy.tBegin();
    for (int i = 0; i < homotopy.getSteps(); i++) {
        elem_t t = homotopy.tBegin() + gap * i;
        trace(lastT, &sol, homotopy.error(lastT, sol));
        Solution * sol2 = sol.copy();
        Vector change = moveOneStep(homotopy, sol2, lastT, t + gap);
        delete sol2;        
        sol.root() += change;
        sol.update(t + gap, &_fitter);
        lastT = t + gap;
    }
    
    for (int j = 0; j < 10; j++) {
        Vector correction = newtonRaphson(sol, homotopy, homotopy.tEnd());
        if (correction.size() != sol.root().size()) {
            break;
        }
        
        sol.root() += correction;
        if (correction.norm() < accuracyHigh()) {
            break;
        }
    }
    
    trace(homotopy.tEnd(), &sol, homotopy.error(homotopy.tEnd(), sol));
}

HomotopyContinuationException::HomotopyContinuationException(elem_t t, elem_t delta, elem_t error) : _t(t), _delta(delta), _error(error) {
    message = "Failed to continue with given accuracy. t=" + ToString(_t)+ ", delta=" + ToString(_delta) + ", error=" + ToString(_error);
}

const char* HomotopyContinuationException::what() const noexcept {
    return message.c_str();
}

