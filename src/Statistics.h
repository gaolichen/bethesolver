#pragma once
#include "common.h"
#include <complex>

template<class T> class Statistics
{
private:
    T _sum;
    T _max;
    T _min;
    T _maxLocation;
    T _minLocation;
    T _squareSum;
    int _num;
public:
    Statistics() {
        clear();
    }
    
    void clear() {
        _sum = (T)0;
        _num = 0;
        _max = (T)0;
        _min = (T)0;
        _squareSum = (T)0;
    }
    
    void addSample(T val, T location) {
        _sum += val;
        _num++;
        if (_num == 1) {
            _max = _min = val;
            _maxLocation = _minLocation = location;
        } else {
            if (val > _max) {
                _max = val;
                _maxLocation = location;
            }
            if (val < _min) { 
                _min = val;
                _minLocation = location;
            }
        }
        _squareSum += val * val;
    }
    
    int num() const {
        return _num;
    }
    
    T max() const {
        return _max;
    }
    
    T maxLocation() const {
        return _maxLocation;
    }
    
    T minLocation() const {
        return _minLocation;
    }
    
    T min() const {
        return _min;
    }
    
    T average() const {
        return _sum / (T)_num;
    }
    
    T deviation() const {
        return sqrt(_squareSum / _num - _sum * _sum / T(_num * _num));
    }
};

template <typename T>
std::ostream& operator<< (std::ostream& out, const Statistics<T>& st)
{
  out << "{range=(" <<  st.min() << ", " << st.max() << "), locations of min/max=(" << st.minLocation() << ", " << st.maxLocation() << "), iteration number=" << st.num() << ", average=" << st.average() << ", deviation=" << st.deviation() << '}';
  return out;
}

