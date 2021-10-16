#pragma once
#include "common.h"

#define CACHE_SIZE 4
#define NEXT(i) ((i+1) % CACHE_SIZE)

class RootFitter
{
private:
    Vector _lastRoots[CACHE_SIZE];
    elem_t _beta[CACHE_SIZE];
    int _begin = 0;
    int _end = 0;
public:    
    void clear();
    
    void addLastRoot(const Vector& lastRoot, elem_t beta);
    
    Vector fit(elem_t beta) const;
    
    void invert(int i);
};
