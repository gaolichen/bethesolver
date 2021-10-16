#pragma once
#include "common.h"

class InfiniteRoots
{
private:
    std::vector<std::vector<var_t> > _roots;
    
    int _L;
    int _M;
    
    void computeRoots(int m);
    
    static std::vector<var_t> equationSolver(std::vector<var_t> &coefs, std::vector<var_t> &out);
public:
    InfiniteRoots(int L, int M);
    
    std::vector<var_t>& getRoot(int m);
};
