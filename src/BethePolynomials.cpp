#include "BethePolynomials.h"
using namespace std::literals::complex_literals;

BethePolynomialBase::BethePolynomialBase(int L, int M) : _L(L), _M(M) {
}

var_t prefactorCoefLeftA[2][2] = {
    // uj + i/2
    {0.5il, 1.0L}, // regular
    // 1 + 0.5i * uj
    {1.0L, 0.5il}, // inverted
};

var_t prefactorCoefRightA[2][2] = {
    // uj - i/2
    {-0.5il, 1.0L}, // regular
    // 1 - 0.5i * uj
    {1.0L, -0.5il}, // inverted
};

var_t coefLeftA[2][2][4] = {
    {
        //{1,  uj,  uk, uj*uk}
        // uj - uk - i
        {-1.0il, 1.0L, -1.0L, 0.0L}, // [regular][regular]
        // uj * uk - 1 - i * uk
        {-1.0L, 0.0L, -1.0il, 1.0L}, // [regular][inverted]
    },
    
    {
        // 1 - uj * uk - i * uj
        {1.0L, -1.0il, 0.0L, -1.0}, // [inverted][regular]
        // uk - uj - i * uj * uk
        {0.0L, -1.0L, 1.0L, -1.0il}, // [inverted][inverted]
    },
};

var_t coefRightA[2][2][4] = {
    {
        //{1,     uj,   uk,  uj*uk}
        // uj - uk + i
        {1.0il, 1.0L, -1.0L, 0.0L}, // [regular][regular]
        // uj * uk - 1 + i * uk
        {-1.0L, 0.0L, 1.0il, 1.0L}, // [regular][inverted]
    },
    
    {
        // 1 - uj * uk + i * uj
        {1.0L, 1.0il, 0.0L, -1.0}, // [inverted][regular]
        // uk - uj + i * uj * uk
        {0.0L, -1.0L, 1.0L, 1.0il}, // [inverted][inverted]
    },
};

