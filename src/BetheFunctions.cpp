#include "BetheFunctions.h"
using namespace std::literals::complex_literals;

BethePolynomial::BethePolynomial(int L, int M, bool isNumerator) {
    _L = L;
    _M = M;
    if (isNumerator) {
        _iUnit = var_t(0, 1.0);
    } else {
        _iUnit = var_t(0, -1.0);
    }
}

//=================================================

var_t prefactorCoefLeft[3][3] = {
    // uj + i/2
    {0.5il, 1.0L, 0.0L}, // regular
    // 1 + 0.5i * omega
    {1.0L, 0.0L, 0.5il}, // largest
    // uj + 0.5i * omega
    {0.0L, 1.0L, 0.5il}  // large
};

var_t prefactorCoefRight[3][3] = {
    // uj - i/2
    {-0.5il, 1.0L, 0.0L}, // regular 
    // 1 - 0.5i omega
    {1.0L, 0.0L, -0.5il}, // largest
    // uj - 0.5i omega
    {0.0L, 1.0L, -0.5il}  // large
};

var_t coefLeft[3][3][6] = {
    {
        //{i,     uj,   uk,  omega, omega*uj, omega*uk}
        // uj - uk - i
        {-1.0il, 1.0L, -1.0L, 0.0L, 0.0L, 0.0L}, // [regular][regular]
        // uj * omega - 1 - i * omega
        {-1.0L, 0.0L, 0.0L, -1.0il, 1.0L, 0.0L}, // [regular][largest]
        // uj * omega - uk - i * omega
        {0.0L, 0.0L, -1.0L, -1.0il, 1.0L, 0.0L}  // [regular][large]
    },
    
    {
        // 1 - uk * omega - i * omega
        {1.0L, 0.0L, 0.0L, -1.0il, 0.0L, -1.0}, // [largest][regular]
        // not exist
        {0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L}, // [largest][large]
        // 1 - uk - i * omega
        {1.0L, 0.0, -1.0, -1.0il, 0.0L, 0.0L}  // [largest][large]
    },
    
    {
        // uj - uk * omega - i * omega
        {0.0L, 1.0L, 0.0L, -1.0il, 0.0L, -1.0L}, // [large][regular]
        // uj - 1 - i * omega
        {-1.0L, 1.0L, 0.0L, -1.0il, 0.0L, 0.0L}, // [large][largest]
        // uj - uk - i * omega
        {0.0L, 1.0L, -1.0L, -1.0il, 0.0L, 0.0L}  // [large][large]
    },
};

var_t coefRight[3][3][6] = {
    {
        //{i,     uj,   uk,  omega, omega*uj, omega*uk}
        // uj - uk + i
        {1.0il, 1.0L, -1.0L, 0.0L, 0.0L, 0.0L}, // [regular][regular]
        // uj * omega - 1 + i * omega
        {-1.0L, 0.0L, 0.0L, 1.0il, 1.0L, 0.0L}, // [regular][largest]
        // uj * omega - uk + i * omega
        {0.0L, 0.0L, -1.0L, 1.0il, 1.0L, 0.0L}  // [regular][large]
    },
    
    {
        // 1 - uk * omega - i * omega
        {1.0L, 0.0L, 0.0L, 1.0il, 0.0L, -1.0}, // [largest][regular]
        // not exist
        {0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L}, // [largest][large]
        // 1 - uk + i * omega
        {1.0L, 0.0, -1.0, 1.0il, 0.0L, 0.0L}  // [largest][large]
    },
    
    {
        // uj - uk * omega - i * omega
        {0.0L, 1.0L, 0.0L, 1.0il, 0.0L, -1.0L}, // [large][regular]
        // uj - 1 + i * omega
        {-1.0L, 1.0L, 0.0L, 1.0il, 0.0L, 0.0L}, // [large][largest]
        // uj - uk - i * omega
        {0.0L, 1.0L, -1.0L, 1.0il, 0.0L, 0.0L}  // [large][large]
    },
};


BetheFactorBase::BetheFactorBase(int L, int M) {
    _L = L;
    _M = M;
}

var_t vectorDot(std::vector<var_t>& a, std::vector<var_t>& b) {
    var_t ret = (elem_t)0.0;
    for (int i = 0; i < a.size(); i++) {
        ret += a[i] * b[i];
    }
    return ret;
}

