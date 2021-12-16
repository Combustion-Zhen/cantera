//! @file TridiagonalMatrix.cpp Tridiagonal matrices.

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/numerics/TridiagonalMatrix.h"
#include "cantera/base/utilities.h"
#include "cantera/base/stringUtils.h"

#if CT_USE_LAPACK
    #include "cantera/numerics/ctlapack.h"
#else
    #if CT_SUNDIALS_VERSION >= 30
        #include "sunlinsol/sunlinsol_band.h"
    #else
        #include "cvodes/cvodes_dense.h"
        #include "cvodes/cvodes_band.h"
    #endif
#endif

#include <cstring>
#include <fstream>

using namespace std;

namespace Cantera
{

TridiagonalMatrix::TridiagonalMatrix() :
    m_n(0)
{
}

TridiagonalMatrix::TridiagonalMatrix(size_t n, double v) :
    m_n(n)
{
    m_d.resize(m_n);
    m_dl.resize(m_n-1);
    m_du.resize(m_n-1);
    fill(m_d.begin(), m_d.end(), v);
    fill(m_dl.begin(), m_dl.end(), v);
    fill(m_du.begin(), m_du.end(), v);
}

void TridiagonalMatrix::resize(size_t n, double v)
{
    m_n = n;

    m_d.resize(m_n);
    m_dl.resize(m_n-1);
    m_du.resize(m_n-1);
    fill(m_d.begin(), m_d.end(), v);
    fill(m_dl.begin(), m_dl.end(), v);
    fill(m_du.begin(), m_du.end(), v);
}

void TridiagonalMatrix::bfill(double v)
{
    fill(m_d.begin(), m_d.end(), v);
    fill(m_dl.begin(), m_dl.end(), v);
    fill(m_du.begin(), m_du.end(), v);
}

int TridiagonalMatrix::solve(vector_fp& b, size_t nrhs, size_t ldb)
{
    int info = 0;

    if (ldb == 0) 
    {
        ldb = nColumns();
    }

#if CT_USE_LAPACK
    ct_dgtsv(size(), nrhs, 
             m_dl.data(), m_d.data(), m_du.data(), b.data(), 
             ldb, info);
#else
    throw Cantera::CanteraError("TridiagonalMatrix::solve", 
                                "Tridiagonal problem solver not set");
#endif

    if (info != 0) {
        throw Cantera::CanteraError("TridiagonalMatrix::solve",
            "Linear solve failed with DGTSV error code {}.", info);
    }
    return info;
}

}
