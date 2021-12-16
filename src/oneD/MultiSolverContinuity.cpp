//! @file MultiSolverContinuity.cpp 

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/oneD/MultiSolverContinuity.h"
#include "cantera/oneD/OneDimConst.h"
#include "cantera/base/utilities.h"

#include <ctime>

using namespace std;

namespace Cantera
{

// ---------------- MultiSolverContinuity methods ----------------

// public

MultiSolverContinuity::MultiSolverContinuity(OneDim& r) :
    TridiagonalMatrix(r.size())
{
    m_resid = &r;
    m_x.resize(m_resid->size());
}

void MultiSolverContinuity::convertFullToVelocity(const vector_fp& full, vector_fp& velocity)
{
    for (size_t i = 0 ; i != m_resid->points() ; i++ )
    {
        // location of the first variable in scalar and full solution vector for point j
        size_t j = m_resid->loc(i);
        size_t k = m_resid->locVelocity(i);
        // take the velocity, which is stored as the first component
        velocity[k] = full[j];
    }
}

void MultiSolverContinuity::convertVelocityToFull(const vector_fp& velocity, vector_fp& full)
{
    for (size_t i = 0 ; i != m_resid->points() ; i++ )
    {
        // location of the first variable in scalar and full solution vector for point j
        size_t j = m_resid->loc(i);
        size_t k = m_resid->locVelocity(i);
        //
        full[j] = velocity[k];
    }
}

int MultiSolverContinuity::newtonSolve(double* x0, double* x1, int loglevel)
{
    copy(x0, x0+m_resid->size(), m_x.begin());

    // get the matrix

    // get the rhs

    // solve the tridiagonal problem

    // int m = this->solve();

    return 0;
}

} // end namespace Cantera
