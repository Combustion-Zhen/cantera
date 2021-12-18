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
    TridiagonalMatrix(r.sizeVelocity())
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

    vector_fp velocity(m_resid->sizeVelocity(), 0.0);
    vector_fp step(m_resid->sizeVelocity(), 0.0);

    copy(x0, x0+m_resid->size(), m_x.begin());

    convertFullToVelocity(m_x, velocity);

    // get the residual and Jacobian from the continuity equation
    m_resid->evalContinuity(m_x, step, m_dl, m_d, m_du);

    if (loglevel>2)
    {
        writelog("\n {:10.4g} {:10.4g} {:10.4g}", 
                step[0], m_d[0], m_du[0]);
        writelog("\n {:10.4g} {:10.4g} {:10.4g}\n", 
                step[m_resid->sizeVelocity()-1], 
                m_d[m_resid->sizeVelocity()-1], 
                m_dl[m_resid->sizeVelocity()-2]);
        writelog("\n {:10.4g} {:10.4g} {:10.4g} {:10.4g}", 
                step[1], m_dl[0], m_d[1], m_du[1]);
    }

    // solve the tridiagonal problem
    //int m = this->solve(step);
    for (size_t i = 0; i != m_resid->sizeVelocity(); i++ )
    {
        //velocity[i] += step[i];
        velocity[i] += step[i]/m_d[i];
    }
    convertVelocityToFull(velocity, m_x);

    if (loglevel>2)
    {
        writelog("\n {:10.4g} {:10.4g}", 
                velocity[0], step[0]);
        writelog("\n {:10.4g} {:10.4g}", 
                velocity[m_resid->sizeVelocity()-1], 
                step[m_resid->sizeVelocity()-1]);
    }

    copy(m_x.begin(), m_x.end(), x1);

    return 0;
}

} // end namespace Cantera
