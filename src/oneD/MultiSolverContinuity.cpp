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
    m_resid(&r)
{
    m_n = m_resid->sizeVelocity();
    m_velocity.resize(m_n);
    m_step.resize(m_n);
    m_dl.resize(m_n-1);
    m_d.resize(m_n);
    m_du.resize(m_n-1);
}

void MultiSolverContinuity::copyFullToVelocity(const double* full, double* velocity)
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

void MultiSolverContinuity::copyVelocityToFull(const double* velocity, double* full)
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
    // get the residual and Jacobian from the continuity equation
    m_resid->evalContinuity(x0, m_step.data(), m_dl.data(), m_d.data(), m_du.data());

    copyFullToVelocity(x0, m_velocity.data());
    // solve the tridiagonal problem
    for (size_t i = 0; i != m_n; i++)
    {
        m_velocity[i] += m_step[i]/m_d[i];
    }

    copy(x0, x0+m_resid->size(), x1);
    copyVelocityToFull(m_velocity.data(), x1);
    
    if (loglevel > 0)
    {
        writelog("\n\n velocity    residual    diagonal    s-diagonals \n");
        writelog("===============left  boundary===============\n");
        writelog(" {:10.4g} {:10.4g} {:10.4g} {:10.4g}\n", 
                 m_velocity[0], m_step[0], m_d[0], m_du[0]);
        writelog(" {:10.4g} {:10.4g} {:10.4g} {:10.4g} {:10.4g}\n", 
                 m_velocity[1], m_step[1], m_d[1], m_dl[0], m_du[1]);
        writelog("===============right boundary===============\n");
        writelog(" {:10.4g} {:10.4g} {:10.4g} {:10.4g} {:10.4g}\n", 
                 m_velocity[m_n-2], m_step[m_n-2], m_d[m_n-2], m_dl[m_n-3], m_du[m_n-2]);
        writelog(" {:10.4g} {:10.4g} {:10.4g} {:10.4g}\n", 
                 m_velocity[m_n-1], m_step[m_n-1], m_d[m_n-1], m_dl[m_n-2]);
    }

    return 0;
}

} // end namespace Cantera
