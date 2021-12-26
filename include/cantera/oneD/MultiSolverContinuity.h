//! @file MultiSolverContinuity.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_MULTISOLVERCONTINUITY_H
#define CT_MULTISOLVERCONTINUITY_H

#include "OneDim.h"

namespace Cantera
{

/**
 * Newton iterator for multi-domain, one-dimensional problems for scalars.
 * Used by class OneDim.
 * @ingroup onedim
 */
class MultiSolverContinuity
{
public:

    //! Constructor
    MultiSolverContinuity(OneDim& r);

    //!
    virtual ~MultiSolverContinuity(){};

    MultiSolverContinuity(const MultiSolverContinuity&) = delete;

    MultiSolverContinuity& operator=(const MultiSolverContinuity&) = delete;

    void copyFullToVelocity(const double* full, double* velocity);

    void copyVelocityToFull(const double* velocity, double* full);

    //! solve the continuity equation to get the velocity
    //! this is a linear problem
    int newtonSolve(double* x0, double* x1, int loglevel);

protected:

    //! Residual evaluator for this Jacobian
    /*!
     * This is a pointer to the residual evaluator. This object isn't owned by
     * this Jacobian object.
     */
    OneDim* m_resid;

    size_t m_n;

    //! Work arrays of size #m_n used in solve().
    vector_fp m_velocity, m_step, m_dl, m_d, m_du;

};
}

#endif
