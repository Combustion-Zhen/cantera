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

    void convertFullToVelocity(const vector_fp& full, vector_fp& velocity);

    void convertVelocityToFull(const vector_fp& velocity, vector_fp& full);

    //! solve the continuity equation to get the velocity
    //! this is a linear problem
    int newtonSolve(double* x0, double* x1, int loglevel);

protected:

    //! Work arrays of size #m_n used in solve().
    vector_fp m_x;

    //! Residual evaluator for this Jacobian
    /*!
     * This is a pointer to the residual evaluator. This object isn't owned by
     * this Jacobian object.
     */
    OneDim* m_resid;

};
}

#endif
