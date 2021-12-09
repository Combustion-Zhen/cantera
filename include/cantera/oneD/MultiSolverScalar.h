//! @file MultiSolverScalar.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_MULTISOLVERSCALAR_H
#define CT_MULTISOLVERSCALAR_H

#include "cantera/numerics/BandMatrix.h"
#include "OneDim.h"

namespace Cantera
{

/**
 * Newton iterator for multi-domain, one-dimensional problems for scalars.
 * Used by class OneDim.
 * @ingroup onedim
 */
class MultiSolverScalar : public BandMatrix
{
public:

    //! Constructor
    MultiSolverScalar(OneDim& r);

    //!
    virtual ~MultiSolverScalar(){};

    MultiSolverScalar(const MultiSolverScalar&) = delete;

    MultiSolverScalar& operator=(const MultiSolverScalar&) = delete;

    /**
     * Find the solution to F(X) = 0 by damped Newton iteration. On entry, x0
     * contains an initial estimate of the solution. On successful return, x1
     * contains the converged solution.
     */
    int newtonSolve(double* x0, double* x1, int loglevel);

    //! Evaluate the Jacobian at m_x
    void evalJac();

    /**
     * On entry, step0 must contain an undamped Newton step for the solution x0.
     * This method attempts to find a damping coefficient such that the next
     * undamped step would have a norm smaller than that of step0. If
     * successful, the new solution after taking the damped step is returned in
     * x1, and the undamped step at x1 is returned in step1.
     */
    int dampStep(double* x1, int loglevel, bool writetitle);

    //! Compute the undamped Newton step.  The residual function is evaluated
    //! at `x`, but the Jacobian is not recomputed.
    void step(double* x, double* step, int loglevel);

    //! Compute the weighted 2-norm of `step`.
    double norm2(const vector_fp& x, const vector_fp& step) const;

    /**
     * Return the factor by which the undamped Newton step 'step0'
     * must be multiplied in order to keep all solution components in
     * all domains between their specified lower and upper bounds.
     */
    double boundStep(const vector_fp& x, const vector_fp& step, int loglevel);

    void convertFullToScalar(const vector_fp& full, vector_fp& scalar);

    void convertScalarToFull(const vector_fp& scalar, vector_fp& full);

    /// Set options.
    void setOptions(int maxJacAge = 5) {
        m_maxJacAge = maxJacAge;
    }

    inline int size() { 
        return m_n;
    }

    /////////// Jacobian

    //! Elapsed CPU time spent computing the Jacobian.
    inline double elapsedTimeJac() const {
        return m_elapsedJac;
    }

    //! Number of Jacobian evaluations.
    inline int nJacEval() const {
        return m_nJacEval;
    }


    //! Increment the Jacobian age.
    inline void incrementJacAge() {
        m_jacAge++;
    }

    //! Number of times 'incrementAge' has been called since the last evaluation
    inline int getJacAge() const {
        return m_jacAge;
    }

    //! Set the Jacobian age.
    inline void setJacAge(int age) {
        m_jacAge = age;
    }

protected:

    //! Work arrays of size #m_n used in solve().
    vector_fp m_x;

    //! Residual evaluator for this Jacobian
    /*!
     * This is a pointer to the residual evaluator. This object isn't owned by
     * this Jacobian object.
     */
    OneDim* m_resid;

    //! number of Jacobian calculations
    int m_jacAge;
    int m_nJacEval;
    int m_maxJacAge;

    //! perturbation parameters
    double m_atol, m_rtol;

    //! calculation time
    double m_elapsedJac;
    double m_elapsedNewton;
};
}

#endif
