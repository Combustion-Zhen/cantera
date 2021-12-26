//! @file MultiSolverScalar.cpp Damped Newton solver for 1D multi-domain scalars

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/oneD/MultiSolverScalar.h"
#include "cantera/oneD/OneDimConst.h"
#include "cantera/base/utilities.h"

#include <ctime>

using namespace std;

namespace Cantera
{

// unnamed-namespace for local helpers
namespace
{

class Indx
{

public:
    Indx(size_t nv, size_t np) : m_nv(nv), m_np(np) {}
    size_t m_nv, m_np;
    size_t operator()(size_t m, size_t j) {
        return j*m_nv + m;
    }
};

/**
 * This function computes the square of a weighted norm of a step vector for one
 * domain.
 *
 * @param x     Solution vector for this domain.
 * @param step  Newton step vector for this domain.
 * @param r     Object representing the domain. Used to get tolerances,
 *              number of components, and number of points.
 *
 * The return value is
 * \f[
 *    \sum_{n,j} \left(\frac{s_{n,j}}{w_n}\right)^2
 * \f]
 * where the error weight for solution component \f$n\f$ is given by
 * \f[
 *     w_n = \epsilon_{r,n} \frac{\sum_j |x_{n,j}|}{J} + \epsilon_{a,n}.
 * \f]
 * Here \f$\epsilon_{r,n} \f$ is the relative error tolerance for component n,
 * and multiplies the average magnitude of solution component n in the domain.
 * The second term, \f$\epsilon_{a,n}\f$, is the absolute error tolerance for
 * component n.
 */
double norm_square(const double* x, const double* step, Domain1D& r)
{
    double sum = 0.0;
    double f2max = 0.0;
    size_t nv = r.nScalars();
    size_t np = r.nPoints();

    for (size_t n = 0; n != nv; n++) 
    {
        double esum = 0.0;
        for (size_t j = 0; j != np; j++) 
        {
            esum += fabs(x[nv*j + n]);
        }
        double ewt = r.rtol(n)*esum/np + r.atol(n);
        for (size_t j = 0; j != np; j++) 
        {
            double f = step[nv*j + n]/ewt;
            sum += f*f;
            f2max = std::max(f*f, f2max);
        }
    }
    return sum;
}

/**
 * Return a damping coefficient that keeps the solution after taking one
 * Newton step between specified lower and upper bounds. This function only
 * considers one domain.
 */
double bound_step(const double* x, const double* step, Domain1D& r, int loglevel)
{
    size_t nv = r.nScalars();
    size_t np = r.nPoints();

    Indx index(nv, np);
    double fbound = 1.0;
    bool writeTitle = false;

    for (size_t m = 0; m != nv; m++) {
        double above = r.upperBoundScalar(m);
        double below = r.lowerBoundScalar(m);

        for (size_t j = 0; j != np; j++) {
            double val = x[index(m,j)];
            if (loglevel > 0 && (val > above + 1.0e-12 || val < below - 1.0e-12)) {
                writelog("\nERROR: solution out of bounds.\n");
                writelog("domain {:d}: {:>20s}({:d}) = {:10.3e} ({:10.3e}, {:10.3e})\n",
                         r.domainIndex(), r.componentName(m), j, val, below, above);
            }

            double newval = val + step[index(m,j)];

            if (newval > above) {
                fbound = std::max(0.0, std::min(fbound,
                                                (above - val)/(newval - val)));
            } else if (newval < below) {
                fbound = std::min(fbound, (val - below)/(val - newval));
            }

            if (loglevel > 1 && (newval > above || newval < below)) {
                if (!writeTitle) {
                    writelog("\nNewton step takes solution out of bounds.\n\n");
                    writelog("  {:>12s}  {:>12s}  {:>4s}  {:>10s}  {:>10s}  {:>10s}  {:>10s}\n",
                             "domain","component","pt","value","step","min","max");
                    writeTitle = true;
                }
                writelog("          {:4d}  {:>12s}  {:4d}  {:10.3e}  {:10.3e}  {:10.3e}  {:10.3e}\n",
                         r.domainIndex(), r.componentName(m), j,
                         val, step[index(m,j)], below, above);
            }
        }
    }
    return fbound;
}

} // end unnamed-namespace

// constants
const double DampFactor = 2.0;
const size_t NDAMP = 8;

// ---------------- MultiSolverScalar methods ----------------

MultiSolverScalar::MultiSolverScalar(OneDim& r) : 
    //BandMatrix(r.sizeScalar(),r.bandwidthScalar(),r.bandwidthScalar()), 
    m_resid(&r),
    m_jacAge(0), m_nJacEval(0), m_maxJacAge(5),
    m_atol(sqrt(std::numeric_limits<double>::epsilon())), m_rtol(1.0e-5),
    m_elapsedJac(0.0), m_elapsedNewton(0.0)
{}

void MultiSolverScalar::resize()
{
    BandMatrix::resize(m_resid->sizeScalar(),
                       m_resid->bandwidthScalar(),
                       m_resid->bandwidthScalar());
    m_scalar0.resize(m_resid->sizeScalar());
    m_scalar1.resize(m_resid->sizeScalar());
    m_step0.resize(m_resid->sizeScalar());
    m_step1.resize(m_resid->sizeScalar());
}

int MultiSolverScalar::dampedNewtonSolve(double* x0, double* x1, int loglevel)
{
    clock_t t0 = clock();

    bool forceNewJac = true;
    bool writeTitle = true;
    int nJacReeval = 0;
    int m = 0;

    while (true) 
    {
        // Check whether the Jacobian should be re-evaluated.
        if (getJacAge() > m_maxJacAge) 
        {
            if (loglevel > 0) 
                writelog("\nMaximum Jacobian age reached ({})\n", m_maxJacAge);

            forceNewJac = true;
        }

        if (forceNewJac) 
        {
            evalJac(x0);
            forceNewJac = false;
        }

        // increment the Jacobian age
        incrementJacAge();

        // damp the Newton step
        m = dampStep(x0, x1, loglevel-1, writeTitle);

        writeTitle = false;

        //if (loglevel > 0)
        //    writelog("\n\n reverse condition number {:10.4g}", rcond(oneNorm()));

        // Successful step, but not converged yet. Take the damped step, and try
        // again.
        if (m == 0) 
        {
            copy(x1, x1 + m_resid->size(), x0);
        } 
        else if (m == 1) 
        {
            break;
        } 
        else if (m < 0) 
        {
            // If dampStep fails, first try a new Jacobian if an old one was
            // being used. If it was a new Jacobian, then return -1 to signify
            // failure.
            if ( getJacAge() > 1 && nJacReeval < 4 ) 
            {
                forceNewJac = true;
                nJacReeval++;
                debuglog("\n\nRe-evaluating Jacobian, since no damping "
                         "coefficient\ncould be found with this Jacobian.\n",
                         loglevel);
            }
            else
            {
                break;
            }
        }
    }

    if (m < 0)
        copy(x0, x0 + m_resid->size(), x1);

    if (m > 0 && nJacEval() == 1)
        m = 100;

    m_elapsedNewton += (clock() - t0)/(1.0*CLOCKS_PER_SEC);

    return m;
}

void MultiSolverScalar::evalJac(double* x)
{
    clock_t t0 = clock();

    bfill(0.0);
    incrementJacEval();

    // evaluate the unperturbed residual
    m_resid->evalScalar(npos, x, m_step0.data(), 0);

    // perturb the full solution vector to obtain the Jacobian of scalars
    for (size_t j = 0; j != m_resid->points(); j++) 
    {
        // the number of scalars
        size_t nv = m_resid->nVarScalar(j);
        // location of the first variable in scalar and full solution vector for point j
        size_t jFull = m_resid->loc(j);
        size_t jScalar = m_resid->locScalar(j);
        // iterate over scalars
        for (size_t n = 0; n < nv; n++) 
        {
            // location of the scalar
            size_t iScalar = jScalar + n;
            // location of the scalar in the full solution vector
            size_t offset = (n==0) ? c_offset_T : n-cOffsetScalarY+c_offset_Y ;
            size_t iFull = jFull + offset;

            // perturb x(n); preserve sign(x(n))
            double tmp = x[iFull];
            double dx = (tmp>=0.0) ? tmp*m_rtol + m_atol : tmp*m_rtol - m_atol; 
            double rdx = 1.0/dx;
            x[iFull] = tmp + dx;

            // calculate perturbed residual
            m_resid->evalScalar(j, x, m_step1.data(), 0);

            // compute nth column of Jacobian
            for (size_t i = j - 1; i != j+2; i++) 
            {
                if (i != npos && i < m_resid->points() ) 
                {
                    size_t mv = m_resid->nVarScalar(i);
                    size_t iloc = m_resid->locScalar(i);
                    for (size_t m = 0; m < mv; m++) 
                    {
                        value(iloc+m,iScalar) = (m_step1[iloc+m]-m_step0[iloc+m])*rdx;
                    }

                }
            }

            // recover the x(n) value
            x[iFull] = tmp;
        }
    }

    m_elapsedJac += double(clock() - t0)/CLOCKS_PER_SEC;
    setJacAge(0);
}

int MultiSolverScalar::dampStep(double* x0, double* x1, int loglevel, bool writeTitle)
{
    double s0, s1;

    // write header
    if (loglevel > 0 && writeTitle) 
    {
        writelog("\n\nDamped Newton iteration:\n");
        writeline('-', 65, false);

        writelog("\n{}  {:>9s}   {:>9s}   {:>9s}   {:>9s}   {:>5s} {:>5s}\n",
                "m", "F_damp", "F_bound", "log10(s0)", "log10(s1)", "N_jac", "Age");
        writeline('-', 65);
    }

    copyFullToScalar(x0, m_scalar0.data());

    // compute the undamped Newton step
    step(x0, m_step0.data(), loglevel-1);

    // compute the weighted norm of the undamped step size step0
    s0 = norm2(m_scalar0, m_step0);

    // compute the multiplier to keep all components in bounds
    double fbound = boundStep(m_scalar0, m_step0, loglevel-1);

    // if fbound is very small, then x0 is already close to the boundary and
    // step0 points out of the allowed domain. In this case, the Newton
    // algorithm fails, so return an error condition.
    if (fbound < 1.e-10)
    {
        debuglog("\nAt limits.\n", loglevel);
        return -3;
    }

    // ---------- Attempt damped step ----------

    // damping coefficient starts at 1.0
    double damp = 1.0;
    size_t m;

    copy(x0, x0 + m_resid->size(), x1);
    for (m = 0; m != NDAMP; m++) {
        double ff = fbound*damp;

        // step the solution by the damped step size
        for (size_t j = 0; j < m_resid->sizeScalar(); j++) {
            m_scalar1[j] = m_scalar0[j] + ff*m_step0[j];
        }
        copyScalarToFull(m_scalar1.data(), x1);

        // compute the next undamped step that would result if x1 is accepted
        step(x1, m_step1.data(), loglevel-1);

        // compute the weighted norm of step1
        s1 = norm2(m_scalar1, m_step1);

        // write log information
        if (loglevel > 0)
        {
            writelog("\n{:d}  {:9.5f}   {:9.5f}   {:9.5f}   {:9.5f}    {:4d}  {:d}/{:d}",
                     m, damp, fbound, log10(s0+SmallNumber), log10(s1+SmallNumber),
                     nJacEval(), getJacAge(), m_maxJacAge);
        }

        // if the norm of s1 is less than the norm of s0, then accept this
        // damping coefficient. Also accept it if this step would result in a
        // converged solution. Otherwise, decrease the damping coefficient and
        // try again.
        if (s1 < 1.0 || s1 < s0)
            break;

        damp /= DampFactor;
    }

    // If a damping coefficient was found, return 1 if the solution after
    // stepping by the damped step would represent a converged solution, and
    // return 0 otherwise. If no damping coefficient could be found, return -2.
    if (m < NDAMP) 
        return (s1 > 1.0) ? 0 : 1 ;
    else 
        return -2;
}

void MultiSolverScalar::step(double* x, double* step, int loglevel)
{
    m_resid->evalScalar(npos, x, step);
    for (size_t n = 0; n != m_resid->sizeScalar(); n++) 
    {
        step[n] = -step[n];
    }

    try 
    {
        this->solve(step, step);
    } 
    catch (CanteraError&) 
    {
        if (this->info() > 0) 
        {
            // Positive value for "info" indicates the row where factorization failed
            size_t row = static_cast<size_t>(this->info() - 1);
            // Find the domain, grid point, and solution component corresponding
            // to this row
            for (size_t n = 0; n < m_resid->nDomains(); n++) 
            {
                Domain1D& dom = m_resid->domain(n);
                size_t nComp = dom.nScalars();
                if (row >= dom.locScalar() && row < dom.locScalar() + dom.sizeScalar()) 
                {
                    size_t offset = row - dom.locScalar();
                    size_t pt = offset / nComp;
                    size_t comp = offset - pt * nComp;
                    throw CanteraError("MultiSolverScalar::step",
                        "Jacobian is singular for domain {}, component {} at point {}\n"
                        "(Matrix row {})",
                        dom.id(), dom.componentName(comp), pt, row);
                }
            }
        }
        throw;
    }
}

double MultiSolverScalar::norm2(const vector_fp& x, const vector_fp& step) const
{
    double sum = 0.0;
    for (size_t n = 0; n != m_resid->nDomains(); n++) 
    {
        sum += norm_square(&x[m_resid->startScalar(n)],
                           &step[m_resid->startScalar(n)],
                           m_resid->domain(n));
    }
    sum /= m_resid->sizeScalar();
    return sqrt(sum);
}

double MultiSolverScalar::boundStep(const vector_fp& x, const vector_fp& step, int loglevel)
{
    double fbound = 1.0;
    for (size_t i = 0; i != m_resid->nDomains(); i++) 
    {
        fbound = std::min(fbound,
                          bound_step(&x[m_resid->startScalar(i)],
                                     &step[m_resid->startScalar(i)],
                                     m_resid->domain(i),
                                     loglevel));
    }
    return fbound;
}

void MultiSolverScalar::copyFullToScalar(const double* full, double* scalar)
{
    for (size_t j = 0 ; j != m_resid->points() ; j++ )
    {
        // location of the first variable in scalar and full solution vector for point j
        size_t jFull = m_resid->loc(j);
        size_t jScalar = m_resid->locScalar(j);
        // the number of scalars
        size_t nScalar = m_resid->nVarScalar(j);
        // iterate over scalars
        for (size_t n = 0; n != nScalar; n++) 
        {
            // location of the scalar
            size_t iScalar = jScalar + n;
            // location of the scalar in the full solution vector
            size_t m = (n==0) ? c_offset_T : n-cOffsetScalarY+c_offset_Y ;
            size_t iFull = jFull + m;

            scalar[iScalar] = full[iFull];
        }
    }
}

void MultiSolverScalar::copyScalarToFull(const double* scalar, double* full)
{
    for (size_t j = 0 ; j != m_resid->points() ; j++ )
    {
        // location of the first variable in scalar and full solution vector for point j
        size_t jFull = m_resid->loc(j);
        size_t jScalar = m_resid->locScalar(j);
        // the number of scalars
        size_t nScalar = m_resid->nVarScalar(j);
        // iterate over scalars
        for (size_t n = 0; n != nScalar; n++) 
        {
            // location of the scalar
            size_t iScalar = jScalar + n;
            // location of the scalar in the full solution vector
            size_t m = (n==0) ? c_offset_T : n-cOffsetScalarY+c_offset_Y ;
            size_t iFull = jFull + m;

            full[iFull] = scalar[iScalar];
        }
    }
}

} // end namespace Cantera
