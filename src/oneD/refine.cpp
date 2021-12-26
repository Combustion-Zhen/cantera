//! @file refine.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/oneD/refine.h"
#include "cantera/oneD/StFlow.h"
#include "cantera/base/global.h"

using namespace std;

namespace Cantera
{
Refiner::Refiner(Domain1D& domain) :
    m_ratio(10.0), m_slope(0.8), m_curve(0.8), m_prune(-0.001),
    m_min_range(0.01), m_domain(&domain), m_npmax(1000),
    m_gridmin(1e-10), m_np(0)
{
    m_nv = m_domain->nComponents();
    m_active.resize(m_nv, true);
    m_thresh = std::sqrt(std::numeric_limits<double>::epsilon());
    m_c.resize(m_nv);
}

void Refiner::setCriteria(doublereal ratio, doublereal slope,
                          doublereal curve, doublereal prune)
{
    if (ratio < 2.0) {
        throw CanteraError("Refiner::setCriteria",
            "'ratio' must be greater than 2.0 ({} was specified).", ratio);
    } else if (slope < 0.0 || slope > 1.0) {
        throw CanteraError("Refiner::setCriteria",
            "'slope' must be between 0.0 and 1.0 ({} was specified).", slope);
    } else if (curve < 0.0 || curve > 1.0) {
        throw CanteraError("Refiner::setCriteria",
            "'curve' must be between 0.0 and 1.0 ({} was specified).", curve);
    } else if (prune > curve || prune > slope) {
        throw CanteraError("Refiner::setCriteria",
            "'prune' must be less than 'curve' and 'slope' ({} was specified).",
            prune);
    }
    m_ratio = ratio;
    m_slope = slope;
    m_curve = curve;
    m_prune = prune;
}

void Refiner::setActives(const vector_int& comp)
{
    for (size_t j = 0; j != m_nv; j++)
    {
        setActive(j, false);
    }
    for (auto i : comp)
    {
        setActive(i, true);
    }
}

int Refiner::analyze(size_t n, const doublereal* z,
                     const doublereal* x)
{
    if (n >= m_npmax) {
        throw CanteraError("Refiner::analyze", "max number of grid points reached ({}).", m_npmax);
    }

    // check consistency
    if (n != m_domain->nPoints()) {
        throw CanteraError("Refiner::analyze", "inconsistent");
    }

    m_v.resize(n);
    m_s.resize(n-1);

    m_keep.resize(n);
    fill(m_keep.begin(), m_keep.end(), 0);
    m_keep[0] = 1;
    m_keep[n-1] = 1;

    m_loc.resize(n);
    fill(m_loc.begin(), m_loc.end(), 0);

    //m_c.clear();
    fill(m_c.begin(), m_c.end(), 0);

    if (m_domain->nPoints() <= 1) {
        return 0;
    }

    m_nv = m_domain->nComponents();
    for (size_t i = 0; i < m_nv; i++) {
        //writelog("Refiner::analyze active compnent {} {}\n",
        //         m_domain->componentName(i),
        //         m_active[i]);
        if (m_active[i]) {
            string name = m_domain->componentName(i);
            // get component i at all points
            for (size_t j = 0; j < n; j++) {
                m_v[j] = value(x, i, j);
            }

            // slope of component i
            for (size_t j = 0; j < n-1; j++) {
                m_s[j] = (value(x, i, j+1) - value(x, i, j))/(z[j+1] - z[j]);
            }

            // find the range of values and slopes
            doublereal vmin = *min_element(m_v.begin(), m_v.end());
            doublereal vmax = *max_element(m_v.begin(), m_v.end());
            doublereal smin = *min_element(m_s.begin(), m_s.end());
            doublereal smax = *max_element(m_s.begin(), m_s.end());

            // max absolute values of v and s
            doublereal aa = std::max(fabs(vmax), fabs(vmin));
            doublereal ss = std::max(fabs(smax), fabs(smin));

            // refine based on component i only if the range of v is
            // greater than a fraction 'min_range' of max |v|. This
            // eliminates components that consist of small fluctuations
            // on a constant background.
            if ((vmax - vmin) > m_min_range*aa) {
                // maximum allowable difference in value between adjacent
                // points.
                doublereal dmax = m_slope*(vmax - vmin) + m_thresh;
                for (size_t j = 0; j < n-1; j++) {
                    double dz = m_domain->dz(j);
                    double r = fabs(m_v[j+1] - m_v[j])/dmax;
                    if (r > 1.0 && dz >= 2 * m_gridmin) {
                        m_loc[j] = 1;
                        //m_c[name] = 1;
                        m_c[i] = 1;
                    }
                    if (r >= m_prune) {
                        m_keep[j] = 1;
                        m_keep[j+1] = 1;
                    } else if (m_keep[j] == 0) {
                        m_keep[j] = -1;
                    }
                }
            }

            // refine based on the slope of component i only if the
            // range of s is greater than a fraction 'min_range' of max
            // |s|. This eliminates components that consist of small
            // fluctuations on a constant slope background.
            if ((smax - smin) > m_min_range*ss) {
                // maximum allowable difference in slope between
                // adjacent points.
                doublereal dmax = m_curve*(smax - smin);
                for (size_t j = 0; j < n-2; j++) {
                    double dz0 = m_domain->dz(j);
                    double dz1 = m_domain->dz(j+1);
                    double r = fabs(m_s[j+1] - m_s[j]) / (dmax + m_thresh/dz0);
                    if (r > 1.0 && dz0 >= 2 * m_gridmin && dz1 >= 2 * m_gridmin) {
                        m_loc[j] = 1;
                        m_loc[j+1] = 1;
                        //m_c[name] = 1;
                        m_c[i] = 1;
                    }
                    if (r >= m_prune) {
                        m_keep[j+1] = 1;
                    } else if (m_keep[j+1] == 0) {
                        m_keep[j+1] = -1;
                    }
                }
            }
        }
    }

    StFlow* fflame = dynamic_cast<StFlow*>(m_domain);

    // Refine based on properties of the grid itself
    for (size_t j = 1; j < n-1; j++) {
        double dz = m_domain->dz(j);
        double dz1 = m_domain->dz(j-1);
        // Add a new point if the ratio with left interval is too large
        if (dz > m_ratio*dz1) {
            m_loc[j] = 1;
            //m_c[fmt::format("point {}", j)] = 1;
            m_keep[j-1] = 1;
            m_keep[j] = 1;
            m_keep[j+1] = 1;
            m_keep[j+2] = 1;
        }

        // Add a point if the ratio with right interval is too large
        if (dz < dz1/m_ratio) {
            m_loc[j-1] = 1;
            //m_c[fmt::format("point {}", j-1)] = 1;
            m_keep[j-2] = 1;
            m_keep[j-1] = 1;
            m_keep[j] = 1;
            m_keep[j+1] = 1;
        }

        // Keep the point if removing would make the ratio with the left
        // interval too large.
        if (j > 1 && m_domain->d2z(j) > m_ratio * m_domain->dz(j-2)) {
            m_keep[j] = 1;
        }

        // Keep the point if removing would make the ratio with the right
        // interval too large.
        if (j < n-2 && m_domain->d2z(j) > m_ratio * m_domain->dz(j+1)) {
            m_keep[j] = 1;
        }

        // Keep the point where the temperature is fixed
        if (fflame && fflame->domainType() == cFreeFlow && z[j] == fflame->m_zfixed) {
            m_keep[j] = 1;
        }
    }

    // Don't allow pruning to remove multiple adjacent grid points
    // in a single pass.
    for (size_t j = 2; j < n-1; j++) {
        if (m_keep[j] == -1 && m_keep[j-1] == -1) {
            m_keep[j] = 1;
        }
    }

    m_np = 0;
    for (size_t j = 0; j != n; j++) {
        m_np += m_loc[j];
    }
    return m_np;
}

double Refiner::value(const double* x, size_t i, size_t j)
{
    return x[m_domain->index(i,j)];
}

void Refiner::show()
{
    //if (!m_loc.empty()) {
    if ( m_np != 0 ) {
        writelog("\n\n");
        writeline('#', 78);
        writelog(string("Refining grid in ") +
                 m_domain->id()+".\n"
                 +"    New points inserted after grid points ");
        for (size_t j = 0; j != m_domain->nPoints(); j++) {
            if (newPointNeeded(j))
                writelog("{} ", j);
        }
        writelog("\n");
        writelog("    to resolve ");
        //for (const auto& c : m_c) {
        //    writelog(string(c.first)+" ");
        //}
        for (size_t j = 0; j != m_domain->nComponents(); j++) {
            if (m_c[j] ==1)
                writelog("{} ", m_domain->componentName(j));
        }
        writelog("\n");
        writeline('#', 78);
    } else if (m_domain->nPoints() > 1) {
        writelog("\n\n");
        writelog("no new points needed in "+m_domain->id()+"\n");
    }
}

int Refiner::getNewGrid(int n, const doublereal* z,
                        int nn, doublereal* zn)
{
    int nnew = static_cast<int>(m_loc.size());
    if (nnew + n > nn) {
        throw CanteraError("Refine::getNewGrid", "array size too small.");
    }

    if (m_loc.empty()) {
        copy(z, z + n, zn);
        return 0;
    }

    int jn = 0;
    for (int j = 0; j < n - 1; j++) {
        zn[jn] = z[j];
        jn++;
        if ( newPointNeeded(j) ) {
            zn[jn] = 0.5*(z[j] + z[j+1]);
            jn++;
        }
    }
    zn[jn] = z[n-1];
    return 0;
}
}
