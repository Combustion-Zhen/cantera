//! @file StFlow.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_STFLOW_H
#define CT_STFLOW_H

#include "Domain1D.h"
#include "cantera/base/Array.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/kinetics/Kinetics.h"

namespace Cantera
{

//------------------------------------------
//   constants
//------------------------------------------

// Offsets of solution components in the solution array.
const size_t c_offset_U = 0; // axial velocity
const size_t c_offset_V = 1; // strain rate
const size_t c_offset_T = 2; // temperature
const size_t c_offset_L = 3; // (1/r)dP/dr
const size_t c_offset_E = 4; // electric poisson's equation
const size_t c_offset_Y = 5; // mass fractions

class Transport;

/**
 *  This class represents 1D flow domains that satisfy the one-dimensional
 *  similarity solution for chemically-reacting, axisymmetric flows.
 *  @ingroup onedim
 */
class StFlow : public Domain1D
{
public:
    //--------------------------------
    // construction and destruction
    //--------------------------------

    //! Create a new flow domain.
    //! @param ph Object representing the gas phase. This object will be used
    //!     to evaluate all thermodynamic, kinetic, and transport properties.
    //! @param nsp Number of species.
    //! @param points Initial number of grid points
    StFlow(ThermoPhase* ph = 0, size_t nsp = 1, size_t points = 1);

    //! Delegating constructor
    StFlow(shared_ptr<ThermoPhase> th, size_t nsp = 1, size_t points = 1) :
        StFlow(th.get(), nsp, points) {
    }

    //! @name Problem Specification
    //! @{

    //! Change the grid size. Called after grid refinement.
    virtual void resize(size_t components, size_t points);

    virtual void setupGrid(size_t n, const doublereal* z);

    virtual void resetBadValues(double* xg);

    //! set the transport manager
    void setTransport(Transport& trans);

    //! Write the initial solution estimate into array x.
    virtual void _getInitialSoln(double* x);

    virtual void _finalize(const doublereal* x);

    virtual std::string componentName(size_t n) const;

    virtual size_t componentIndex(const std::string& name) const;

    // Zhen Lu 210917
    //! Return the type of flow domain being represented.
    //! @see setFreeFlow setAxisymmetricFlow
    virtual std::string flowType() const;

    // Zhen Lu 210920
    //! Set the Cartesian coordinates
    void setCartesian();

    //! Set the Cylindrical coordinates
    void setCylindrical();

    //! Set the Spherical coordinates
    void setSpherical();

    //! Print the solution.
    virtual void showSolution(const doublereal* x);

    virtual void restore(const XML_Node& dom, doublereal* soln,
                         int loglevel);

    //! Save the current solution for this domain into an XML_Node
    /*!
     *  @param o    XML_Node to save the solution to.
     *  @param sol  Current value of the solution vector. The object will pick
     *              out which part of the solution vector pertains to this
     *              object.
     *
     * @deprecated The XML output format is deprecated and will be removed in
     *     Cantera 3.0.
     */
    virtual XML_Node& save(XML_Node& o, const doublereal* const sol);

    void solveEnergyEqn(size_t j=npos);

    void fixTemperature(size_t j=npos);

    //! Set the emissivities for the boundary values
    /*!
     * Reads the emissivities for the left and right boundary values in the
     * radiative term and writes them into the variables, which are used for the
     * calculation.
     */
    void setBoundaryEmissivities(double e_left, double e_right);

    //! Set the gas object state to be consistent with the solution at point j.
    void setGas(const doublereal* x, size_t j);

    //! Set the gas state to be consistent with the solution at the midpoint
    //! between j and j + 1.
    void setGasAtMidpoint(const doublereal* x, size_t j);

    inline ThermoPhase& phase() {
        return *m_thermo;
    }
    inline Kinetics& kinetics() {
        return *m_kin;
    }

    /**
     * Set the thermo manager. Note that the flow equations assume
     * the ideal gas equation.
     */
    inline void setThermo(IdealGasPhase& th) {
        m_thermo = &th;
    }

    //! Set the kinetics manager. The kinetics manager must
    inline void setKinetics(Kinetics& kin) {
        m_kin = &kin;
    }

    //! Enable thermal diffusion, also known as Soret diffusion.
    //! Requires that multicomponent transport properties be
    //! enabled to carry out calculations.
    inline void enableSoret(bool withSoret) {
        m_do_soret = withSoret;
    }
    inline bool withSoret() const {
        return m_do_soret;
    }

    //! Set the pressure. Since the flow equations are for the limit of small
    //! Mach number, the pressure is very nearly constant throughout the flow.
    inline void setPressure(doublereal p) {
        m_press = p;
    }

    //! The current pressure [Pa].
    inline doublereal pressure() const {
        return m_press;
    }

    //! Sometimes it is desired to carry out the simulation using a specified
    //! temperature profile, rather than computing it by solving the energy
    //! equation. This method specifies this profile.
    void setFixedTempProfile(vector_fp& zfixed, vector_fp& tfixed) {
        m_zfix = zfixed;
        m_tfix = tfixed;
    }

    /*!
     * Set the temperature fixed point at grid point j, and disable the energy
     * equation so that the solution will be held to this value.
     */
    void setTemperature(size_t j, doublereal t) {
        m_fixedtemp[j] = t;
        m_do_energy[j] = false;
    }

    //! The fixed temperature value at point j.
    inline doublereal T_fixed(size_t j) const {
        return m_fixedtemp[j];
    }

    // @}

    //! Set flow configuration for freely-propagating flames, using an internal
    //! point with a fixed temperature as the condition to determine the inlet
    //! mass flux.
    void setFreeFlow() {
        m_type = cFreeFlow;
        m_dovisc = false;
    }

    //! Set flow configuration for axisymmetric counterflow or burner-stabilized
    //! flames, using specified inlet mass fluxes.
    void setAxisymmetricFlow() {
        m_type = cAxisymmetricStagnationFlow;
        m_dovisc = true;
    }

    // Zhen Lu 210917
    //! Set flow configuration for one-dimensional radial propagating flames, 
    //! using specific inlet mass fluxes
    void setRadialFlow() {
        m_type = cRadialFlow;
        m_dovisc = false;
    }

    //! Set flow configuration for radial stagnation flames, 
    //! using specific inlet mass fluxes
    void setTubularFlow() {
        m_type = cTubularFlow;
        m_dovisc = true;
    }

    //! Set steady state
    void setSteady() {
        m_ttype = cSteady;
    }

    //! Set transient state
    void setTransient() {
        m_ttype = cTransient;
    }

    //! Turn radiation on / off.
    /*!
     *  The simple radiation model used was established by Y. Liu and B. Rogg
     *  [Y. Liu and B. Rogg, Modelling of thermally radiating diffusion flames
     *  with detailed chemistry and transport, EUROTHERM Seminars, 17:114-127,
     *  1991]. This model considers the radiation of CO2 and H2O.
     */
    void enableRadiation(bool doRadiation) {
        m_do_radiation = doRadiation;
    }

    // not used
    //! Returns `true` if the radiation term in the energy equation is enabled
    bool radiationEnabled() const {
        return m_do_radiation;
    }

    // not used
    //! Return radiative heat loss at grid point j
    double radiativeHeatLoss(size_t j) const {
        return m_qdotRadiation[j];
    }

    // not used
    //! Return emissivitiy at left boundary
    double leftEmissivity() const { return m_epsilon_left; }

    // not used
    //! Return emissivitiy at right boundary
    double rightEmissivity() const { return m_epsilon_right; }

    inline bool doEnergy(size_t j) const {
        return m_do_energy[j];
    }

    inline doublereal density(size_t j) const {
        return m_rho[j];
    }

    inline virtual bool fixed_mdot() const {
        return (domainType() != cFreeFlow);
    }

    // in Outlet1D::init
    void setViscosityFlag(bool dovisc) {
        m_dovisc = dovisc;
    }

    //! Index of the species on the left boundary with the largest mass fraction
    inline size_t leftExcessSpecies() const {
        return m_kExcessLeft;
    }

    //! Index of the species on the right boundary with the largest mass fraction
    inline size_t rightExcessSpecies() const {
        return m_kExcessRight;
    }

protected:
    /*!
     *  Evaluate the residual function for axisymmetric stagnation flow. If
     *  j == npos, the residual function is evaluated at all grid points.
     *  Otherwise, the residual function is only evaluated at grid points
     *  j-1, j, and j+1. This option is used to efficiently evaluate the
     *  Jacobian numerically.
     */
    virtual void eval(size_t j, doublereal* x, doublereal* r,
                      integer* mask, doublereal rdt);

    //! Evaluate the residual function. This function is called in eval
    //! after updateProperties is called.
    virtual void evalResidual(double* x, double* rsd, int* diag,
                              double rdt, size_t jmin, size_t jmax);

    // Zhen Lu 210920
    //! Evaluate all residual components at the left boundary.
    virtual void evalLeftBoundary(double* x, double* res, int* diag,
                                  double rdt);

    //! Evaluate all residual components at the right boundary.
    virtual void evalRightBoundary(double* x, double* res, int* diag,
                                   double rdt);

    //! Evaluate the residual corresponding to the continuity equation at all
    //! interior grid points.
    virtual void evalContinuity(size_t j, double* x, double* r,
                                int* diag, double rdt);

    //! Evaluate the residual corresponding to the V equation at all
    //! interior grid points.
    virtual void evalRadialMomentum(size_t j, double* x, double* r,
                                    int* diag, double rdt);

    //! Evaluate the residual corresponding to the species equation at all
    //! interior grid points.
    virtual void evalSpecies(size_t j, double* x, double* r,
                             int* diag, double rdt);

    //! Evaluate the residual corresponding to the energy equation at all
    //! interior grid points.
    virtual void evalEnergy(size_t j, double* x, double* r,
                            int* diag, double rdt);

    //! Update the properties (thermo, transport, and diffusion flux).
    //! This function is called in eval after the points which need
    //! to be updated are defined.
    virtual void updateProperties(size_t jg, double* x, size_t jmin, size_t jmax);

    /**
     * Update the thermodynamic properties from point j0 to point j1
     * (inclusive), based on solution x.
     */
    void updateThermo(const doublereal* x, size_t j0, size_t j1);

    //! Update the transport properties at grid points in the range from `j0`
    //! to `j1`, based on solution `x`.
    virtual void updateTransport(doublereal* x, size_t j0, size_t j1);

    //! Update the diffusive mass fluxes.
    virtual void updateDiffFluxes(const doublereal* x, size_t j0, size_t j1);

    //! Calculate divergence of the stress
    doublereal shear(const doublereal* x, size_t j) const;

    //! Calculate divergence of the diffusive flux
    doublereal divDiffFlux(size_t k, size_t j) const;

    //! Calculate divergence of the heat flux
    doublereal divHeatFlux(const doublereal* x, size_t j) const;

    //! @name convective spatial derivatives.
    //! These use upwind differencing, assuming u(z) is negative
    //! @{
    doublereal dVdz(const doublereal* x, size_t j) const;

    doublereal dYdz(const doublereal* x, size_t k, size_t j) const;

    doublereal dTdz(const doublereal* x, size_t j) const;
    //! @}

    //! Write the net production rates at point `j` into array `m_wdot`
    void getWdot(doublereal* x, size_t j) {
        setGas(x,j);
        m_kin->getNetProductionRates(&m_wdot(0,j));
    }

    //! @name Solution components
    //! @{

    inline doublereal wdot(size_t k, size_t j) const {
        return m_wdot(k,j);
    }

    inline doublereal T(const doublereal* x, size_t j) const {
        return x[index(c_offset_T, j)];
    }
    inline doublereal& T(doublereal* x, size_t j) {
        return x[index(c_offset_T, j)];
    }
    inline doublereal T_prev(size_t j) const {
        return prevSoln(c_offset_T, j);
    }

    inline doublereal rho_u(const doublereal* x, size_t j) const {
        return m_rho[j]*x[index(c_offset_U, j)];
    }

    inline doublereal u(const doublereal* x, size_t j) const {
        return x[index(c_offset_U, j)];
    }

    inline doublereal V(const doublereal* x, size_t j) const {
        return x[index(c_offset_V, j)];
    }
    inline doublereal V_prev(size_t j) const {
        return prevSoln(c_offset_V, j);
    }

    inline doublereal lambda(const doublereal* x, size_t j) const {
        return x[index(c_offset_L, j)];
    }

    inline doublereal Y(const doublereal* x, size_t k, size_t j) const {
        return x[index(c_offset_Y + k, j)];
    }
    inline doublereal& Y(doublereal* x, size_t k, size_t j) {
        return x[index(c_offset_Y + k, j)];
    }
    inline doublereal Y_prev(size_t k, size_t j) const {
        return prevSoln(c_offset_Y + k, j);
    }

    inline doublereal X(const doublereal* x, size_t k, size_t j) const {
        return m_wtm[j]*Y(x,k,j)/m_wt[k];
    }

    //! @}

    inline size_t mindex(size_t k, size_t j, size_t m) const {
        return m*m_nsp*m_nsp + m_nsp*j + k;
    }

    //---------------------------------------------------------
    //             member data
    //---------------------------------------------------------

    doublereal m_press; // pressure

    // mixture thermo properties
    vector_fp m_rho;
    vector_fp m_wtm;

    // species thermo properties
    vector_fp m_wt;
    vector_fp m_cp;

    // transport properties
    vector_fp m_visc;
    vector_fp m_tcon;
    vector_fp m_diff;
    vector_fp m_multidiff;
    Array2D m_dthermal;
    Array2D m_flux;

    // production rates
    Array2D m_wdot;

    size_t m_nsp;

    IdealGasPhase* m_thermo;
    Kinetics* m_kin;
    Transport* m_trans;

    // boundary emissivities for the radiation calculations
    doublereal m_epsilon_left;
    doublereal m_epsilon_right;

    //! Indices within the ThermoPhase of the radiating species. First index is
    //! for CO2, second is for H2O.
    std::vector<size_t> m_kRadiating;

    // flags
    std::vector<bool> m_do_energy;
    bool m_do_soret;
    std::vector<bool> m_do_species;
    bool m_do_multicomponent;

    //! flag for the radiative heat loss
    bool m_do_radiation;

    //! radiative heat loss vector
    vector_fp m_qdotRadiation;

    // fixed T and Y values
    vector_fp m_fixedtemp;
    vector_fp m_zfix;
    vector_fp m_tfix;

    //! Index of species with a large mass fraction at each boundary, for which
    //! the mass fraction may be calculated as 1 minus the sum of the other mass
    //! fractions
    size_t m_kExcessLeft;
    size_t m_kExcessRight;

    bool m_dovisc;

public:
    //! Location of the point where temperature is fixed
    double m_zfixed;

    //! Temperature at the point used to fix the flame location
    double m_tfixed;

private:
    vector_fp m_ybar;
};

}

#endif
