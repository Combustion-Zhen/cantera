//! @file StFlow.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/oneD/StFlow.h"
#include "cantera/oneD/refine.h"
#include "cantera/base/ctml.h"
#include "cantera/transport/TransportBase.h"
#include "cantera/numerics/funcs.h"
#include "cantera/base/global.h"
#include "cantera/zerodim.h"
#include "cantera/base/Solution.h"

#include <set>

using namespace std;

namespace Cantera
{

// public member functions

StFlow::StFlow(ThermoPhase* ph, size_t nsp, size_t points) :
    Domain1D(nsp+c_offset_Y, points),
    m_press(-1.0),
    m_nsp(nsp),
    m_thermo(0), m_kin(0), m_trans(0),
    m_epsilon_left(0.0), m_epsilon_right(0.0),
    m_do_soret(false),
    m_do_multicomponent(false),
    m_do_radiation(false),
    m_do_reaction(true),
    m_kExcessLeft(0), m_kExcessRight(0),
    m_zfixed(Undef), m_tfixed(-1.),
    m_do_ignition(false),
    m_ign_energy(0.0), m_ign_radius(0.0), m_ign_time(0.0)
{
    if (ph->type() == "IdealGas") {
        m_thermo = static_cast<IdealGasPhase*>(ph);
    } else {
        throw CanteraError("StFlow::StFlow",
                           "Unsupported phase type: need 'IdealGasPhase'");
    }
    m_type = cFlowType;
    m_points = points;

    if (ph == 0) {
        return; // used to create a dummy object
    }

    size_t nsp2 = m_thermo->nSpecies();
    if (nsp2 != m_nsp) {
        m_nsp = nsp2;
        Domain1D::resize(m_nsp+c_offset_Y, points);
    }

    // make a local copy of the species molecular weight vector
    m_wt = m_thermo->molecularWeights();

    // the species mass fractions are the last components in the solution
    // vector, so the total number of components is the number of species
    // plus the offset of the first mass fraction.
    m_nv = m_nsp + c_offset_Y;

    m_nc = m_nsp + 1;
    m_maxScalar.resize(m_nc, 0.0);
    m_minScalar.resize(m_nc, 0.0);

    // enable all species equations by default
    m_do_species.resize(m_nsp, true);

    // but turn off the energy equation at all points
    m_do_energy.resize(m_points,false);

    m_diff.resize(m_nsp*m_points);
    m_multidiff.resize(m_nsp*m_nsp*m_points);
    m_flux.resize(m_nsp,m_points);
    m_wdot.resize(m_nsp,m_points, 0.0);
    m_ybar.resize(m_nsp);
    m_qdotRadiation.resize(m_points, 0.0);

    //-------------- default solution bounds --------------------
    setBounds(0, -1e20, 1e20); // no bounds on u
    setBounds(1, -1e20, 1e20); // V
    setBounds(2, 200.0, 2*m_thermo->maxTemp()); // temperature bounds
    setBounds(3, -1e20, 1e20); // lambda should be negative

    // mass fraction bounds
    for (size_t k = 0; k < m_nsp; k++) {
        setBounds(c_offset_Y+k, -1.0e-7, 1.0e5);
    }

    setBoundsScalar(0, 200.0, 2*m_thermo->maxTemp()); // temperature bounds
    for (size_t k = 0; k < m_nsp; k++) {
        setBoundsScalar(cOffsetScalarY+k, -1.0e-7, 1.0e5);
    }

    //-------------------- grid refinement -------------------------
    m_refiner->setActive(c_offset_U, false);
    m_refiner->setActive(c_offset_V, false);
    m_refiner->setActive(c_offset_T, false);
    m_refiner->setActive(c_offset_L, false);

    vector_fp gr;
    for (size_t ng = 0; ng < m_points; ng++) {
        gr.push_back(1.0*ng/m_points);
    }
    setupGrid(m_points, gr.data());

    // Find indices for radiating species
    m_kRadiating.resize(2, npos);
    m_kRadiating[0] = m_thermo->speciesIndex("CO2");
    m_kRadiating[1] = m_thermo->speciesIndex("H2O");
}

void StFlow::resize(size_t ncomponents, size_t points)
{
    Domain1D::resize(ncomponents, points);
    m_rho.resize(m_points, 0.0);
    m_wtm.resize(m_points, 0.0);
    m_cp.resize(m_points, 0.0);
    m_visc.resize(m_points, 0.0);
    m_tcon.resize(m_points, 0.0);
    // Zhen Lu 211005
    m_rho_last.resize(m_points, 0.0);

    m_diff.resize(m_nsp*m_points);
    if (m_do_multicomponent) {
        m_multidiff.resize(m_nsp*m_nsp*m_points);
        m_dthermal.resize(m_nsp, m_points, 0.0);
    }
    m_flux.resize(m_nsp,m_points);
    m_wdot.resize(m_nsp,m_points, 0.0);
    m_do_energy.resize(m_points,false);
    m_qdotRadiation.resize(m_points, 0.0);
    m_fixedtemp.resize(m_points);

    m_z.resize(m_points);
}

void StFlow::setupGrid(size_t n, const doublereal* z)
{
    resize(m_nv, n);

    m_z[0] = z[0];
    for (size_t j = 1; j < m_points; j++) {
        if (z[j] <= z[j-1]) {
            throw CanteraError("StFlow::setupGrid",
                               "grid points must be monotonically increasing");
        }
        m_z[j] = z[j];
    }
}

void StFlow::resetBadValues(double* xg)
{
    double* x = xg + loc();
    for (size_t j = 0; j < m_points; j++) {
        double* Y = x + m_nv*j + c_offset_Y;
        m_thermo->setMassFractions(Y);
        m_thermo->getMassFractions(Y);
    }
}

void StFlow::setTransport(Transport& trans)
{
    m_trans = &trans;
    m_do_multicomponent = (m_trans->transportType() == "Multi" || m_trans->transportType() == "CK_Multi");

    m_diff.resize(m_nsp*m_points);
    if (m_do_multicomponent) {
        m_multidiff.resize(m_nsp*m_nsp*m_points);
        m_dthermal.resize(m_nsp, m_points, 0.0);
    }
}

void StFlow::_getInitialSoln(double* x)
{
    for (size_t j = 0; j < m_points; j++) {
        T(x,j) = m_thermo->temperature();
        m_thermo->getMassFractions(&Y(x, 0, j));
    }
}

void StFlow::_finalize(const doublereal* x)
{
    if (!m_do_multicomponent && m_do_soret) {
        throw CanteraError("StFlow::_finalize",
            "Thermal diffusion (the Soret effect) is enabled, and requires "
            "using a multicomponent transport model.");
    }

    size_t nz = m_zfix.size();
    bool e = m_do_energy[0];
    for (size_t j = 0; j < m_points; j++) {
        if (e || nz == 0) {
            m_fixedtemp[j] = T(x, j);
        } else {
            double zz = (z(j) - z(0))/(z(m_points - 1) - z(0));
            double tt = linearInterp(zz, m_zfix, m_tfix);
            m_fixedtemp[j] = tt;
        }
    }
    if (e) {
        solveEnergyEqn();
    }

    if (domainType() == cFreeFlow) {
        // If the domain contains the temperature fixed point, make sure that it
        // is correctly set. This may be necessary when the grid has been modified
        // externally.
        if (m_tfixed != Undef) {
            for (size_t j = 0; j < m_points; j++) {
                if (z(j) == m_zfixed) {
                    return; // fixed point is already set correctly
                }
            }

            for (size_t j = 0; j < m_points - 1; j++) {
                // Find where the temperature profile crosses the current
                // fixed temperature.
                if ((T(x, j) - m_tfixed) * (T(x, j+1) - m_tfixed) <= 0.0) {
                    m_tfixed = T(x, j+1);
                    m_zfixed = z(j+1);
                    return;
                }
            }
        }
    }
}

string StFlow::componentName(size_t n) const
{
    switch (n) {
    case 0:
        return "velocity";
    case 1:
        return "spread_rate";
    case 2:
        return "T";
    case 3:
        return "lambda";
    case 4:
        return "eField";
    default:
        if (n >= c_offset_Y && n < (c_offset_Y + m_nsp)) {
            return m_thermo->speciesName(n - c_offset_Y);
        } else {
            return "<unknown>";
        }
    }
}

size_t StFlow::componentIndex(const std::string& name) const
{
    if (name=="velocity") {
        return 0;
    } else if (name=="spread_rate") {
        return 1;
    } else if (name=="T") {
        return 2;
    } else if (name=="lambda") {
        return 3;
    } else if (name == "eField") {
        return 4;
    } else {
        for (size_t n=c_offset_Y; n<m_nsp+c_offset_Y; n++) {
            if (componentName(n)==name) {
                return n;
            }
        }
        throw CanteraError("StFlow1D::componentIndex",
                           "no component named " + name);
    }
}

bool StFlow::componentActive(size_t n) const
{
    switch (n) {
    case c_offset_V: // spread_rate
        return m_type != cFreeFlow;
    case c_offset_L: // lambda
        return m_type != cFreeFlow;
    case c_offset_E: // eField
        return false;
    default:
        return true;
    }
}

// Zhen Lu 210917
string StFlow::flowType() const
{
    switch (m_type)
    {
    case cFreeFlow:
        return "Free Flame";
    case cAxisymmetricStagnationFlow:
        return "Axisymmetric Stagnation";
    case cRadialFlow:
        return "Radial Flame";
    case cPolarFlow:
        return "Free Polar Flame";
    case cTubularFlow:
        return "Tubular Flame";
    default:
        throw CanteraError("StFlow::flowType", "Unknown value for 'm_type'");
    }
}

// Zhen Lu 210920
void StFlow::setCartesian()
{
    if (m_type == cTubularFlow)
        throw CanteraError
        (
            "StFlow::setCartesian",
            "TubularFlow works in the Cylindrical coordinates"
        );

    m_ctype = cCartesian;
}

void StFlow::setCylindrical()
{
    if (m_type == cAxisymmetricStagnationFlow)
        throw CanteraError
        (
            "StFlow::setCylindrical",
            "AxisymmetricStagnationFlow works in the Cartesian coordinates"
        );

    m_ctype = cCylindrical;
}

void StFlow::setSpherical()
{
    if (m_type == cTubularFlow)
        throw CanteraError
        (
            "StFlow::setSpherical",
            "TubularFlow works in the Cylindrical coordinates"
        );
    else if (m_type == cAxisymmetricStagnationFlow)
        throw CanteraError
        (
            "StFlow::setSpherical",
            "AxisymmetricStagnationFlow works in the Cartesian coordinates"
        );

    m_ctype = cSpherical;
}

void StFlow::initTimeInteg(doublereal dt, doublereal* x0)
{
    Domain1D::initTimeInteg(dt, x0);
    updateThermo(x0+loc(), 0, m_points-1);
    std::copy(m_rho.begin(), m_rho.end(), m_rho_last.begin());
}

void StFlow::showSolution(const doublereal* x)
{
    writelog("    Pressure:  {:10.4g} Pa\n", m_press);

    Domain1D::showSolution(x);

    if (m_do_radiation) {
        writeline('-', 79, false, true);
        writelog("\n          z      radiative heat loss");
        writeline('-', 79, false, true);
        for (size_t j = 0; j < m_points; j++) {
            writelog("\n {:10.4g}        {:10.4g}", m_z[j], m_qdotRadiation[j]);
        }
        writelog("\n");
    }
}

void StFlow::restore(const XML_Node& dom, doublereal* soln, int loglevel)
{
    Domain1D::restore(dom, soln, loglevel);
    vector<string> ignored;
    size_t nsp = m_thermo->nSpecies();
    vector_int did_species(nsp, 0);

    vector<XML_Node*> str = dom.getChildren("string");
    for (size_t istr = 0; istr < str.size(); istr++) {
        const XML_Node& nd = *str[istr];
        writelog(nd["title"]+": "+nd.value()+"\n");
    }

    double pp = getFloat(dom, "pressure", "pressure");
    setPressure(pp);
    vector<XML_Node*> d = dom.child("grid_data").getChildren("floatArray");
    vector_fp x;
    size_t np = 0;
    bool readgrid = false, wrote_header = false;
    for (size_t n = 0; n < d.size(); n++) {
        const XML_Node& fa = *d[n];
        string nm = fa["title"];
        if (nm == "z") {
            getFloatArray(fa,x,false);
            np = x.size();
            if (loglevel >= 2) {
                writelog("Grid contains {} points.\n", np);
            }
            readgrid = true;
            setupGrid(np, x.data());
        }
    }
    if (!readgrid) {
        throw CanteraError("StFlow::restore",
                           "domain contains no grid points.");
    }

    debuglog("Importing datasets:\n", loglevel >= 2);
    for (size_t n = 0; n < d.size(); n++) {
        const XML_Node& fa = *d[n];
        string nm = fa["title"];
        getFloatArray(fa,x,false);
        if (nm == "u") {
            debuglog("axial velocity   ", loglevel >= 2);
            if (x.size() != np) {
                throw CanteraError("StFlow::restore",
                                   "axial velocity array size error");
            }
            for (size_t j = 0; j < np; j++) {
                soln[index(c_offset_U,j)] = x[j];
            }
        } else if (nm == "z") {
            ; // already read grid
        } else if (nm == "V") {
            debuglog("radial velocity   ", loglevel >= 2);
            if (x.size() != np) {
                throw CanteraError("StFlow::restore",
                                   "radial velocity array size error");
            }
            for (size_t j = 0; j < np; j++) {
                soln[index(c_offset_V,j)] = x[j];
            }
        } else if (nm == "T") {
            debuglog("temperature   ", loglevel >= 2);
            if (x.size() != np) {
                throw CanteraError("StFlow::restore",
                                   "temperature array size error");
            }
            for (size_t j = 0; j < np; j++) {
                soln[index(c_offset_T,j)] = x[j];
            }

            // For fixed-temperature simulations, use the imported temperature
            // profile by default.  If this is not desired, call
            // setFixedTempProfile *after* restoring the solution.
            vector_fp zz(np);
            for (size_t jj = 0; jj < np; jj++) {
                zz[jj] = (grid(jj) - zmin())/(zmax() - zmin());
            }
            setFixedTempProfile(zz, x);
        } else if (nm == "L") {
            debuglog("lambda   ", loglevel >= 2);
            if (x.size() != np) {
                throw CanteraError("StFlow::restore",
                                   "lambda array size error");
            }
            for (size_t j = 0; j < np; j++) {
                soln[index(c_offset_L,j)] = x[j];
            }
        } else if (m_thermo->speciesIndex(nm) != npos) {
            debuglog(nm+"   ", loglevel >= 2);
            if (x.size() == np) {
                size_t k = m_thermo->speciesIndex(nm);
                did_species[k] = 1;
                for (size_t j = 0; j < np; j++) {
                    soln[index(k+c_offset_Y,j)] = x[j];
                }
            }
        } else {
            ignored.push_back(nm);
        }
    }

    if (loglevel >=2 && !ignored.empty()) {
        writelog("\n\n");
        writelog("Ignoring datasets:\n");
        size_t nn = ignored.size();
        for (size_t n = 0; n < nn; n++) {
            writelog(ignored[n]+"   ");
        }
    }

    if (loglevel >= 1) {
        for (size_t ks = 0; ks < nsp; ks++) {
            if (did_species[ks] == 0) {
                if (!wrote_header) {
                    writelog("Missing data for species:\n");
                    wrote_header = true;
                }
                writelog(m_thermo->speciesName(ks)+" ");
            }
        }
    }

    if (dom.hasChild("energy_enabled")) {
        getFloatArray(dom, x, false, "", "energy_enabled");
        if (x.size() == nPoints()) {
            for (size_t i = 0; i < x.size(); i++) {
                m_do_energy[i] = (x[i] != 0);
            }
        } else if (!x.empty()) {
            throw CanteraError("StFlow::restore", "energy_enabled is length {}"
                               "but should be length {}", x.size(), nPoints());
        }
    }

    if (dom.hasChild("species_enabled")) {
        getFloatArray(dom, x, false, "", "species_enabled");
        if (x.size() == m_nsp) {
            for (size_t i = 0; i < x.size(); i++) {
                m_do_species[i] = (x[i] != 0);
            }
        } else if (!x.empty()) {
            // This may occur when restoring from a mechanism with a different
            // number of species.
            if (loglevel > 0) {
                warn_user("StFlow::restore", "species_enabled is "
                    "length {} but should be length {}. Enabling all species "
                    "equations by default.", x.size(), m_nsp);
            }
            m_do_species.assign(m_nsp, true);
        }
    }

    if (dom.hasChild("refine_criteria")) {
        XML_Node& ref = dom.child("refine_criteria");
        refiner().setCriteria(getFloat(ref, "ratio"), getFloat(ref, "slope"),
                              getFloat(ref, "curve"), getFloat(ref, "prune"));
        refiner().setGridMin(getFloat(ref, "grid_min"));
    }

    if (domainType() == cFreeFlow) {
        getOptionalFloat(dom, "t_fixed", m_tfixed);
        getOptionalFloat(dom, "z_fixed", m_zfixed);
    }
}

XML_Node& StFlow::save(XML_Node& o, const doublereal* const sol)
{
    Array2D soln(m_nv, m_points, sol + loc());
    XML_Node& flow = Domain1D::save(o, sol);
    flow.addAttribute("type",flowType());

    XML_Node& gv = flow.addChild("grid_data");
    addFloat(flow, "pressure", m_press, "Pa", "pressure");

    addFloatArray(gv,"z",m_z.size(), m_z.data(),
                  "m","length");
    vector_fp x(soln.nColumns());

    soln.getRow(c_offset_U, x.data());
    addFloatArray(gv,"u",x.size(),x.data(),"m/s","velocity");

    soln.getRow(c_offset_V, x.data());
    addFloatArray(gv,"V",x.size(),x.data(),"1/s","rate");

    soln.getRow(c_offset_T, x.data());
    addFloatArray(gv,"T",x.size(),x.data(),"K","temperature");

    soln.getRow(c_offset_L, x.data());
    addFloatArray(gv,"L",x.size(),x.data(),"N/m^4");

    for (size_t k = 0; k < m_nsp; k++) {
        soln.getRow(c_offset_Y+k, x.data());
        addFloatArray(gv,m_thermo->speciesName(k),
                      x.size(),x.data(),"","massFraction");
    }
    if (m_do_radiation) {
        addFloatArray(gv, "radiative_heat_loss", m_z.size(),
            m_qdotRadiation.data(), "W/m^3", "specificPower");
    }
    vector_fp values(nPoints());
    for (size_t i = 0; i < nPoints(); i++) {
        values[i] = m_do_energy[i];
    }
    addNamedFloatArray(flow, "energy_enabled", nPoints(), &values[0]);

    values.resize(m_nsp);
    for (size_t i = 0; i < m_nsp; i++) {
        values[i] = m_do_species[i];
    }
    addNamedFloatArray(flow, "species_enabled", m_nsp, &values[0]);

    XML_Node& ref = flow.addChild("refine_criteria");
    addFloat(ref, "ratio", refiner().maxRatio());
    addFloat(ref, "slope", refiner().maxDelta());
    addFloat(ref, "curve", refiner().maxSlope());
    addFloat(ref, "prune", refiner().prune());
    addFloat(ref, "grid_min", refiner().gridMin());
    if (m_zfixed != Undef) {
        addFloat(flow, "z_fixed", m_zfixed, "m");
        addFloat(flow, "t_fixed", m_tfixed, "K");
    }
    return flow;
}

void StFlow::restore(const AnyMap& state, double* soln, int loglevel)
{
    Domain1D::restore(state, soln, loglevel);
    m_press = state["pressure"].asDouble();
    setupGrid(nPoints(), state["grid"].asVector<double>(nPoints()).data());

    for (size_t i = 0; i < nComponents(); i++) {
        if (!componentActive(i)) {
            continue;
        }
        std::string name = componentName(i);
        if (state.hasKey(name)) {
            const vector_fp& data = state[name].asVector<double>(nPoints());
            for (size_t j = 0; j < nPoints(); j++) {
                soln[index(i,j)] = data[j];
            }
        } else if (loglevel) {
            warn_user("StFlow::restore", "Saved state does not contain values for "
                "component '{}' in domain '{}'.", name, id());
        }
    }

    if (state.hasKey("energy-enabled")) {
        const AnyValue& ee = state["energy-enabled"];
        if (ee.isScalar()) {
            m_do_energy.assign(nPoints(), ee.asBool());
        } else {
            m_do_energy = ee.asVector<bool>(nPoints());
        }
    }

    if (state.hasKey("Soret-enabled")) {
        m_do_soret = state["Soret-enabled"].asBool();
    }

    if (state.hasKey("species-enabled")) {
        const AnyValue& se = state["species-enabled"];
        if (se.isScalar()) {
            m_do_species.assign(m_thermo->nSpecies(), se.asBool());
        } else {
            m_do_species = se.asVector<bool>(m_thermo->nSpecies());
        }
    }

    if (state.hasKey("radiation-enabled")) {
        m_do_radiation = state["radiation-enabled"].asBool();
        if (m_do_radiation) {
            m_epsilon_left = state["emissivity-left"].asDouble();
            m_epsilon_right = state["emissivity-right"].asDouble();
        }
    }

    if (state.hasKey("refine-criteria")) {
        const AnyMap& criteria = state["refine-criteria"].as<AnyMap>();
        double ratio = criteria.getDouble("ratio", m_refiner->maxRatio());
        double slope = criteria.getDouble("slope", m_refiner->maxDelta());
        double curve = criteria.getDouble("curve", m_refiner->maxSlope());
        double prune = criteria.getDouble("prune", m_refiner->prune());
        m_refiner->setCriteria(ratio, slope, curve, prune);

        if (criteria.hasKey("grid-min")) {
            m_refiner->setGridMin(criteria["grid-min"].asDouble());
        }
        if (criteria.hasKey("max-points")) {
            m_refiner->setMaxPoints(criteria["max-points"].asInt());
        }
    }

    if (state.hasKey("fixed-point")) {
        m_zfixed = state["fixed-point"]["location"].asDouble();
        m_tfixed = state["fixed-point"]["temperature"].asDouble();
    }
}

AnyMap StFlow::serialize(const double* soln) const
{
    AnyMap state = Domain1D::serialize(soln);
    state["type"] = flowType();
    state["pressure"] = m_press;

    state["phase"]["name"] = m_thermo->name();
    AnyValue source = m_thermo->input().getMetadata("filename");
    state["phase"]["source"] = source.empty() ? "<unknown>" : source.asString();

    state["radiation-enabled"] = m_do_radiation;
    if (m_do_radiation) {
        state["radiative-heat-loss"] = m_qdotRadiation;
        state["emissivity-left"] = m_epsilon_left;
        state["emissivity-right"] = m_epsilon_right;
    }

    std::set<bool> energy_flags(m_do_energy.begin(), m_do_energy.end());
    if (energy_flags.size() == 1) {
        state["energy-enabled"] = m_do_energy[0];
    } else {
        state["energy-enabled"] = m_do_energy;
    }

    state["Soret-enabled"] = m_do_soret;

    std::set<bool> species_flags(m_do_species.begin(), m_do_species.end());
    if (species_flags.size() == 1) {
        state["species-enabled"] = m_do_species[0];
    } else {
        for (size_t k = 0; k < m_nsp; k++) {
            state["species-enabled"][m_thermo->speciesName(k)] = m_do_species[k];
        }
    }

    state["refine-criteria"]["ratio"] = m_refiner->maxRatio();
    state["refine-criteria"]["slope"] = m_refiner->maxDelta();
    state["refine-criteria"]["curve"] = m_refiner->maxSlope();
    state["refine-criteria"]["prune"] = m_refiner->prune();
    state["refine-criteria"]["grid-min"] = m_refiner->gridMin();
    state["refine-criteria"]["max-points"] =
        static_cast<long int>(m_refiner->maxPoints());

    if (m_zfixed != Undef) {
        state["fixed-point"]["location"] = m_zfixed;
        state["fixed-point"]["temperature"] = m_tfixed;
    }

    state["grid"] = m_z;
    vector_fp data(nPoints());
    for (size_t i = 0; i < nComponents(); i++) {
        if (componentActive(i)) {
            for (size_t j = 0; j < nPoints(); j++) {
                data[j] = soln[index(i,j)];
            }
            state[componentName(i)] = data;
        }
    }

    return state;
}

void StFlow::solveEnergyEqn(size_t j)
{
    bool changed = false;
    if (j == npos) {
        for (size_t i = 0; i < m_points; i++) {
            if (!m_do_energy[i]) {
                changed = true;
            }
            m_do_energy[i] = true;
        }
    } else {
        if (!m_do_energy[j]) {
            changed = true;
        }
        m_do_energy[j] = true;
    }
    m_refiner->setActive(c_offset_U, true);
    if (domainType()==cAxisymmetricStagnationFlow)
        m_refiner->setActive(c_offset_V, true);
    m_refiner->setActive(c_offset_T, true);
    if (changed) {
        needJacUpdate();
    }
}

void StFlow::fixTemperature(size_t j)
{
    bool changed = false;
    if (j == npos) {
        for (size_t i = 0; i < m_points; i++) {
            if (m_do_energy[i]) {
                changed = true;
            }
            m_do_energy[i] = false;
        }
    } else {
        if (m_do_energy[j]) {
            changed = true;
        }
        m_do_energy[j] = false;
    }
    m_refiner->setActive(c_offset_U, false);
    m_refiner->setActive(c_offset_V, false);
    m_refiner->setActive(c_offset_T, false);
    if (changed) {
        needJacUpdate();
    }
}

void StFlow::setBoundaryEmissivities(doublereal e_left, doublereal e_right)
{
    if (e_left < 0 || e_left > 1) {
        throw CanteraError("StFlow::setBoundaryEmissivities",
            "The left boundary emissivity must be between 0.0 and 1.0!");
    } else if (e_right < 0 || e_right > 1) {
        throw CanteraError("StFlow::setBoundaryEmissivities",
            "The right boundary emissivity must be between 0.0 and 1.0!");
    } else {
        m_epsilon_left = e_left;
        m_epsilon_right = e_right;
    }
}

void StFlow::setGas(const doublereal* x, size_t j)
{
    m_thermo->setTemperature(T(x,j));
    const doublereal* yy = x + m_nv*j + c_offset_Y;
    m_thermo->setMassFractions_NoNorm(yy);
    m_thermo->setPressure(m_press);
}

void StFlow::setGasAtMidpoint(const doublereal* x, size_t j)
{
    m_thermo->setTemperature(0.5*(T(x,j)+T(x,j+1)));
    const doublereal* yyj = x + m_nv*j + c_offset_Y;
    const doublereal* yyjp = x + m_nv*(j+1) + c_offset_Y;
    for (size_t k = 0; k < m_nsp; k++) {
        m_ybar[k] = 0.5*(yyj[k] + yyjp[k]);
    }
    m_thermo->setMassFractions_NoNorm(m_ybar.data());
    m_thermo->setPressure(m_press);
}

void StFlow::advanceChemistry(double* xg, double dt)
{
    double* x = xg + loc();

    for (size_t j=0; j!=nPoints(); j++)
    {
        setGas(x, j);

        IdealGasConstPressureReactor combustor;
        combustor.setThermoMgr(phase());
        combustor.setKineticsMgr(kinetics());

        // set simulation
        ReactorNet sim;
        sim.addReactor(combustor);

        // integrate
        sim.advance(dt);

        // update solution vector
        setScalars(x, j, combustor.contents());
    }

}

void StFlow::setScalars(double* x, size_t j, const ThermoPhase& ph)
{
    T(x, j) = ph.temperature();
    for (size_t i=0; i!=m_nsp; i++)
    {
        Y(x, i, j) = ph.massFraction(i);
    }
}

// protected member functions

void StFlow::eval(size_t jg, doublereal* xg,
                  doublereal* rg, integer* diagg, doublereal rdt)
{
    // if evaluating a Jacobian, and the global point is outside the domain of
    // influence for this domain, then skip evaluating the residual
    if (jg != npos && (jg + 1 < firstPoint() || jg > lastPoint() + 1)) {
        return;
    }

    // start of local part of global arrays
    doublereal* x = xg + loc();
    doublereal* rsd = rg + loc();
    integer* diag = diagg + loc();

    size_t jmin, jmax;
    if (jg == npos) { // evaluate all points
        jmin = 0;
        jmax = m_points - 1;
    } else { // evaluate points for Jacobian
        size_t jpt = (jg == 0) ? 0 : jg - firstPoint();
        jmin = std::max<size_t>(jpt, 1) - 1;
        jmax = std::min(jpt+1,m_points-1);
    }

    updateProperties(jg, x, jmin, jmax);
    evalResidual(x, rsd, diag, rdt, jmin, jmax);
}

void StFlow::evalResidual(double* x, double* rsd, int* diag,
                          double rdt, size_t jmin, size_t jmax)
{
    //----------------------------------------------------
    // evaluate the residual equations at all required
    // grid points
    //----------------------------------------------------

    for (size_t j = jmin; j <= jmax; j++) 
    {
        if (j == 0) 
        {
            //----------------------------------------------
            //         left boundary
            //----------------------------------------------
            evalLeftBoundary(x, rsd, diag, rdt);
        } 
        else if (j == m_points - 1) 
        {
            //----------------------------------------------
            //         right boundary
            //----------------------------------------------
            evalRightBoundary(x, rsd, diag, rdt);
        } 
        else 
        {
            //----------------------------------------------
            //         interior points
            //----------------------------------------------
            evalContinuity(j, x, rsd, diag, rdt);

            evalRadialMomentum(j, x, rsd, diag, rdt);

            evalSpecies(j, x, rsd, diag, rdt);

            evalEnergy(j, x, rsd, diag, rdt);

            rsd[index(c_offset_L, j)] = lambda(x,j) - lambda(x,j-1);
            diag[index(c_offset_L, j)] = 0;
        }
        // set residual of poisson's equ to zero
        rsd[index(c_offset_E, j)] = x[index(c_offset_E, j)];
    }
}

void StFlow::evalLeftBoundary(double* x, double* rsd, int* diag, double rdt)
{
    // Continuity. This propagates information right-to-left, since
    // rho_u at point 0 is dependent on rho_u at point 1, but not on
    // mdot from the inlet.

    size_t m = coordinatesType();

    rsd[index(c_offset_U,0)]
    = 
    -
    (
        rho_u(x,1) * pow(z(1), m)
        -
        rho_u(x,0) * pow(z(0), m)
    ) / dz(0) / pow(z(0), m);

    if (domainType() == cAxisymmetricStagnationFlow) {
        rsd[index(c_offset_U,0)] 
        -= 
        (density(1)*V(x,1) + density(0)*V(x,0));
    } 

    // the inlet (or other) object connected to this one will modify
    // these equations by subtracting its values for V, T, and mdot. As
    // a result, these residual equations will force the solution
    // variables to the values for the boundary object

    rsd[index(c_offset_L,0)] = -rho_u(x,0);

    rsd[index(c_offset_V,0)] = V(x,0);

    if (doEnergy(0)) {
        rsd[index(c_offset_T,0)] = T(x,0);
    } else {
        rsd[index(c_offset_T,0)] = T(x,0) - T_fixed(0);
    }

    // The default boundary condition for species is zero flux. However,
    // the boundary object may modify this.
    double sum = 0.0;
    for (size_t k = 0; k < m_nsp; k++) {
        sum += Y(x,k,0);
        rsd[index(c_offset_Y+k,0)] = - m_flux(k,0) - rho_u(x,0)*Y(x,k,0);
    }
    rsd[index(c_offset_Y + leftExcessSpecies(),0)] = 1.0 - sum;
}

void StFlow::evalRightBoundary(double* x, double* rsd, int* diag, double rdt)
{
    size_t j = m_points - 1;

    rsd[index(c_offset_L,j)] = lambda(x,j) - lambda(x,j-1);

    // the boundary object connected to the right of this one may modify or
    // replace these equations. The default boundary conditions are zero u, V,
    // and T, and zero diffusive flux for all species.

    rsd[index(c_offset_U,j)] = rho_u(x,j);

    // continuity for the stationary flame
    size_t m = coordinatesType();

    if (domainType() == cFreeFlow) {
        rsd[index(c_offset_U,j)] 
        = 
        -
        (
            rho_u(x,j) * pow(z(j), m)
            -
            rho_u(x,j-1) * pow(z(j-1), m)
        ) / dz(j-1) / pow(z(j), m);
    }

    rsd[index(c_offset_V,j)] = V(x,j);

    // Zhen Lu 210924
    if (doEnergy(j)) {
        rsd[index(c_offset_T,j)] = T(x,j);
    } else {
        rsd[index(c_offset_T, j)] = T(x,j) - T_fixed(j);
    }

    doublereal sum = 0.0;
    for (size_t k = 0; k < m_nsp; k++) {
        sum += Y(x,k,j);
        rsd[index(c_offset_Y+k,j)] = m_flux(k,j-1) + rho_u(x,j)*Y(x,k,j);
    }
    rsd[index(c_offset_Y + rightExcessSpecies(),j)] = 1.0 - sum;
}

void StFlow::evalContinuity(size_t j, double* x, double* rsd, int* diag, double rdt)
{
    //----------------------------------------------
    //    Continuity equation
    //----------------------------------------------

    size_t m = coordinatesType();

    //algebraic constraint
    diag[index(c_offset_U, j)] = 0;

    if (domainType() == cAxisymmetricStagnationFlow) {
        //----------------------------------------------
        //    d(\rho u)/dz + 2\rho V = 0
        //----------------------------------------------
        // Note that this propagates the mass flow rate information to the left
        // (j+1 -> j) from the value specified at the right boundary. The
        // lambda information propagates in the opposite direction.
        rsd[index(c_offset_U,j)] =
            -(rho_u(x,j+1) - rho_u(x,j))/dz(j)
            -(density(j+1)*V(x,j+1) + density(j)*V(x,j));
    } else if (domainType() == cFreeFlow) {
        if (z(j) == m_zfixed) {
            if (m_do_energy[j]) {
                rsd[index(c_offset_U,j)] = T(x,j) - m_tfixed;
            } else {
                rsd[index(c_offset_U,j)] = rho_u(x,j) - m_rho[0]*0.3;
            }
        } else {
            if (z(j) > m_zfixed) {
                rsd[index(c_offset_U,j)] 
                = 
                -
                (
                    rho_u(x,j) * pow(z(j), m)
                    -
                    rho_u(x,j-1) * pow(z(j-1), m)
                ) / dz(j-1) / pow(z(j), m);
            } else if (z(j) < m_zfixed) {
                rsd[index(c_offset_U,j)] 
                = 
                -
                (
                    rho_u(x,j+1) * pow(z(j+1), m)
                    -
                    rho_u(x,j) * pow(z(j), m)
                ) / dz(j) / pow(z(j), m);
            }
        } 
    } 
    else 
    {
        rsd[index(c_offset_U,j)] 
        = 
        -
        (
            rho_u(x,j+1) * pow(z(j+1), m)
            -
            rho_u(x,j) * pow(z(j), m)
        ) / dz(j) / pow(z(j), m)
        -
        rdt * (m_rho[j] - m_rho_last[j]);
    }
}

void StFlow::evalRadialMomentum(size_t j, double* x, double* rsd, int* diag, double rdt)
{
    //------------------------------------------------
    //    Radial momentum equation
    //
    //    \rho dV/dt + \rho u dV/dz + \rho V^2
    //       = d(\mu dV/dz)/dz - lambda
    //-------------------------------------------------
    rsd[index(c_offset_V,j)]
    = 
    (
        shear(x,j) 
        -
        lambda(x,j)
        -
        rho_u(x,j) * dVdz(x,j)
        - m_rho[j] * V(x,j) * V(x,j)
    ) / m_rho[j]
    - rdt * (V(x,j) - V_prev(j));

    diag[index(c_offset_V, j)] = 1;
}

void StFlow::evalSpecies(size_t j, double* x, double* rsd, int* diag, double rdt)
{
    //-------------------------------------------------
    //    Species equations
    //
    //   \rho dY_k/dt + \rho u dY_k/dz + dJ_k/dz
    //   = M_k\omega_k
    //-------------------------------------------------

    for (size_t k = 0; k < m_nsp; k++) 
    {
        rsd[index(c_offset_Y + k, j)] = - rho_u(x,j) * dYdz(x,k,j)
                                        - divDiffFlux(k,j);
        if ( m_do_reaction )
        {
            rsd[index(c_offset_Y + k, j)] += m_wt[k] * wdot(k,j);
        }
        rsd[index(c_offset_Y + k, j)] /= m_rho[j];
        rsd[index(c_offset_Y + k, j)] -= rdt * (Y(x,k,j) - Y_prev(k,j));

        diag[index(c_offset_Y + k, j)] = 1;
    }
}

void StFlow::evalEnergy(size_t j, double* x, double* rsd, int* diag, double rdt)
{
    //-----------------------------------------------
    //    energy equation
    //
    //    \rho c_p dT/dt + \rho c_p u dT/dz
    //    = d(k dT/dz)/dz
    //      - sum_k(\omega_k h_k_ref)
    //      - sum_k(J_k c_p_k / M_k) dT/dz
    //-----------------------------------------------
    if (m_do_energy[j]) {
        setGas(x,j);

        // heat release term
        const vector_fp& h_RT = m_thermo->enthalpy_RT_ref();
        const vector_fp& cp_R = m_thermo->cp_R_ref();

        double hrr = 0.0;
        double sum_flux = 0.0;

        for (size_t k = 0; k < m_nsp; k++) {
            hrr += wdot(k,j) * h_RT[k];

            double flxk = (
                m_flux(k,j) * dz(j-1)
                +
                m_flux(k,j-1) * dz(j)
            ) / d2z(j);
            sum_flux += flxk * cp_R[k] / m_wt[k];
        }

        hrr *= GasConstant * T(x,j);

        double dtdzj = dTdz(x,j);
        sum_flux *= GasConstant * dtdzj;

        // convection and diffusion
        rsd[index(c_offset_T, j)] = - m_cp[j]*rho_u(x,j)*dtdzj
                                    - divHeatFlux(x,j) - sum_flux;
        // heat release
        if ( m_do_reaction )
        {
            rsd[index(c_offset_T, j)] -= hrr;
        }
        // Zhen Lu 211027 ignition
        if ( m_do_ignition ) 
        {
            rsd[index(c_offset_T, j)] += ignEnergy(j);
        }
        rsd[index(c_offset_T, j)] /= (m_rho[j]*m_cp[j]);
        rsd[index(c_offset_T, j)] -= rdt*(T(x,j) - T_prev(j));
        rsd[index(c_offset_T, j)] -= (m_qdotRadiation[j] / (m_rho[j] * m_cp[j]));
        diag[index(c_offset_T, j)] = 1;
    } else {
        // residual equations if the energy equation is disabled
        rsd[index(c_offset_T, j)] = T(x,j) - T_fixed(j);
        diag[index(c_offset_T, j)] = 0;
    }
}

void StFlow::evalContinuityResidualJacobian(double* xg, double* rg, 
                                            double* dlg, double* dg, double* dug,
                                            double rdt)
{
    size_t jmin = 0;
    size_t jmax = nPoints()-1;

    double* x = xg + loc();

    // update density
    updateThermo(x, jmin, jmax);

    size_t iloc = locVelocity();
    double* r = rg + iloc;
    double* d = dg + iloc;
    double* dl = dlg + iloc;
    double* du = dug + iloc;

    // coordinates type
    size_t m = coordinatesType();

    // left boundary
    // density
    setGasAtMidpoint(x, jmin);
    double rho = m_thermo->density();
    setGasAtMidpoint(m_slast.data(), jmin);
    double rhoPrev = m_thermo->density();

    r[jmin] = -(rho_u(x,jmin+1)*pow(z(jmin+1)/zm(jmin),m)
                -rho_u(x,jmin)*pow(z(jmin)/zm(jmin),m))
              - dz(jmin) * rdt * (rho-rhoPrev);
    // diagonals
    d[jmin] = -density(jmin) * pow(z(jmin)/zm(jmin), m);
    du[jmin] = density(jmin+1) * pow(z(jmin+1)/zm(jmin), m);
    //writelog("\n {:4d} {:10.4g} {:10.4g} {:10.4g}", 
    //         jmin, rg[iloc+jmin], dg[iloc+jmin], dug[iloc+jmin]);

    // interior points
    for (size_t j = jmin+1; j != nPoints(); j++) 
    {
        // residual
        if ( divScheme() == 1)
        {
            // central difference
            // density
            setGasAtMidpoint(x, j-1);
            double rho = m_thermo->density();
            setGasAtMidpoint(m_slast.data(), j-1);
            double rhoPrev = m_thermo->density();

            r[j] = -(rho_u(x,j)*pow(z(j)/zm(j-1),m)-rho_u(x,j-1)*pow(z(j-1)/zm(j-1),m))
                   -dz(j-1)*rdt*(rho-rhoPrev);

            // diagonals
            d[j] = density(j) * pow(z(j)/zm(j-1), m);
            dl[j-1] = - density(j-1) * pow(z(j-1)/zm(j-1), m);
        }
        else
        {
            // upwind
            r[j] = -(rho_u(x,j)-rho_u(x,j-1) * pow(z(j-1)/z(j),m))
                   -dz(j-1)*rdt*(m_rho[j]-m_rho_last[j]);

            // diagonals
            d[j] = density(j);
            dl[j-1] = - density(j-1) * pow(z(j-1)/z(j), m);
        }
        if (j != jmin+1)
            du[j-1] = 0.0;
        //writelog("\n {:4d} {:10.4g} {:10.4g} {:10.4g} {:10.4g}", 
        //         j, rg[iloc+j], dlg[iloc+j-1], dg[iloc+j], dug[iloc+j]);
    }
}

void StFlow::evalScalarLeftBoundary(double* x, double* rsd, double dt)
{
    // the inlet (or other) object connected to this one will modify
    // these equations by subtracting its values for V, T, and mdot. As
    // a result, these residual equations will force the solution
    // variables to the values for the boundary object

    if (doEnergy(0)) {
        rsd[indexScalar(cOffsetScalarT,0)] = -T(x,0);
    } else {
        rsd[indexScalar(cOffsetScalarT,0)] = -T(x,0)+T_fixed(0);
    }

    // The default boundary condition for species is zero flux. However,
    // the boundary object may modify this.
    double sum = 0.0;
    for (size_t k = 0; k < m_nsp; k++) {
        sum += Y(x,k,0);
        rsd[indexScalar(cOffsetScalarY+k,0)] = - m_flux(k,0) - rho_u(x,0)*Y(x,k,0);
    }
    rsd[indexScalar(cOffsetScalarY+leftExcessSpecies(),0)] = 1.0 - sum;
}

void StFlow::evalScalarRightBoundary(double* x, double* rsd, double dt)
{
    size_t j = nPoints() - 1;

    if (doEnergy(j)) {
        rsd[indexScalar(cOffsetScalarT, j)] = -T(x,j);
    } else {
        rsd[indexScalar(cOffsetScalarT, j)] = -T(x,j)+T_fixed(j);
    }

    doublereal sum = 0.0;
    for (size_t k = 0; k < m_nsp; k++) {
        sum += Y(x,k,j);
        rsd[indexScalar(cOffsetScalarY+k,j)] = m_flux(k,j-1) + rho_u(x,j)*Y(x,k,j);
    }
    rsd[indexScalar(cOffsetScalarY+rightExcessSpecies(),j)] = 1.0 - sum;
}

void StFlow::evalScalarSpecies(size_t j, double *x, double* r, double dt)
{
    double* Yc = x + index(c_offset_Y, j);
    size_t kExcess = distance(Yc, max_element(Yc, Yc + m_nsp));
    double sum = 0.0;

    for (size_t k = 0; k < m_nsp; k++) 
    {
        r[indexScalar(cOffsetScalarY+k, j)] = - rho_u(x,j) * dYdz(x,k,j)
                                              - divDiffFlux(k,j);
        if ( m_do_reaction )
        {
            r[indexScalar(cOffsetScalarY+k, j)] += m_wt[k] * wdot(k,j);
        }
        r[indexScalar(cOffsetScalarY+k, j)] /= m_rho[j];
        r[indexScalar(cOffsetScalarY+k, j)] *= dt;
        r[indexScalar(cOffsetScalarY+k, j)] -= (Y(x,k,j) - Y_prev(k,j));

        sum += Y(x,k,j);
    }
    r[indexScalar(cOffsetScalarY+kExcess, j)] = 1.0 - sum;
}

void StFlow::evalScalarTemperature(size_t j, double *x, double* r, double dt)
{
    if (m_do_energy[j])
    {
        setGas(x,j);

        // heat release term
        const vector_fp& h_RT = m_thermo->enthalpy_RT_ref();
        const vector_fp& cp_R = m_thermo->cp_R_ref();

        double hrr = 0.0;
        double sum_flux = 0.0;

        for (size_t k = 0; k < m_nsp; k++) {
            hrr += wdot(k,j) * h_RT[k];

            double flxk = (
                m_flux(k,j) * dz(j-1)
                +
                m_flux(k,j-1) * dz(j)
            ) / d2z(j);
            sum_flux += flxk * cp_R[k] / m_wt[k];
        }

        hrr *= GasConstant * T(x,j);

        double dtdzj = dTdz(x,j);
        sum_flux *= GasConstant * dtdzj;

        // convection and diffusion
        r[indexScalar(cOffsetScalarT, j)] = - m_cp[j]*rho_u(x,j)*dtdzj;
        r[indexScalar(cOffsetScalarT, j)] -= divHeatFlux(x, j);
        r[indexScalar(cOffsetScalarT, j)] -= sum_flux;
        // heat release
        if ( m_do_reaction )
            r[indexScalar(cOffsetScalarT, j)] -= hrr;
        // ignition
        if ( m_do_ignition ) 
            r[indexScalar(cOffsetScalarT, j)] += ignEnergy(j);
        r[indexScalar(cOffsetScalarT, j)] -= m_qdotRadiation[j];
        r[indexScalar(cOffsetScalarT, j)] /= (m_rho[j]*m_cp[j]);
        r[indexScalar(cOffsetScalarT, j)] *= dt;
        r[indexScalar(cOffsetScalarT, j)] -= (T(x,j) - T_prev(j));
    }
    else
    {
        // residual equations if the energy equation is disabled
        r[indexScalar(cOffsetScalarT, j)] = -T(x,j) + T_fixed(j);
    }
}

void StFlow::evalScalarResidual(double* x, double* rsd,
                                double dt, size_t jmin, size_t jmax)
{
    //----------------------------------------------------
    // evaluate the residual equations at all required
    // grid points
    //----------------------------------------------------

    for (size_t j = jmin; j <= jmax; j++) 
    {
        if (j == 0) 
        {
            //----------------------------------------------
            //         left boundary
            //----------------------------------------------
            evalScalarLeftBoundary(x, rsd, dt);
        } 
        else if (j == nPoints()-1) 
        {
            //----------------------------------------------
            //         right boundary
            //----------------------------------------------
            evalScalarRightBoundary(x, rsd, dt);
        } 
        else 
        {
            //----------------------------------------------
            //         interior points
            //----------------------------------------------
            evalScalarSpecies(j, x, rsd, dt);

            evalScalarTemperature(j, x, rsd, dt);
        }
    }
}

void StFlow::evalScalar(size_t jg, double* xg, double* rg, double dt)
{
    // if evaluating a Jacobian, and the global point is outside the domain of
    // influence for this domain, then skip evaluating the residual
    if (jg != npos && (jg + 1 < firstPoint() || jg > lastPoint() + 1)) {
        return;
    }

    // start of local part of global arrays
    double* x = xg + loc();
    double* rsd = rg + locScalar();

    size_t jmin, jmax;
    if (jg == npos) { // evaluate all points
        jmin = 0;
        jmax = m_points - 1;
    } else { // evaluate points for Jacobian
        size_t jpt = (jg == 0) ? 0 : jg - firstPoint();
        jmin = std::max<size_t>(jpt, 1) - 1;
        jmax = std::min(jpt+1,m_points-1);
    }

    updateProperties(jg, x, jmin, jmax);
    evalScalarResidual(x, rsd, dt, jmin, jmax);
}

double StFlow::evalMaxCFL(double* xg, double dt)
{
    double maxCFL = SmallNumber;
    double* x = xg + loc();

    size_t jmin = 0;
    size_t jmax = nPoints() - 1;

    updateThermo(x, jmin, jmax);
    updateTransport(x, jmin, jmax);
    updateDiffFluxes(x, jmin, jmax);

    for ( size_t j = jmin; j != jmax; j++)
    {
        // convection
        double convCFL = dt*abs(u(x,j))/dz(j);
        maxCFL = std::max(maxCFL, convCFL);
        // diffusion of species
        for ( size_t k = 0; k != m_nsp; k++)
        {
            // approximation
            setGasAtMidpoint(x,j);
            double rho = m_thermo->density();
            double diffCFL = dt 
                            *abs(m_flux(k,j))
                            /(rho*0.5*(Y(x,k,j)+Y(x,k,j+1))*dz(j));
            maxCFL = std::max(maxCFL, diffCFL);
        }
    }

    return maxCFL;
}

void StFlow::updateProperties(size_t jg, double* x, size_t jmin, size_t jmax)
{
    // properties are computed for grid points from j0 to j1
    size_t j0 = std::max<size_t>(jmin, 1) - 1;
    size_t j1 = std::min(jmax+1,m_points-1);

    updateThermo(x, j0, j1);
    if (jg == npos || m_force_full_update) {
        // update transport properties only if a Jacobian is not being
        // evaluated, or if specifically requested
        updateTransport(x, j0, j1);
    }
    if (jg == npos) {
        double* Yleft = x + index(c_offset_Y, jmin);
        m_kExcessLeft = distance(Yleft, max_element(Yleft, Yleft + m_nsp));
        double* Yright = x + index(c_offset_Y, jmax);
        m_kExcessRight = distance(Yright, max_element(Yright, Yright + m_nsp));
    }

    // update the species diffusive mass fluxes whether or not a
    // Jacobian is being evaluated
    updateDiffFluxes(x, j0, j1);

    updateRadiation(x, j0, j1);

    updateReaction(x, j0, j1);
}

void StFlow::updateThermo(const doublereal* x, size_t j0, size_t j1) {
    for (size_t j = j0; j <= j1; j++) {
        setGas(x,j);
        m_rho[j] = m_thermo->density();
        m_wtm[j] = m_thermo->meanMolecularWeight();
        m_cp[j] = m_thermo->cp_mass();
    }
}

void StFlow::updateTransport(doublereal* x, size_t j0, size_t j1)
{
     if (m_do_multicomponent) {
        for (size_t j = j0; j < j1; j++) {
            setGasAtMidpoint(x,j);
            doublereal wtm = m_thermo->meanMolecularWeight();
            doublereal rho = m_thermo->density();
            m_visc[j] = (m_dovisc ? m_trans->viscosity() : 0.0);
            m_trans->getMultiDiffCoeffs(m_nsp, &m_multidiff[mindex(0,0,j)]);

            // Use m_diff as storage for the factor outside the summation
            for (size_t k = 0; k < m_nsp; k++) {
                m_diff[k+j*m_nsp] = m_wt[k] * rho / (wtm*wtm);
            }

            m_tcon[j] = m_trans->thermalConductivity();
            if (m_do_soret) {
                m_trans->getThermalDiffCoeffs(m_dthermal.ptrColumn(0) + j*m_nsp);
            }
        }
    } else { // mixture averaged transport
        for (size_t j = j0; j < j1; j++) {
            setGasAtMidpoint(x,j);
            m_visc[j] = (m_dovisc ? m_trans->viscosity() : 0.0);
            m_trans->getMixDiffCoeffs(&m_diff[j*m_nsp]);
            m_tcon[j] = m_trans->thermalConductivity();
        }
    }
}

void StFlow::updateDiffFluxes(const doublereal* x, size_t j0, size_t j1)
{
    if (m_do_multicomponent) {
        for (size_t j = j0; j < j1; j++) {
            double dz = z(j+1) - z(j);
            for (size_t k = 0; k < m_nsp; k++) {
                doublereal sum = 0.0;
                for (size_t m = 0; m < m_nsp; m++) {
                    sum += m_wt[m] * m_multidiff[mindex(k,m,j)] * (X(x,m,j+1)-X(x,m,j));
                }
                m_flux(k,j) = sum * m_diff[k+j*m_nsp] / dz;
            }
        }
    } else {
        for (size_t j = j0; j < j1; j++) {
            double sum = 0.0;
            // Zhen Lu 210928
            // use midpoint properties
            setGasAtMidpoint(x,j);
            double rho = m_thermo->density();
            double wtm = m_thermo->meanMolecularWeight();
            for (size_t k = 0; k < m_nsp; k++) {
                m_flux(k,j) = m_wt[k] * (rho * m_diff[k+m_nsp*j] / wtm);
                m_flux(k,j) *= (X(x,k,j) - X(x,k,j+1))/dz(j);
                sum -= m_flux(k,j);
            }
            // correction flux to insure that \sum_k Y_k V_k = 0.
            for (size_t k = 0; k < m_nsp; k++) {
                m_flux(k,j) += sum * 0.5 * (Y(x,k,j) + Y(x,k,j+1));
            }
        }
    }

    if (m_do_soret) {
        for (size_t m = j0; m < j1; m++) {
            double gradlogT = 2.0 * (T(x,m+1) - T(x,m)) /
                              ((T(x,m+1) + T(x,m)) * (z(m+1) - z(m)));
            for (size_t k = 0; k < m_nsp; k++) {
                m_flux(k,m) -= m_dthermal(k,m)*gradlogT;
            }
        }
    }
}

void StFlow::updateRadiation(const double* x, size_t jmin, size_t jmax)
{
    // calculation of qdotRadiation

    // The simple radiation model used was established by Y. Liu and B. Rogg [Y.
    // Liu and B. Rogg, Modelling of thermally radiating diffusion flames with
    // detailed chemistry and transport, EUROTHERM Seminars, 17:114-127, 1991].
    // This model uses the optically thin limit and the gray-gas approximation
    // to simply calculate a volume specified heat flux out of the Planck
    // absorption coefficients, the boundary emissivities and the temperature.
    // The model considers only CO2 and H2O as radiating species. Polynomial
    // lines calculate the species Planck coefficients for H2O and CO2. The data
    // for the lines is taken from the RADCAL program [Grosshandler, W. L.,
    // RADCAL: A Narrow-Band Model for Radiation Calculations in a Combustion
    // Environment, NIST technical note 1402, 1993]. The coefficients for the
    // polynomials are taken from [http://www.sandia.gov/TNF/radiation.html].

    if (m_do_radiation) {
        // variable definitions for the Planck absorption coefficient and the
        // radiation calculation:
        doublereal k_P_ref = 1.0*OneAtm;

        // polynomial coefficients:
        const doublereal c_H2O[6] = {-0.23093, -1.12390, 9.41530, -2.99880,
                                     0.51382, -1.86840e-5};
        const doublereal c_CO2[6] = {18.741, -121.310, 273.500, -194.050,
                                     56.310, -5.8169};

        // calculation of the two boundary values
        double boundary_Rad_left = m_epsilon_left * StefanBoltz * pow(T(x, 0), 4);
        double boundary_Rad_right = m_epsilon_right * StefanBoltz * pow(T(x, m_points - 1), 4);

        // loop over all grid points
        for (size_t j = jmin; j < jmax; j++) {
            // helping variable for the calculation
            double radiative_heat_loss = 0;

            // calculation of the mean Planck absorption coefficient
            double k_P = 0;
            // absorption coefficient for H2O
            if (m_kRadiating[1] != npos) {
                double k_P_H2O = 0;
                for (size_t n = 0; n <= 5; n++) {
                    k_P_H2O += c_H2O[n] * pow(1000 / T(x, j), (double) n);
                }
                k_P_H2O /= k_P_ref;
                k_P += m_press * X(x, m_kRadiating[1], j) * k_P_H2O;
            }
            // absorption coefficient for CO2
            if (m_kRadiating[0] != npos) {
                double k_P_CO2 = 0;
                for (size_t n = 0; n <= 5; n++) {
                    k_P_CO2 += c_CO2[n] * pow(1000 / T(x, j), (double) n);
                }
                k_P_CO2 /= k_P_ref;
                k_P += m_press * X(x, m_kRadiating[0], j) * k_P_CO2;
            }

            // calculation of the radiative heat loss term
            radiative_heat_loss = 2 * k_P *(2 * StefanBoltz * pow(T(x, j), 4)
            - boundary_Rad_left - boundary_Rad_right);

            // set the radiative heat loss vector
            m_qdotRadiation[j] = radiative_heat_loss;
        }
    }
}

void StFlow::updateReaction(const double* x, size_t j0, size_t j1)
{
    for (size_t j = j0; j != j1; j++)
    {
        setGas(x, j);
        m_kin->getNetProductionRates(&m_wdot(0,j));
    }
}

doublereal StFlow::shear(const doublereal* x, size_t j) const 
{
    doublereal stress_l = m_visc[j-1]*(V(x,j) - V(x,j-1))/dz(j-1);
    doublereal stress_r = m_visc[j]*(V(x,j+1) - V(x,j))/dz(j);

    return 2.0 * (stress_r - stress_l) / d2z(j);
}

doublereal StFlow::divDiffFlux(size_t k, size_t j) const
{
    int m = coordinatesType();
    double div = ( m_flux(k,j)*pow(zm(j)/z(j),m)
                  -m_flux(k,j-1)*pow(zm(j-1)/z(j),m) )/d2z(j)*2.0;
    return div;
}

doublereal StFlow::divHeatFlux(const doublereal* x, size_t j) const
{
    double flux_l = m_tcon[j-1] * ( T(x,j) - T(x,j-1) ) / dz(j-1);
    double flux_r = m_tcon[j] * ( T(x,j+1) - T(x,j) ) / dz(j);
    int m = coordinatesType();
    double div = ( flux_r*pow(zm(j)/z(j),m)
                  -flux_l*pow(zm(j-1)/z(j),m) )/d2z(j)*2.0;
    return - div;
}

doublereal StFlow::dVdz(const doublereal* x, size_t j) const 
{
    size_t jloc = (u(x,j) > 0.0 ? j : j + 1);
    return (V(x,jloc) - V(x,jloc-1))/dz(jloc-1);
}

doublereal StFlow::dYdz(const doublereal* x, size_t k, size_t j) const 
{
    size_t jloc = (u(x,j) > 0.0 ? j : j + 1);
    return (Y(x,k,jloc) - Y(x,k,jloc-1))/dz(jloc-1);
}

doublereal StFlow::dTdz(const doublereal* x, size_t j) const 
{
    size_t jloc = (u(x,j) > 0.0 ? j : j + 1);
    return (T(x,jloc) - T(x,jloc-1))/dz(jloc-1);
}

double StFlow::scalarGradient(const vector_fp& s, const double v, size_t j) const
{
    switch (m_convectiveScheme)
    {
    case 2:
        return scalarGradientGamma(s, v, j);
    case 1:
        return scalarGradientLinear(s, v, j);
    default:
        return scalarGradientUpwind(s, v, j);
    }
}

double StFlow::scalarGradientUpwind(const vector_fp& s, const double v, size_t j) const
{
    int k = (v > 0.0 ? 0 : 1);
    return (s[k+1] - s[k])/dz(j+k-1);
}

double StFlow::scalarGradientLinear(const vector_fp& s, const double v, size_t j) const
{
    return (s[2] - s[0])/d2z(j);
}

double StFlow::scalarGradientGamma(const vector_fp& s, const double v, size_t j) const
{
    double grad_UD = scalarGradientUpwind(s, v, j);
    double grad_CD = scalarGradientLinear(s, v, j);

    double phi = grad_UD / (2*grad_CD);

    if (phi <= 0.0 || phi >= 1.0) 
    {
        return grad_UD;
    } 
    else if ( phi <= m_gammaSchemeBeta ) 
    {
        return phi/m_gammaSchemeBeta*grad_CD + (1-phi/m_gammaSchemeBeta)*grad_UD;
    } 
    else 
    {
        return grad_CD;
    }
}

double StFlow::ignEnergy(size_t j) const
{

    if (m_time >= m_ign_time) return 0.0;

    if (z(j) >= 4.0*m_ign_radius)
    {
        return 0.0;
    }
    else
    {
        double q = m_ign_energy 
                  *exp(-pow(z(j)/m_ign_radius, 2.0))
                  /(pow(Pi, 1.5)*pow(m_ign_radius, 3.0)*m_ign_time);
        return q;
    }
}

} // namespace