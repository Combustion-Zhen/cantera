//! @file vcs_species_thermo.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef VCS_SPECIES_THERMO_H
#define VCS_SPECIES_THERMO_H

#include "cantera/equil/vcs_defs.h"

namespace Cantera
{

class vcs_VolPhase;

// Models for the species standard state Naught temperature dependence
#define VCS_SS0_NOTHANDLED -1
#define VCS_SS0_CONSTANT 0
//#define VCS_SS0_NASA_POLY 1
#define VCS_SS0_CONSTANT_CP 2

// Models for the species standard state extra pressure dependence
#define VCS_SSSTAR_NOTHANDLED -1
#define VCS_SSSTAR_CONSTANT 0
#define VCS_SSSTAR_IDEAL_GAS 1

/*!
 * Identifies the thermo model for the species. This structure is shared by
 * volumetric and surface species. However, each will have its own types of
 * thermodynamic models. These quantities all have appropriate units.
 */
class VCS_SPECIES_THERMO
{
    /*
     * All objects are public for ease of development
     */
public:
    VCS_SPECIES_THERMO() = default;

    //! Index of the phase that this species belongs to.
    size_t IndexPhase = 0;

    //! Index of this species in the current phase.
    size_t IndexSpeciesPhase = 0;

    //! Pointer to the owning phase object.
    vcs_VolPhase* OwningPhase = nullptr;

    //! Integer representing the models for the species standard state Naught
    //! temperature dependence. They are listed above and start with VCS_SS0_...
    int SS0_Model = VCS_SS0_CONSTANT;

    //! Internal storage of the last calculation of the reference naught Gibbs
    //! free energy at SS0_TSave. (always in units of Kelvin)
    double SS0_feSave = 0.0;

    //! Internal storage of the last temperature used in the calculation of the
    //! reference naught Gibbs free energy. units = kelvin
    double SS0_TSave = -90.0;

    //! Base temperature used in the VCS_SS0_CONSTANT_CP model
    double SS0_T0 = 273.15;

    //! Base enthalpy used in the VCS_SS0_CONSTANT_CP model
    double SS0_H0 = 0.0;

    //! Base entropy used in the VCS_SS0_CONSTANT_CP model
    double SS0_S0 = 0.0;

    //! Base heat capacity used in the VCS_SS0_CONSTANT_CP model
    double SS0_Cp0 = 0.0;

    //! Value of the pressure for the reference state.
    double SS0_Pref = OneAtm;

    //! Integer value representing the star state model.
    int SSStar_Model = VCS_SSSTAR_CONSTANT;

    //! Models for the standard state volume of each species
    int SSStar_Vol_Model = VCS_SSVOL_IDEALGAS;

    //! parameter that is used in the VCS_SSVOL_CONSTANT model.
    double SSStar_Vol0 = -1;
};

}

#endif
