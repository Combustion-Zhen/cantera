/**
 * @file OneDimConst.h
 *
 * Constants for one-dimensional simulations.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_ONEDIMCONST_H
#define CT_ONEDIMCONST_H

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

const int cOffsetScalarT = 0;
const int cOffsetScalarY = 1;

// domain types
const int cFlowType = 50;
const int cFreeFlow = 51;
const int cAxisymmetricStagnationFlow = 52;
//  Zhen Lu
const int cRadialFlow = 53;
const int cPolarFlow = 54;
const int cTubularFlow = 55;

const int cConnectorType = 100;
const int cSurfType = 102;
const int cInletType = 104;
const int cSymmType = 105;
const int cOutletType = 106;
const int cEmptyType = 107;
const int cOutletResType = 108;
const int cPorousType = 109;

// Zhen Lu 210916
// coordinate types
const int cCartesian = 0;
const int cCylindrical = 1;
const int cSpherical = 2;
// temporal types
const int cSteady = 0;
const int cTransient = 1;

}

#endif
