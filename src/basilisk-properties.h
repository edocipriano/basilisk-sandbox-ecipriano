/**
# Basilisk Properties

We implement the calculation of thermodynamic properties, without
using external libraries.
*/

#include "variable-properties.h"

/**
## Properties Functions

Functions for the update of the density, given the thermodynamic
state.
*/

/**
### *opensmoke_gasprop_density()*: gas phase density according to the ideal gas low
*/

double gasprop_density (void * p) {
  return 1.;
}

/**
## Thermodynamic Properties

We create the instance of two structures with the
thermodynamic properties, *tp1* for the liquid phase
and *tp2* for the gas phase. The same nomenclature is used
for the thermodynamic states.
*/

ThermoProps tp1, tp2;
ThermoState ts1, ts2;

/**
## Initialization

We set the thermodynamic properties function pointers
to the specific opensmoke functions declared above.
*/

event defaults (i = 0) {

  /**
  We set the thermodynamic properties pointers to the
  correct opensmoke functions that compute material
  properties. */

  tp1.rhov    = NULL;
  tp1.muv     = NULL;
  tp1.lambdav = NULL;
  tp1.cpv     = NULL;
  tp1.pvap    = NULL;
  tp1.dhev    = NULL;
  tp1.diff    = NULL;

  tp2.rhov    = NULL;
  tp2.muv     = NULL;
  tp2.lambdav = NULL;
  tp2.cpv     = NULL;
  tp2.diff    = NULL;
}

