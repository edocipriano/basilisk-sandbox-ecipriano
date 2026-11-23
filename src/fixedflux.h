/**
# Fixed Flux Phase Change Model

This phase change model does not compute the vaporization
rate from a mass or energy balance, instead, it uses a constant
value *mEvapVal* defined by the user.
*/

#include "phasechange.h"

double mEvapVal = 0.;

/**
## Fields

We define the field *mEvap* with the vaporization rate
per unit of surface in every interfacial cell, and we
decleare the *mEvapList* used by [evaporation.h](evaporation.h).
*/

scalar mEvap[];
scalar * mEvapList = {mEvap};

/**
## Vaporization Rate

We loop over the interfacial cells and we assign the
value of *mEvapVal* in the field *mEvap*.
*/

void phasechange_rate (void) {
  foreach() {
    mEvap[] = 0.;
    if (f[] > F_ERR && f[] < 1.-F_ERR)
      mEvap[] = mEvapVal;
  }
}

/**
## Order of integration

The phase change step precedes the transport of the interface and of the
tracers: phase change, VOF, advection, diffusion, Navier--Stokes (see
[phasechange.h](phasechange.h)). */

event vof (i++) {
  phasechange_reset();
  phasechange_properties();
  phasechange_chemistry();
  phasechange_rate();
  phasechange_expansion();
  phasechange_divergence();
  phasechange_advection();
}

event tracer_diffusion (i++) {
  phasechange_diffusion();
  phasechange_reset_sources();
}

