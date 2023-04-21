/**
# Fixed Flux Phase Change Model

This phase change model does not compute the vaporization
rate from a mass or energy balance, but it uses a constant
value *mEvapVal* defined by the user.
*/

#include "curvature.h"

extern double mEvapVal;

scalar mEvap[];
scalar * mEvapList = {mEvap};

event phasechange (i++)
{
  foreach() {
    mEvap[] = 0.;
    if (interfacial(point, f))
      mEvap[] = mEvapVal;
  }
}

