#include "phase.h"

int NS;
double rhov = 1., muv = 1., lambda = 1., cp = 1.;
double T0 = 300., Pref = 101325., Dmix = 1.;

Phase * phase;

event defaults (i = 0) {
  phase = new_phase ("", NS, true);
}
