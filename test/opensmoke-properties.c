
#include "run.h"
#include "utils.h"
#include "opensmoke-properties.h"

int main (void) {
  kinfolder = "skeletal/methanol";
  run();
}

/**
## Properties

We update the material properties according to the
thermodynamic state of the mixture.
*/

event properties (i++) {

  /**
  We first set the thermodynamic state of the mixture: */

  double x[OpenSMOKE_NumberOfSpecies()];
  double q[OpenSMOKE_NumberOfLiquidSpecies()];
  for (int jj=0; jj<OpenSMOKE_NumberOfSpecies(); jj++)
    x[jj] = 0.;
  x[OpenSMOKE_IndexOfSpecies("O2")] = 0.21;
  x[OpenSMOKE_IndexOfSpecies("N2")] = 0.79;
  q[0] = 1.;

  ts1.T = 300., ts1.P = 101325., ts1.x = q;
  ts2.T = 300., ts2.P = 101325., ts2.x = x;

  /**
  We update the properties of the variable fields. */

  fprintf (stderr, "Gas Properties\n");
  fprintf (stderr, "rho2 = %f\n", tp2.rho (&ts2));
  fprintf (stderr, "mu2 = %f\n", tp2.mu (&ts2));
  fprintf (stderr, "lambda2 = %f\n", tp2.lambda (&ts2));
  fprintf (stderr, "cp2 = %f\n", tp2.cp (&ts2));

  fprintf (stderr, "\n");

  fprintf (stderr, "Liquid Properties\n");
  fprintf (stderr, "rho1 = %f\n", tp1.rho (&ts1));
  fprintf (stderr, "mu1 = %f\n", tp1.mu (&ts1));
  fprintf (stderr, "lambda1 = %f\n", tp1.lambda (&ts1));
  fprintf (stderr, "cp1 = %f\n", tp1.cp (&ts1));

  for (unsigned int i=0; i<OpenSMOKE_NumberOfLiquidSpecies(); i++) {
    fprintf (stderr, "pvap[%s] = %f\n",
        OpenSMOKE_NamesOfLiquidSpecies (i), tp1.pvap (&ts1, i));
    fprintf (stderr, "dhev[%s] = %f\n",
        OpenSMOKE_NamesOfLiquidSpecies (i), tp1.dhev (&ts1, i));
  }
}

