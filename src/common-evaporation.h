/**
# Common Evaporation

We define tolerances and functions useful
for evaporation models. */


#include "curvature.h"
#include "adapt_wavelet_leave_interface.h"

#ifndef F_ERR
# define F_ERR 1.e-10
#endif
#ifndef T_ERR
# define T_ERR 0.
#endif
#include "mass_refine_prolongation.h"

#ifndef I_TOL
# define I_TOL 1.e-8
#endif

/**
We define a custom for that simplifies
loops over the chemical species. */

#define foreach_elem(list, index) \
for (int index=0; index<list_len (list); index++)

/**
## Compute the normal in an interfacial cell
*/

coord normal (Point point, scalar c)
{
  coord n = mycs (point, c);
  double nn = 0.;
  foreach_dimension()
    nn += sq(n.x);
  nn = sqrt(nn);
  foreach_dimension()
    n.x /= nn;
  return n;
}

/**
## *mass2molefrac()*: Compute mole fractions from mass fractions:

* *X*: vector filled with mole fractions
* *W*: vector with the mass fractions
* *MW*: vector with the molecular weights of each species
* *NS*: total number of species (vectors length)
*/

void mass2molefrac (double * X, const double * W, const double * MW, const int NS)
{
  double MWmix = 0.;
  for (int i=0; i<NS; i++) {
    MWmix += W[i]/MW[i];
  }
  for (int i=0; i<NS; i++) {
    X[i] = W[i]/MW[i]/(MWmix + 1.e-10);
  }
}

/**
## *mole2massfrac()*: Compute mass fractions from mole fractions:

* *W*: vector filled with mole fractions
* *X*: vector with the mass fractions
* *MW*: vector with the molecular weights of each species
* *NS*: total number of species (vectors length)
*/

void mole2massfrac (double * W, const double * X, const double * MW, const int NS)
{
  double MWmix = 0.;
  for (int i=0; i<NS; i++) {
    MWmix += X[i]*MW[i];
  }
  for (int i=0; i<NS; i++) {
    W[i] = X[i]*MW[i]/(MWmix + 1.e-10);
  }
}

/**
## *avg_neighbor()*: Compute the average value of a scalar field Y in a 3x3 stencil around the current cell. */

double avg_neighbor (Point point, scalar Y, scalar f) {
  double fYnei = 0., fnei = 0.;
  foreach_neighbor(1) {
    fYnei += f[]*Y[];
    fnei  += f[];
  }
  return fYnei / fnei;
}

/**
## *avg_interface()*: Compute the average value of a scalar field Y in a region around the interface. */

double avg_interface (scalar Y, scalar f) {
  double Yavg = 0.;
  int counter = 0;
  foreach(reduction(+:Yavg) reduction(+:counter)) {
    if (f[] > F_ERR && f[] < 1.-F_ERR) {
      counter++;
      Yavg += Y[];
    }
  }
  return Yavg / counter;
}

