/**
# Pinning a Liquid Droplet - AXI

We want to test the module [pinning.h](../src/pinning.h), to suspend
a liquid droplet at a specific point of the domain. The aim is
to compute the equilibrium contact angle and to compare it with
the analytical solution:

$$
  \theta = \arccos \left( \dfrac{4\rho_l g r_d^3}{3 \sigma d_f}\right)
$$

given by the balance between gravity and surface tension force.

*/

#include "grid/multigrid.h"
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "pinning.h"
#include "tension.h"
#include "reduced.h"
#include "view.h"

/**
We set the diameter of the solid fiber which suspends
the liquid droplet. */

double df = 0.3;

int main()
{
  /**
  We set the length of the domain, and we move the origin
  according to the fiber diameter along y (specific to AXI).
  */

  size (2);
  origin (-1, 0.5*df);

  /**
  We set a small density and viscosity ratio. */

  mu1 = 1.;
  mu2 = 0.1;
  rho2 = 0.1;

  /**
  We set the surface tension coefficient and run for two position of
  the pinning point. */

  f.sigma = 1;
  pinning.ap = 0.75;

  /**
  We run the simulation at different values of gravity .*/

  for (G.x = 0.; G.x >= -1.8; G.x -= 0.2)
    run();
}

/**
We initialize half liquid droplet on the bottom boundary. */

event init (t = 0)
{
  fraction (f, - (sq(x - 0.25) + sq(y) - sq(0.5)));
}

/**
At equilibrium (t = 10 seems sufficient), we compute the
numerical contact angle. */

event end (t = 10)
{
  Point point = locate (pinning.ap, 0.5*df);
  double ca = atan (1./(h.x[0,-1] - h.x[]));
  //double exact = acos (4.*pow(0.5,3)*rho1*(-G.x)/(3.*f.sigma*df));

  fprintf (stderr, "%g %g\n", fabs (G.x), fabs (ca*180./pi));
}

/**
## Results

We plot the comparison between the equilibrium contact angle
and the one obtained from the numerical simulation.

~~~gnuplot
reset
set xrange[0:2]

f(x) = acos(4.*0.5**3*1.*x/(3.*1*0.3))*180/pi

plot f(x) lw 2 t "Theoretical", \
     "log" u 1:2 w p pt 8 ps 1.5 t "Numerical"
~~~
*/
