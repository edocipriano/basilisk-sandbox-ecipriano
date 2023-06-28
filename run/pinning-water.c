/**
# Pinning a 1mm Water Droplet

We want to extend the simulation in [pinning.c](pinning.c)
trying to pin a 1mm water droplet to a specific
point of the domain.

What we see is that, when a stronger density and viscosity
ratio is imposed, the droplet is deformed in a strange manner.
Trying to figure out what's wrong I noticed the following
things:

1. If the boundary condition for *h* is not imposed
the problem persists;
2. If we comment the line `f.height = h`, the problem
disappears; of course in that case the droplet is not
pinned, but it does not show strange deformations. The
difference, in this case, is that [contact.h](/src/contact.h)
does not update height.

From the animation, it can be observed that the droplet
remains fixed at the pinning point but the interface is
aggressively deformed.

![Evolution of the volume fraction field](pinning-water/movie.mp4)
*/

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "pinning.h"
#include "tension.h"
#include "reduced.h"
#include "view.h"

double scale = 1.e-3;

int main()
{
  L0 = 2.*scale;
  X0 = -0.5*L0;

  /**
  We use air-water properties. */

  rho1 = 1000., rho2 = 1.;
  mu1 = 1.e-3, mu2 = 1.e-5;

  /**
  We set the surface tension of water and the normal gravity value. */

  f.sigma = 0.073;
  G.x = -9.81;
  DT = 1.e-6;

  /**
  We set the coordinate of the pinning point. */

  pinning.ap = 0.75*scale;

  run();
}

/**
We initialize half liquid droplet on the bottom boundary. */

event init (t = 0)
{
  fraction (f, - (sq(x - 0.25*scale) + sq(y) - sq(0.5*scale)));
}

/**
At equilibrium (t = 10 seems sufficient), we output the interface
shape. */

event end (t = 0.01)
{
  char name[80];
  sprintf (name, "out-%.1f", fabs(G.x));

  FILE * fp = fopen (name, "w");
  output_facets (f, fp);
  fclose (fp);
}

/**
We write the animation with the evolution of the volume
fraction field and the gas-liquid interface. */

event movie (t += 0.0001)
{
  clear();
  view (ty = -0.5);
  box();
  draw_vof ("f", lw = 1.5);
  squares ("f", min = 0., max = 1.);
  save ("movie.mp4");
}

