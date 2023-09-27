/**
# Fixed Flux Droplet Evaporation

Evaporation of a liquid droplet with a fixed vaporization
flowrate. In these conditions, the droplet consumption
is obtained from a mass balance on the liquid droplet and
the analytic solution can be easily derived:
$$
  \dfrac{dR}{dt} = \dfrac{\dot{m}}{\rho_l}
$$

We want to test and compare the droplet consumption with
the exact solution, considering the presence of the Stefan
flow. The simulation setup used here was adapted from
[Malan et al, 2021](#malan2021geometric).
The animation shows the consuption of the liquid droplet,
that maintains a perfectly spherical shape throughtout the
entire lifetime, as well as the radial velocity field
caused by the phase change.
![Evolution of the interface position and the velocity field](fixedflux/movie.mp4)(width="400" height="400")
*/

/**
## Phase Change Setup

We divide the mass balance by the density of the gas phase,
according with the non-dimensional treatment reported in
the [paper](#malan2021geometric). There is no need to shift
the volume expansion term because we do not transport the
temperature or mass fractions in this test case.
*/

#define VARPROP

/**
## Simulation Setup

We use the centered Navier-Stokes equations solver with
the evaporation source term in the projection step. The
double pressure velocity coupling approach is adopted to
obtain a divergence-free velocity which can be used for
the VOF advection equation. The evporation model is
combined with the fixed flux mechanism, that imposes
a constant vaporization flowrate.
*/

//#include "grid/multigrid.h"
#include "axi.h"
#include "navier-stokes/centered-evaporation.h"
#define ufext uf
#include "two-phase-varprop.h"
#include "navier-stokes/conserving-varprop.h"
#include "basilisk-properties.h"
#include "tension.h"
#include "icentripetal.h"
#include "evaporation.h"
#include "thermalnew.h"
#include "view.h"

/**
### Boundary Conditions

Outflow boundary conditions are imposed on every
domain boundary since the droplet is initially placed
at the center of the domain. Because the
centered-doubled method is used, the boundary conditions
must be imposed also for the *extended* fields.
*/

double lambda1, lambda2, cp1, cp2, dhev, TIntVal, TL0, TG0;

u.n[top] = neumann (0.);
u.t[top] = neumann (0.);
p[top] = dirichlet (0.);
q1.n[top] = neumann (0.);
q2.n[top] = neumann (0.);

//u.n[bottom] = neumann (0.);
//u.t[bottom] = neumann (0.);
//p[bottom] = dirichlet (0.);
uf.n[bottom] = 0.;
uf.t[bottom] = dirichlet(0.);

u.n[left] = neumann (0.);
u.t[left] = neumann (0.);
p[left] = dirichlet (0.);
q1.n[left] = neumann (0.);
q2.n[left] = neumann (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);
q1.n[right] = neumann (0.);
q2.n[right] = neumann (0.);

/**
### Model Data

We set the initial radius of the droplet, and
the value of the vaporization rate.
*/

int maxlevel = 6, minlevel = 5;
double R0 = 1.e-4;
double radius0 = 0.;
double mEvapVal = -0.;
bool is_diffusion = true;
double rhoR = 0.080470865772519;
double tend = 1.;

int main(void) {
  //kinfolder = "materials/n-heptane";
  L0 = 8.*R0; X0 = -0.5*L0;

  /**
  We set the density and viscosity values. */

  rho1 = 975.91, rho2 = rhoR*rho1;
  mu1 = 3.7e-4, mu2 = mu1*1e-2;
  lambda1 = 0.6705, lambda2 = 0.03153;
  cp1 = 4285., cp2 = 1063.; dhev = 0.;
  TG0 = 175., TL0 = 350., TIntVal = TL0;
  Pref = 5e+6;

  /**
  The surface tension is set to zero in order to assess
  the sphericity of the droplet without being influenced
  by the surface tension. */

  f.sigma = 0.;
  f.gradient = zero;

#ifdef CENTRIPETAL
  sfm.sigma = 0.001;
  sfm.eps = 0.;
#endif

  /**
  We perform the simulation at different levels of
  refinement. */

  TOLERANCE = 1.e-6;
  tend *= rho2*cp2/lambda2*sq(R0);

  for (maxlevel = 7; maxlevel <= 8; maxlevel++) {
    init_grid (1 << maxlevel);
    run();
  }
}

#define circle(x, y, R) (sq(R) - sq(x) - sq(y))

/**
We initialize the volume fraction field.
*/

event init (i = 0) {
  fraction (f, circle (x, y, R0));

#if AXI
  radius0 = pow(3./2.*statsf(f).sum, 1./3.);
#else
  radius0 = sqrt (statsf(f).sum/pi);
#endif

  foreach() {
    f[] = 1. - f[];
    TL[] = f[]*TL0;
    TG[] = (1. - f[])*TG0;
  }
}

//event stability (i++) {
//  dtmax = rhoR*cp2*sq(L0/pow (2,maxlevel))/lambda2/2.;
//}

#if TREE
event adapt (i++) {
  adapt_wavelet_leave_interface ({T,u.x,u.y}, {f},
      (double[]){1.e-3,1.e-3,1.e-3}, maxlevel, minlevel, 1);
}
#endif

/**
## Post-Processing

The following lines of code are for post-processing purposes.
*/

/**
### Output Files

We write the droplet volume from the numerical simulation
as well as the exact solution.
*/

event output (i++) {
  char name[80];
  sprintf (name, "OutputData-%d", maxlevel);

  double dropvol = 0.;
  foreach(reduction(+:dropvol))
    dropvol += (1. - f[])*dv();

#if AXI
  double radius = pow(3./2.*dropvol, 1./3.);
#else
  double radius = sqrt (dropvol/pi);
#endif

  double tau = sq(1.e-4)*(rhoR*rho1*cp2/lambda2);

  static FILE * fp = fopen (name, "w");
  fprintf (fp, "%g %g %g %g\n", t, t/tau, radius, radius/radius0);
  fflush (fp);
}

/**
### Movie

We write the animation with the volume fraction field
and the velocity vectors.
*/

event movie (t += 0.001; t <= tend) {
  clear();
  view (ty = -0.5);
  draw_vof ("f", lw = 1.5);
  squares ("T", linear = true, min = min (TL0, TG0),
      max = max (TL0, TG0), map = blue_white_red);
  save ("movie.mp4");
}

