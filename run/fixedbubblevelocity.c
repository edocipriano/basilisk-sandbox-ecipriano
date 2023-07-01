/**
# Fixed Flux Bubble Expansion

Expansion of a bubble with a fixed vaporization rate.
This test case was proposed by [Tanguy et al., 2014](#tanguy2014benchmarks)
and it is interesting because it can be used to test
the velocity profile with the Stefan flow, and how the
velocity changes when the volume expansion term in the
continuity equation is shifted and/or distributed.

We perform this test case using the *Velocity Potential*
approach, and we compare the results both with the
theoretical velocity profile and with the numerical results
obtained by [Tanguy et al., 2014](#tanguy2014benchmarks).

The animation shows the growth of the bubble in time.
![Evolution of the interface position and the volume fraction field](fixedbubblevelocity/movie.mp4)(width="500" height="500")
*/

/**
## Phase Change Setup

We divide the mass balance by the density of the gas phase,
and we avoid the stefan flow shifting procedure.
*/

#define NOSHIFTING
#define BYRHOGAS

/**
## Simulation Setup

We use the centered Navier-Stokes equations solver with
the evaporation source term in the projection step. The
velocity potential approach is used to obtain a
divergence-free velocity which can be used for
the VOF advection equation. The evporation model is
combined with the fixed flux mechanism, that imposes
a constant vaporization flowrate. */

#include "navier-stokes/centered-evaporation.h"
#include "navier-stokes/velocity-potential.h"
#include "two-phase.h"
#include "tension.h"
#include "evaporation.h"
#include "fixedflux.h"
#include "view.h"

/**
### Boundary Conditions

Outflow boundary conditions are imposed on every
domain boundary since the bubble is initially placed
at the center of the domain. We set the BCs also for
the velocity potential *ps*. */

u.n[top] = neumann (0.);
u.t[top] = neumann (0.);
p[top] = dirichlet (0.);
ps[top] = dirichlet (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);
ps[right] = dirichlet (0.);

u.n[left] = neumann (0.);
u.t[left] = neumann (0.);
p[left] = dirichlet (0.);
ps[left] = dirichlet (0.);

u.n[bottom] = neumann (0.);
u.t[bottom] = neumann (0.);
p[bottom] = dirichlet (0.);
ps[bottom] = dirichlet (0.);

/**
### Model Data

We set the initial radius of the bubble,
the value of the vaporization rate, and
the maximum level of refinement. */

int maxlevel = 7;
double D0 = 0.002, mEvapVal = -0.1;
double effective_radius0;

int main (void) {
  /**
  We set the density and viscosity values. */

  rho1 = 1000., rho2 = 1.;
  mu1 = 0.001, mu2 = 1.78e-5;

  /**
  We set the surface tension coefficient. */

  f.sigma = 0.07;

  /**
  We change the dimension of the domain and
  we shift the origin. */

  L0 = 0.008;
  origin (-0.5*L0, -0.5*L0);

  /**
  We setup the grid and run the simulation. */

  init_grid (1 << maxlevel);
  run();
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

/**
We initialize the volume fraction field, and from that
field we compute the numerical initial radius. */

event init (i = 0) {
  fraction (f, -circle(x,y,0.5*D0));

  double gasvol = 0.;
  foreach(reduction(+:gasvol))
    gasvol += (1. - f[])*dv();
  effective_radius0 = sqrt (gasvol/pi);
}

/**
We adapt the grid according to the velocity and the
volume fraction fields. */

#if TREE
event adapt (i++) {
  adapt_wavelet_leave_interface ({u.x,u.y}, {f},
      (double[]){1.e-3,1.e-3}, maxlevel, 1);
}
#endif

/**
## Post-Processing

The following lines of code are for post-processing purposes.
*/

/**
### Exact Solution

We write a function that computes the analytic
solution for the bubble radius.
*/

double exact (double t) {
  return effective_radius0 - mEvapVal/rho2*t;
}

/**
### Output Files

We write the bubble radius from the numerical simulation,
the exact radius and the relative error. */

event output_data (i++) {
  char name[80];
  sprintf (name, "OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, "w");

  double gasvol = 0.;
  foreach(reduction(+:gasvol))
    gasvol += (1. - f[])*dv();
  double effective_radius = sqrt (gasvol/pi);
  double relerr = fabs (exact(t) - effective_radius) / exact(t);
  fprintf (fp, "%g %g %g %g\n", t, effective_radius, exact(t), relerr);
  fflush (fp);
}

/**
### Velocity Profile

We write on a file the velocity profile at a specific
time step. */

event profiles (t = 0.005) {
  char name[80];
  sprintf (name, "Profiles-%d", maxlevel);
  FILE * fp = fopen (name, "w");

  /**
  We create an array with the mass fraction profiles
  for each processor. */

  Array * arruf = array_new();
  Array * arrus = array_new();
  for (double x = 0.; x < 0.5*L0; x += 0.05*L0/(1 << maxlevel)) {
    double val_uf = interpolate (uf.x,  x, 0.);
    double val_us = interpolate (ufs.x, x, 0.);
    val_uf = (val_uf == nodata) ? 0. : val_uf;
    val_us = (val_us == nodata) ? 0. : val_us;
    array_append (arruf, &val_uf, sizeof(double));
    array_append (arrus, &val_us, sizeof(double));
  }
  double * val_ufp = (double *)arruf->p;
  double * val_usp = (double *)arrus->p;

  /**
  We sum each element of the arrays in every processor. */

  @if _MPI
  int size = arruf->len/sizeof(double);
  MPI_Allreduce (MPI_IN_PLACE, val_ufp, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, val_usp, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  @endif

  /**
  The master node writes the profiles on a file. */

  if (pid() == 0) {
    int count = 0;
    for (double x = 0.; x < 0.5*L0; x += 0.05*L0/(1 << maxlevel)) {
      fprintf (fp, "%g %g %g\n", x, val_ufp[count], val_usp[count]);
      count++;
    }
    fprintf (fp, "\n\n");
    fflush (fp);
  }
  fclose (fp);
  array_free (arruf);
  array_free (arrus);
}

/**
### Movie

We write the animation with the evolution of the
volume fraction field and the gas-liquid interface. */

event movie (t += 0.0001; t <= 0.01) {
  clear();
  draw_vof ("f", lw = 1.5);
  squares ("f", min = 0., max = 1., linear = true);
  save ("movie.mp4");
}

/**
## Results

We compare the velocity profile with the theretical
solution and with the results obtained by [Tanguy et al., 2014](#tanguy2014benchmarks)
(Fig. 5 (b), (c)).

~~~gnuplot Radial velocity profile
reset
set xlabel "radius [m]"
set ylabel "Velocity [m/s]"
set key top right
set size square
set grid

plot "../data/tanguy-fixedbubblevelocity-theoretical.csv" u 1:2 w p pt 6 t "Theoretical", \
     "../data/tanguy-fixedbubblevelocity-b.csv" u 1:2 w p pt 6 t "Tanguy et al., 2014 b", \
     "../data/tanguy-fixedbubblevelocity-c.csv" u 1:2 w p pt 6 t "Tanguy et al., 2014 c", \
     "Profiles-7" u 1:2 w l lw 1.2 t "Field Velocity", \
     "Profiles-7" u 1:3 w l lw 1.2 t "Stefan Velocity"

~~~

~~~gnuplot Evolution of the radius in time
reset
set xlabel "t [s]"
set ylabel "Bubble Radius [m]"
set key top left
set size square
set grid

plot "OutputData-7" u 1:2 w l lw 2 t "LEVEL 7", \
     "OutputData-7" u 1:3 every 10 w p pt 8 ps 1.5 t "Analytical"
~~~

## References

~~~bib
@article{tanguy2014benchmarks,
  title={Benchmarks and numerical methods for the simulation of boiling flows},
  author={Tanguy, S{\'e}bastien and Sagan, Micha{\"e}l and Lalanne, Benjamin and Couderc, Fr{\'e}d{\'e}ric and Colin, Catherine},
  journal={Journal of Computational Physics},
  volume={264},
  pages={1--22},
  year={2014},
  publisher={Elsevier}
}
~~~
*/
