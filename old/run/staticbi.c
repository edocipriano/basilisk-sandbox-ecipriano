/**
# Isothermal Evaporation of a Static Binary Droplet

A binary liquid droplet is placed on the lower-left edge of
the domain. The two chemical species in liquid phase have
the same physical properties, but different volatility. The
relative volatility between the heavy and the light species
is equal to 0.5. Therefore, we expect the light component to
start to evaporate first, increasing its mass fractions in
gas phase and decreasing its concentration in the liquid
phase. The heavy component accumulates in the liquid phase
in response to the evaporation of the light species. The
gas phase is initially full of an inert compound, important
in combustion simulations, which always remains in gas
phase.

![Evolution of the mass fraction fields of the light, heavy, and inert components, and the grid refinement](staticbi/movie.mp4)(height=600 width=600)
*/

/**
## Phase Change Setup

We let the default settings of the evaporation model: the
Stefan flow is shifted toward the liquid phase, and the
consistent phase is the liquid, which is advected using the
extended velocity. The multicomponent model requires the
number of gas and liquid species to be set as compiler
variables. We don't need to solve the temperature field
because the vapor pressure is set to a constant value,
different for each chemical species. Using GSL at level
1 we can activate the coupled solution of the interface jump
condition. */

#define NGS 3
#define NLS 2

/**
## Simulation Setup

We use the centered solver with the evaporation source term
in the projection step. The extended velocity is obtained
from the doubled pressure-velocity coupling. We use the
evaporation model together with the multiomponent phase
change mechanism. */

#include "navier-stokes/centered-evaporation.h"
#include "navier-stokes/centered-doubled.h"
#include "two-phase.h"
#include "tension.h"
#include "evaporation.h"
#include "multicomponent.h"
#include "balances.h"
#include "view.h"

/**
### Data for multicomponent model

We define the lists with the names of the chemical species
in gas and in liquid phase. The initial mass fractions are
defined for each component, as well as the diffusion
coefficients and the thermodynamic equilibrium constant.
Since we set the vapor pressure to a constant value, we
don't need to solve the temperature field. */

char* gas_species[NGS] = {"A", "B", "C"};
char* liq_species[NLS] = {"A", "B"};
char* inert_species[1] = {"C"};
double gas_start[NGS] = {0.0, 0.0, 1.0};
double liq_start[NLS] = {0.5, 0.5};
double inDmix1[NLS] = {4.e-6, 4.e-6};
double inDmix2[NGS] = {8.e-5, 8.e-5, 8.e-5};
double inKeq[NLS] = {0.8, 0.4};

/**
### Boundary conditions

Outflow boundary conditions are set at the top and right
sides of the domain. */

u.n[top] = neumann (0.);
u.t[top] = neumann (0.);
p[top] = dirichlet (0.);
uext.n[top] = neumann (0.);
uext.t[top] = neumann (0.);
pext[top] = dirichlet (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);
uext.n[right] = neumann (0.);
uext.t[right] = neumann (0.);
pext[right] = dirichlet (0.);

/**
### Simulation Data

We declare the maximum and minimum levels of refinement,
the initial radius and diameter, and the radius from the
numerical simulation. */

int maxlevel, minlevel = 5;
double D0 = 0.4e-3, effective_radius0;

int main (void) {
  /**
  We set the material properties of the fluids. */

  rho1 = 10.; rho2 = 1.;
  mu1 = 1.e-4; mu2 = 1.e-5;

  /**
  We change the dimension of the domain as a function
  of the initial diameter of the droplet. */

  D0 = 0.4e-3; L0 = 4.*D0;

  /**
  We change the surface tension coefficient. */

  f.sigma = 0.03;

  /**
  We run the simulation at three different levels
  of refinement. */

  for (maxlevel = 7; maxlevel <= 7; maxlevel++) {
    //CFL = 0.1;
    init_grid (1 << maxlevel);
    run ();
  }
}

#define circle(x, y, R) (sq(R) - sq(x) - sq(y))

/**
We initialize the volume fraction field and we compute the
initial radius of the droplet. We don't use the value D0
because for small errors of initialization the squared
diameter decay would not start from 1. */

event init (i = 0) {
  fraction (f, circle (x, y, 0.5*D0));
  effective_radius0 = sqrt (4./pi*statsf(f).sum);

#ifdef BALANCES
  mb.liq_species = liq_species;
  mb.gas_species = gas_species;
  mb.YLList = YLList;
  mb.YGList = YGList;
  mb.mEvapList = mEvapList;
  mb.liq_start = liq_start;
  mb.gas_start = gas_start;
  mb.rho1 = rho1;
  mb.rho2 = rho2;
  mb.inDmix1 = inDmix1;
  mb.inDmix2 = inDmix2;
  mb.maxlevel = maxlevel;
  mb.boundaries = true;
#endif
}

/**
We overwrite the boundary conditions for the chemical
species mass fractions in the *bcs* event. */

event bcs (i = 0) {
  scalar YLA = YLList[0], YLB = YLList[1];
  scalar YGA = YGList[0], YGB = YGList[1], YGC = YGList[2];

  YLA[top] = dirichlet (0.);
  YLB[top] = dirichlet (0.);
  YLA[right] = dirichlet (0.);
  YLB[right] = dirichlet (0.);
  YGA[top] = dirichlet (0.);
  YGB[top] = dirichlet (0.);
  YGC[top] = dirichlet (1.);
  YGA[right] = dirichlet (0.);
  YGB[right] = dirichlet (0.);
  YGC[right] = dirichlet (1.);
}

/**
We adapt the grid according to the mass fractions of the
species A and B, the velocity and the interface position. */

#if TREE
event adapt (i++) {
  scalar YA = YList[0], YB = YList[1];
  adapt_wavelet_leave_interface ({YA,YB,u.x,u.y}, {f},
      (double[]){1.e-4,1.e-4,1.e-3,1.e-3}, maxlevel, minlevel, 1);
}
#endif

/**
## Post-Processing

The following lines of code are for post-processing purposes. */

/**
### Output Files

We write on a file the squared diameter decay and the
dimensionless time. */

event output_data (i++) {
  char name[80];
  sprintf (name, "OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, "w");

  double effective_radius = sqrt (4./pi*statsf(f).sum);
  double tad = t*inDmix2[0]/sq (2.*effective_radius0);
  double d_over_d0  = effective_radius / effective_radius0;
  double d_over_d02 = sq (d_over_d0);

  fprintf (fp, "%g %g %g %g\n", t, tad, d_over_d0, d_over_d02);
  fflush (fp);
}

/**
### Movie

We write the animation with the evolution of the
chemical species, the interface position and the
grid refinement. */

event movie (t += 2.e-5; t <= 0.005) {
  clear();
  draw_vof ("f", lw = 1.5);
  squares ("B", linear = true, min = 0., max = 0.56);
  mirror ({0,1}) {
    draw_vof ("f", lw = 1.5);
    squares ("C", linear = true, min = 0., max = 1.);
  }
  mirror ({1,0}) {
    draw_vof ("f", lw = 1.5);
    squares ("A", linear = true, min = 0., max = liq_start[0]);
  }
  mirror ({1,1}) {
    cells ();
    draw_vof ("f", lw = 1.5);
  }
  save ("movie.mp4");
}

#if DUMP
event snapshot (t += 2.e-4) {
  char name[80];
  sprintf (name, "snapshot-%g", t);
  dump (name);
}
#endif

/**
## Results

~~~gnuplot Squared Diameter Decay
reset
set xlabel "t [s]"
set ylabel "(D/D_0)^2"
set key top right
set size square
set grid

plot "OutputData-7" u 2:4 w l lw 2 t "LEVEL 7"
~~~

The conservation tests compare the mass of the chemical species in
liquid phase with the total amount of the same species that
evaporates. If the global conservation is considered, the volume
fraction is used instead of the mass fraction field. See
[balances.h](../src/balances.h) for details.

~~~gnuplot Liquid Phase Mass Conservation
reset
set xlabel "t [s]"
set ylabel "(m_L - m_L^0) [kg]"
set key top right
set size square
set grid

plot "balances-7" every 500 u 1:10 w p ps 1.2 lc 1 title "Evaporated Mass Species A", \
     "balances-7" every 500 u 1:11 w p ps 1.2 lc 2 title "Evaporated Mass Species B", \
     "balances-7" every 500 u 1:4  w p ps 1.2 lc 3 title "Evaporated Mass Total", \
     "balances-7" u 1:(-$5) w l lw 2 lc 1 title "Variation Mass Species A", \
     "balances-7" u 1:(-$6) w l lw 2 lc 2 title "Variation Mass Species B", \
     "balances-7" u 1:(-$2) w l lw 2 lc 3 title "Variation Mass Total"
~~~

~~~gnuplot Gas Phase Mass Conservation
reset
set xlabel "t [s]"
set ylabel "(m_G - m_G^0) [kg]"
set key top left
set size square
set grid

plot "balances-7" every 500 u 1:(-$10) w p ps 1.2 lc 1 title "Evaporated Mass Species A", \
     "balances-7" every 500 u 1:(-$11) w p ps 1.2 lc 2 title "Evaporated Mass Species B", \
     "balances-7" every 500 u 1:(-$4)  w p ps 1.2 lc 3 title "Evaporated Mass Total", \
     "balances-7" u 1:7 w l lw 2 lc 1 title "Variation Mass Species A", \
     "balances-7" u 1:8 w l lw 2 lc 2 title "Variation Mass Species B", \
     "balances-7" u 1:3 w l lw 2 lc 3 title "Variation Mass Total"
~~~
*/

