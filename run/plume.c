/**
# Thermal Expansion of a Liquid Droplet

In this test case, we study the thermal expansion of a liquid droplet of
n-heptane, at different ambient temperatures and at different levels of
refinement. The aim is to evaluate the convergence of the
[multicomponent.h](multicomponent.h) solver with variable properties.

The evaporation module is used suppressing the phase change, in order to focus
on the thermal expansion only, and to avoid evaporation.

![Evolution of the temperature field](expansion/movie.mp4)
*/

/**
## Default Simulation Data

The following data can be overwritten using compilation flags
in order to study the sensitivity to these parameters running
different simulations in parallel. */

/**
## Phase Change Setup

We define the number of gas and liquid species in the domain,
and we initialize all the properties necessary for the multicomponent phase
change model. The properties are set to null values because they are
overwritten by the variable properties formulation, which computes all the
physical properties as a function of the thermodynamic state of the mixture. */

#include "navier-stokes/low-mach.h"
#include "navier-stokes/perfs.h"
#include "opensmoke/properties.h"
#include "combustion.h"
#include "view.h"
#include "custom-cmaps.h"

/**
### Boundary conditions

Outflow boundary conditions are set at the top and right
sides of the domain. */

u.n[top] = neumann (0.);
u.t[top] = neumann (0.);
p[top] = dirichlet (0.);

/**
### Simulation Data

We declare the maximum and minimum levels of refinement,
the initial radius and diameter, and the radius from the
numerical simulation, and additional data for post-processing. */

int maxlevel = 8, minlevel = 2;
double D0 = 0.025, TMAX = 500.;

face vector av[];

scalar d[];

int main (void) {
  /**
  We set the kinetics folder, which defines the species
  of the simulation, and it is used by OpenSMOKE++ for the
  calculation of the thermodynamic and transport properties. */

  kinetics ("evaporation/n-heptane-in-nitrogen", &NS);
  gas_species = new_species_names (NS);

  /**
  We set additional data for the simulation. */

  rhoval = 1., muval = 1.e-2;
  lambdaval = 0.06, cpval = 1050.;
  Dmixval = 1.e-2, MWval = 28.;

  T0 = 300., Pref = 101325.;

  /**
  We change the dimension of the domain as a function
  of the initial diameter of the droplet. */

  L0 = 4.*D0;
  a = av;

  pcm.isomassfrac = true;
  pcm.divergence = true;

  /**
  We run the simulation at different maximum
  levels of refinement. */

  init_grid (1 << maxlevel);
  run();

  free_species_names (NS, gas_species), gas_species = NULL;
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

/**
We initialize the volume fraction field and we compute the
initial radius of the droplet. We don't use the value D0
because for small errors of initialization the squared
diameter decay would not start from 1. */

event init (i = 0) {
  ThermoState ts;
  ts.T = T0, ts.P = Pref, ts.x = (double[]){0.,1.};
  phase_set_thermo_state (gas, &ts);

  phase_set_properties (gas, MWs = (double[]){100.2,29.});

  scalar f1[], T = gas->T;
  fraction (f1, circle (x - 0.5*L0, y - 0.25*L0, 0.25*D0));
  foreach()
    T[] = T0*(1. - f1[]) + TMAX*f1[];

  T.gradient = zero;
}

event logger (i++) {
  fprintf (stdout, "%g %g %g\n", t, dt, statsf(u.y).max);
}

event acceleration (i++) {
  coord G = {0, -10., 0.};
  coord Z = {0.,0.,0.};
  scalar rho = gas->rho;
  foreach_face() {
    coord o = {x,y,z};
    double phia = 0.;
    foreach_dimension()
      phia += (o.x - Z.x)*G.x;
    a.x[] -= alpha.x[]/(fm.x[] + SEPS)*phia*(rho[] - rho[-1])/Delta;
  }
}

/**
We adapt the grid according to the mass fractions of the
mass fraction of n-heptane, the temperature, and the
velocity field. */

#if TREE
event adapt (i++) {
  scalar Y1 = gas->YList[0];
  adapt_wavelet ({Y1,gas->T,u.x,u.y},
      (double[]){1.e-2,1.e-1,1.e-1,1.e-1,1.e-1}, maxlevel);
}
#endif

/**
## Post-Processing

The following lines of code are for post-processing purposes. */

event movie (t += 0.01; t <= 1.6) {
  clear();
  view (tx = -0.5, ty = -0.5);
  squares ("T", min = (2.*T0-TMAX), max = TMAX,
      linear = true, map = blue_white_red);
  isoline ("d", lc = {255/255., 127./255., 14/255.});
  save ("movie.mp4");
}

