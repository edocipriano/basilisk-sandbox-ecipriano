/**
# Expansion of a Liquid Droplet

In this test case, we study the thermal expansion of a liquid droplet of
n-heptane, at different ambient temperatures and at different levels of
refinement. The aim is to evaluate the convergence of the
[multicomponent.h](multicomponent.h) solver with variable properties.

The evaporation module is used suppressing the phase change, in order to focus
on the thermal expansion only, and to avoid evaporation.

![Evolution of the temperature field](expansion/movie.mp4)
*/

/**
## Phase Change Setup

We define the number of gas and liquid species in the domain,
and we initialize all the properties necessary for the multicomponent phase
change model. The properties are set to null values because they are
overwritten by the variable properties formulation, which computes all the
physical properties as a function of the thermodynamic state of the mixture. */

#define NGS 2
#define NLS 1

char* gas_species[NGS] = {"NC7H16", "N2"};
char* liq_species[NLS] = {"NC7H16"};
char* inert_species[1] = {"N2"};
double gas_start[NGS] = {0., 1.};
double liq_start[NLS] = {1.};
double inDmix1[NLS] = {0.};
double inDmix2[NGS] = {0., 0.};
double inKeq[NLS] = {0.};

double lambda1 = 0.;
double lambda2 = 0.;
double dhev = 0.;
double cp1 = 0.;
double cp2 = 0.;

/**
We set the initial temperature of the liquid and of the gas phase. */

double TL0 = 300.;
double TG0 = 350.;

/**
We solve the temperature field, with variable properties, and we reduce the
tolerance for the calculation of the variable properties. The interfacial
temperature is not computed from the jump conditon, it is just set to the gas
phase temperature value. */

#define SOLVE_TEMPERATURE
#define NO_ADVECTION_DIV 1
#define T_PROP 1.e-6
#define FIXED_INTERFACE_TEMPERATURE TG0

/**
## Simulation Setup

We use the centered solver with the evaporation source term
in the projection step. The extended velocity is obtained
from the doubled pressure-velocity coupling. We use the
evaporation model together with the multiomponent phase
change mechanism.

We use the centered solver with the divergence source term in the projection
step. The calculation of the extended velocity can be skipped, because no phase
change is present. OpenSMOKE++ is used for the variable properties calculation. */

#include "axi.h"
#include "navier-stokes/centered-evaporation.h"
#include "navier-stokes/centered-doubled.h"
#include "opensmoke-properties.h"
#include "two-phase.h"
#include "tension.h"
#include "evaporation.h"
#include "multicomponent-varprop.h"
#include "view.h"

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
numerical simulation, and additional data for post-processing. */

int maxlevel, minlevel = 2;
double D0 = 1e-3, effective_radius0, d_over_d02 = 1.;
double volumecorr = 0., volume0 = 0., rho0 = 0., rhof = 0.;

int main (void) {
  /**
  We set the kinetics folder, which defines the species
  of the simulation, and it is used by OpenSMOKE++ for the
  calculation of the thermodynamic and transport properties. */

  kinfolder = "evaporation/n-heptane-in-nitrogen";

  /**
  We set additional data for the simulation. */

  rho1 = 1.; rho2 = 1.;
  mu1 = 1.; mu2 = 1.;
  Pref = 10*101325.;

  /**
  We change the dimension of the domain as a function
  of the initial diameter of the droplet. */

  L0 = 1.5*D0;

  /**
  We change the surface tension coefficient. */

  f.sigma = 0.03;

  /**
  We run the simulation at different maximum
  levels of refinement. */

  for (maxlevel = 5; maxlevel <= 8; maxlevel++) {
    init_grid (1 << maxlevel);
    run();
  }
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

/**
We initialize the volume fraction field and we compute the
initial radius of the droplet. We don't use the value D0
because for small errors of initialization the squared
diameter decay would not start from 1. */

event init (i = 0) {
  fraction (f, circle (x, y, 0.5*D0));

  /**
  We compute initial variables useful for post-processing. */

  effective_radius0 = pow (3.*statsf(f).sum, 1./3.);
  volume0 = 4./3.*pi*pow (effective_radius0, 3.);

  foreach (reduction(+:mLiq0))
    mLiq0 += rho1v[]*f[]*dv();

  /**
  We set the molecular weights of the chemial species
  involved in the simulation (by default inMW=1). */

  foreach_elem (YGList, jj) {
    inMW[jj] = OpenSMOKE_MW (jj);
    fprintf (stdout, "inMW[%d] = %g\n", jj, inMW[jj]);
  }

  /**
  We compute variables for the analytical steady-state solution. */

  ThermoState ts0;
  ts0.T = TL0;
  ts0.P = Pref;
  ts0.x = (double[]){1.};

  ThermoState tsf;
  tsf.T = TG0;
  tsf.P = Pref;
  tsf.x = (double[]){1.};

  rho0 = tp1.rhov (&ts0);
  rhof = tp1.rhov (&tsf);
}

/**
We use the same boundary conditions used by
[Pathak at al., 2018](#pathak2018steady). */

event bcs (i = 0) {
  scalar NC7H16 = YGList[OpenSMOKE_IndexOfSpecies ("NC7H16")];
  scalar N2    = YGList[OpenSMOKE_IndexOfSpecies ("N2")];

  NC7H16[top] = dirichlet (0.);
  NC7H16[right] = dirichlet (0.);

  N2[top] = dirichlet (1.);
  N2[right] = dirichlet (1.);

  TG[top] = dirichlet (TG0);
  TG[right] = dirichlet (TG0);
}

/**
We adapt the grid according to the mass fractions of the
mass fraction of n-heptane, the temperature, and the
velocity field. */

#if TREE
event adapt (i++) {
  adapt_wavelet_leave_interface ({T,u.x,u.y}, {f},
      (double[]){1.e-2,1.e-1,1.e-1}, maxlevel, minlevel, 1);
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

  double effective_radius = pow (3.*statsf(f).sum, 1./3.);
  d_over_d02 = sq (effective_radius / effective_radius0);

  double mLiq = 0.;
  foreach(reduction(+:mLiq))
    mLiq += rho1v[]*f[]*dv();

  double exact_volume = volume0*rho0/rhof;
  double exact_radius = pow (3./4./pi*exact_volume, 1./3.);
  double exact_d_over_d02 = sq (exact_radius/effective_radius0);
  double relerr = fabs (effective_radius - exact_radius)/effective_radius;

  fprintf (fp, "%g %g %g %g %g %g %g\n", t, t/sq(D0*1e3), effective_radius,
      d_over_d02, mLiq/mLiq0, relerr, exact_d_over_d02);
}

/**
### Movie

We write the animation with the evolution of the
n-heptane mass fraction, the interface position
and the temperature field. */

event pictures (t = 0.1) {
  clear();
  view (ty = -0.5);
  draw_vof ("f", filled = -1, fc = {1.,1.,1.});
  squares ("T", min = TL0, max = TG0, linear = true);
  isoline ("T", n = 12., min = TL0, max = TG0);
  mirror ({1.,0.}) {
    draw_vof ("f", filled = -1, fc = {1.,1.,1.});
    squares ("T", min = TL0, max = TG0, linear = true);
    isoline ("T", n = 12., min = TL0, max = TG0);
  }
  save ("temperature.png");

  scalar drhodtplot[];
  foreach()
    drhodtplot[] = drhodtext[]*f[];

  clear();
  view (ty = -0.5);
  draw_vof ("f", filled = -1, fc = {1.,1.,1.});
  squares ("drhodtplot", spread = -1, linear = true);
  isoline ("drhodtplot", n = 12.,
      min = statsf(drhodtplot).min, max = statsf(drhodtplot).max);
  mirror ({1.,0.}) {
    draw_vof ("f", filled = -1, fc = {1.,1.,1.});
    squares ("drhodtplot", spread = -1, linear = true);
    isoline ("drhodtplot", n = 12.,
        min = statsf(drhodtplot).min, max = statsf(drhodtplot).max);
  }
  save ("divergence.png");

  scalar un[];
  foreach()
    un[] = norm (u)*f[];
  clear();
  view (ty = -0.5);
  draw_vof ("f", filled = -1, fc = {1.,1.,1.});
  squares ("un", spread = -1, linear = true);
  isoline ("un", n = 12., min = statsf(un).min, max = statsf(un).max);
  mirror ({1.,0.}) {
    draw_vof ("f", filled = -1, fc = {1.,1.,1.});
    squares ("un", spread = -1, linear = true);
    isoline ("un", n = 12., min = statsf(un).min, max = statsf(un).max);
  }
  save ("velocity.png");
}

event movie (t += 0.01; t <= 3) {
  if (maxlevel == 6) {
    clear();
    box();
    view (tx = -0.5, ty = -0.5);
    draw_vof ("f");
    squares ("T", min = TL0, max = statsf(T).max, linear = true);
    save ("movie.mp4");
  }
}

/**
## Results

~~~gnuplot Expansion of the Square Diameter
reset
set xr[0:3]
set grid
set key bottom right

plot "../data/expansion.sol1" u 2:5 every 2 w p pt 2 t "1D Model", \
     "OutputData-7" u 2:7 w l lc black dt 2 t "Exact", \
     "OutputData-5" u 2:4 w l t "LEVEL 5", \
     "OutputData-6" u 2:4 w l t "LEVEL 6", \
     "OutputData-7" u 2:4 w l t "LEVEL 7", \
     "OutputData-8" u 2:4 w l t "LEVEL 8"
~~~

~~~gnuplot Relative Errors
reset

stats "OutputData-5" using (last5=$6) nooutput
stats "OutputData-6" using (last6=$6) nooutput
stats "OutputData-7" using (last7=$6) nooutput
stats "OutputData-8" using (last8=$6) nooutput

set print "Errors.csv"

print sprintf ("%d %.12f", 2**5, last5)
print sprintf ("%d %.12f", 2**6, last6)
print sprintf ("%d %.12f", 2**7, last7)
print sprintf ("%d %.12f", 2**8, last8)

unset print

reset
set xlabel "N"
set ylabel "Relative Error"

set logscale x 2
set logscale y

set xr[2**4:2**9]
#set yr[1e-4:10]

set size square
set grid

plot "Errors.csv" w p pt 8 ps 2 title "Results", \
  50*x**(-1) lw 2 title "1^{st} order", \
  1000*x**(-2) lw 2 title "2^{nd} order"
~~~
*/

