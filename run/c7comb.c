/**
# Non-Isothermal Evaporation of a n-Heptane Droplet

In this test case we want to simulate the evaporation of a
n-heptane droplet in nitrogen. The droplet is pure in a non
isothermal environment. Therefore, we want to solve both
the chemical species mass fractions and the temperature
field. At the beginning of the simulation, the droplet
is initialized at 363K, while the ambient temperature is 565K.
The heat conduction from the environment heats up the liquid
droplet increasing the vapor pressure value and, therefore,
increasing the evaporation rate of the droplet. The
interface temperature tends to a plateau, given by the
interplay between the heat conduction from the environment
and the evaporation process that cools down the interface.

![Evolution of the temperature field (left) and the n-heptane mass fraction (right)](c7pathak/movie.mp4)(height=400 width=900)
*/

#define NGS 6
#define NLS 1

/**
## Phase Change Setup

We define the number of gas and liquid species in the domain,
we filter the density field to reduce problems related with
the strong density ratio. The interface temperature must be
obtained from the root-finding procedure with the *fsolve()*
function, using the GSL interface. The Antoine equation is
used according to the thermodynamic equilibrium implemented
in [Pathak et al., 2018](#pathak2018steady), which proposed
this test case. */

#define VARPROP
#define SOLVE_TEMPERATURE
#define USE_GSL 0
#define USE_ANTOINE
#define POSSIBLE_BOILING

/**
## Simulation Setup

We use the centered solver with the evaporation source term
in the projection step. The extended velocity is obtained
from the doubled pressure-velocity coupling. We use the
evaporation model together with the multiomponent phase
change mechanism. */

#include "axi.h"
#include "navier-stokes/centered-evaporation.h"
#include "navier-stokes/centered-doubled.h"
#include "two-phase-varprop.h"
#include "opensmoke.h"
#include "opensmoke-properties.h"
#include "tension.h"
#include "evaporation-varprop.h"
#include "multicomponent-varprop.h"
//#include "chemistry.h"
//#include "spark.h"
#include "view.h"

/**
### Data for multicomponent model

We define the data required by the multicomponent phase
change mechanism, including the solution of the temperature
field. The equilibrium constant *inKeq* is ignored when
*USE_ANTOINE* or *USE_CLAPEYRON* is used. */

char* gas_species[NGS] = {"NC7H16", "O2", "CO2", "CO", "H2O", "N2"};
char* liq_species[NLS] = {"NC7H16"};
char* inert_species[1] = {"N2"};
double gas_start[NGS] = {0., 0.21, 0., 0., 0., 0.79};
double liq_start[NLS] = {1.};
double inDmix1[NLS] = {0.};
double inDmix2[NGS] = {6.77e-7,6.77e-7,6.77e-7,6.77e-7,6.77e-7,6.77e-7};
double inKeq[NLS] = {0.1};
double lambda1 = 0.1121;
double lambda2 = 0.04428;
double dhev = 3.23e5;
double cp1 = 2505.;
double cp2 = 1053.;
double TL0 = 300.;
double TG0 = 3000.;

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

TG[top] = dirichlet (TG0);
TG[right] = dirichlet (TG0);

/**
### Simulation Data

We declare the maximum and minimum levels of refinement,
the initial radius and diameter, and the radius from the
numerical simulation. */

int maxlevel, minlevel = 2;
double D0 = 0.7e-3, effective_radius0;

int main (void) {
  kinfolder = "two-step/n-heptane";

  /**
  We set the material properties of the fluids. The
  properties correspond to the n-heptane/nitrogen system
  at 28.6 bar. */

  rho1 = 626.7; rho2 = 17.51;
  mu1 = 1.e-3; mu2 = 1.e-5;
  Pref = 1e+5;

  /**
  We change the dimension of the domain as a function
  of the initial diameter of the droplet. */

  L0 = 10.*D0;

  /**
  We change the surface tension coefficient. and we
  decrease the tolerance of the Poisson solver. */

  f.sigma = 0.03;

  /**
  We run the simulation at different maximum
  levels of refinement. */

  for (maxlevel = 8; maxlevel <= 8; maxlevel++) {
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
  effective_radius0 = pow (3.*statsf(f).sum, 1./3.);

  /**
  We set the molecular weights of the chemial species
  involved in the simulation (by default inMW=1). */

  for (int jj=0; jj<NGS; jj++) {
    inMW[jj] = OpenSMOKE_MW (jj);
  }
  assert (NGS == OpenSMOKE_NumberOfSpecies());
  assert (NLS == OpenSMOKE_NumberOfLiquidSpecies());

  /**
  The proper Antoine equation function must be set
  to the attribute *antoine* of the liquid phase mass
  fraction fields. */

#ifdef USE_ANTOINE
  scalar YL_C7 = YLList[0];
  YL_C7.antoine = antoine_heptane;
#endif

#ifdef RADIATION
  divq_rad = optically_thin;
#endif

#ifdef SPARK
  spark.T = TG;
  spark.position = (coord){0.8*D0, 0.8*D0};
  spark.diameter = 0.1*D0;
  //spark.time = 0.05*sq(D0*1e3);
  spark.time = 0.;
  spark.duration = 0.02;
  spark.temperature = 4000.;
#endif
}

#define radial(r,r1,r2,T1,T2)(-r1*r2/(r*(r1-r2))*(T1 - T2) + (r1*T1 - r2*T2)/(r1-r2))

event end_init (i = 0) {
  scalar C7 = YGList[0];
  foreach() {
    double T1 = TL0, T2 = TG0;
    double r1 = 0.5*D0, r2 = 0.25*L0;
    double r = sqrt (sq(x) + sq(y));
    TG[] = (r <= r2) ? radial (r, r1, r2, T1, T2) : T2;
    C7[] = (r <= r2) ? radial (r, r1, r2, 1., 0.) : 0.;
    TG[] *= (1. - f[]);
    C7[] *= (1. - f[]);
  }
}

/**
We adapt the grid according to the mass fractions of the
mass fraction of n-heptane, the temperature, and the
velocity field. */

#if TREE
event adapt (i++) {
  scalar C7 = YList[0];
  adapt_wavelet_leave_interface ({C7,T,u.x,u.y}, {f},
      (double[]){1.e-2,1.e-2,1.e-2,1.e-2}, maxlevel, minlevel, 1);
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

  scalar YInt_c7 = YGIntList[0];
  scalar Y_c7 = YList[0];
  double effective_radius = pow (3.*statsf(f).sum, 1./3.);
  double d_over_d02 = sq (effective_radius / effective_radius0);

  double TIntavg = avg_interface (TInt, f);
  double YIntavg = avg_interface (YInt_c7, f);
  double Tavg = avg_interface (T, f);
  double Yavg = avg_interface (Y_c7, f);

  fprintf (fp, "%g %g %g %g %g %g %g\n",
      t/sq(2.*effective_radius0*1e3), effective_radius, d_over_d02, TIntavg, YIntavg, Tavg, Yavg);
  fflush (fp);
}

/**
### Movie

We write the animation with the evolution of the
n-heptane mass fraction, the interface position
and the temperature field. */

event movie (t += 0.01) {
  clear();
  box();
  view (tx = -0.5, ty = -0.5);
  draw_vof ("f");
  squares ("NC7H16", min = 0., max = 1.0, linear = true);
  save ("movie.mp4");
}

//#if DUMP
event snapshots (t += 0.001) {
  char name[80];
  sprintf (name, "snapshot-%g", t);
  dump (name);
}
//#endif

#if TRACE > 1
event profiling (i += 20) {
  static FILE * fp = fopen ("profiling", "w");
  trace_print (fp, 1);
}
#endif

event stop (t = 6) {
}

/**
## Results

The numerical results are compared with the results obtained
by [Pathak et al., 2018](#pathak2018steady) using a radially
symmetric model, integrated using an ODE solver.

~~~gnuplot Evolution of the squared diameter decay
reset
set xlabel "t [s]"
set ylabel "(D/D_0)^2"
set key top right
set grid

plot "../data/pathak-heptane-T563-diam.csv" w p ps 2 t "Pathank et al., 2018", \
     "OutputData-6" u 1:3 w l lw 2 t "LEVEL 6"
~~~

~~~gnuplot Evolution of the interface temperature
reset
set xlabel "t [s]"
set ylabel "Interface Temperature [K]"
set key bottom right
set grid

plot "../data/pathak-heptane-T563-temp.csv" w p ps 2 t "Pathank et al., 2018", \
     "OutputData-6" u 1:4 w l lw 2 t "LEVEL 6"
~~~

~~~gnuplot Evolution of the temperature profiles
reset
set xlabel "radius \mu m"
set ylabel "Temperature [K]"
set key bottom right
set grid

plot "../data/pathak-heptane-T563-Tprofile-329e-6.csv" w p pt 8 lc 1 t "time = 3.29x10^{-6} s", \
     "../data/pathak-heptane-T563-Tprofile-3e-5.csv"   w p pt 8 lc 2 t "time = 3.00x10^{-5} s", \
     "../data/pathak-heptane-T563-Tprofile-105e-4.csv" w p pt 8 lc 3 t "time = 1.05x10^{-4} s", \
     "../data/pathak-heptane-T563-Tprofile-150e-4.csv" w p pt 8 lc 4 t "time = 1.50x10^{-4} s", \
     "Profiles-6" index 0 u ($1*1e+6):2 w l lw 2 lc 1 notitle, \
     "Profiles-6" index 1 u ($1*1e+6):2 w l lw 2 lc 2 notitle, \
     "Profiles-6" index 2 u ($1*1e+6):2 w l lw 2 lc 3 notitle, \
     "Profiles-6" index 3 u ($1*1e+6):2 w l lw 2 lc 4 notitle
~~~

~~~gnuplot Evolution of the n-heptane mass fraction profiles
reset
set xlabel "radius \mu m"
set ylabel "Mass Fraction [-]"
set key top right
set grid

plot "../data/pathak-heptane-T563-Yprofile-329e-6.csv" w p pt 8 lc 1 t "time = 3.29x10^{-6} s", \
     "../data/pathak-heptane-T563-Yprofile-3e-5.csv"   w p pt 8 lc 2 t "time = 3.00x10^{-5} s", \
     "../data/pathak-heptane-T563-Yprofile-105e-4.csv" w p pt 8 lc 3 t "time = 1.05x10^{-4} s", \
     "../data/pathak-heptane-T563-Yprofile-150e-4.csv" w p pt 8 lc 4 t "time = 1.50x10^{-4} s", \
     "Profiles-6" index 0 u ($1*1e+6):3 w l lw 2 lc 1 notitle, \
     "Profiles-6" index 1 u ($1*1e+6):3 w l lw 2 lc 2 notitle, \
     "Profiles-6" index 2 u ($1*1e+6):3 w l lw 2 lc 3 notitle, \
     "Profiles-6" index 3 u ($1*1e+6):3 w l lw 2 lc 4 notitle
~~~

## References

~~~bib
@article{pathak2018steady,
  title={Steady-state and transient solutions to drop evaporation in a finite domain: Alternative benchmarks to the d2 law},
  author={Pathak, Ashish and Raessi, Mehdi},
  journal={International Journal of Heat and Mass Transfer},
  volume={127},
  pages={1147--1158},
  year={2018},
  publisher={Elsevier}
}
~~~
*/
