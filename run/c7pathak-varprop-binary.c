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

#define NGS 3
#define NLS 2

#define VARPROP
#define FILTERED
#define SOLVE_TEMPERATURE
#define USE_GSL 0
#define USE_ANTOINE

/**
## Simulation Setup

We use the centered solver with the evaporation source term
in the projection step. The extended velocity is obtained
from the doubled pressure-velocity coupling. We use the
evaporation model together with the multiomponent phase
change mechanism. */

//#include "grid/multigrid.h"
#include "axi.h"
#include "navier-stokes/centered-evaporation.h"
#include "navier-stokes/centered-doubled.h"
#include "two-phase-varprop.h"
#include "opensmoke-properties.h"
#include "tension.h"
#include "icentripetal.h"
#include "reduced.h"
#include "evaporation-varprop.h"
#include "multicomponent-varprop.h"
#include "view.h"

/**
### Data for multicomponent model

We define the data required by the multicomponent phase
change mechanism, including the solution of the temperature
field. The equilibrium constant *inKeq* is ignored when
*USE_ANTOINE* or *USE_CLAPEYRON* is used. */

char* gas_species[NGS] = {"NC7H16", "NC16H34", "N2"};
char* liq_species[NLS] = {"NC7H16", "NC16H34"};
char* inert_species[1] = {"N2"};
double gas_start[NGS] = {0., 0., 1.};
double liq_start[NLS] = {0.5, 0.5};
double inDmix1[NLS] = {6.e-8, 6.e-8};
double inDmix2[NGS] = {6.77e-7, 6.77e-7,6.77e-7};
double inKeq[NLS] = {0.,0.};

double lambda1 = 0.1121;
double lambda2 = 0.04428;
double dhev = 3.23e5;
double cp1 = 2505.;
double cp2 = 1053.;
double TL0 = 300.;
double TG0 = 773.;

/**
### Boundary conditions

Outflow boundary conditions are set at the top and right
sides of the domain. */

u.n[left] = neumann (0.);
u.t[left] = neumann (0.);
p[left] = dirichlet (0.);
uext.n[left] = neumann (0.);
uext.t[left] = neumann (0.);
pext[left] = dirichlet (0.);

u.n[right] = dirichlet (0.);
u.t[right] = dirichlet (0.);
p[right] = neumann (0.);
uext.n[right] = dirichlet (0.);
uext.t[right] = dirichlet (0.);
pext[right] = neumann (0.);

/**
### Simulation Data

We declare the maximum and minimum levels of refinement,
the initial radius and diameter, and the radius from the
numerical simulation. */

int maxlevel, minlevel = 2;
//double D0 = 5.e-6, effective_radius0;
double D0 = 1.e-3, effective_radius0;

int main (void) {
  kinfolder = "evaporation/n-heptane-hexadecane-in-nitrogen";

  /**
  We set the material properties of the fluids. The
  properties correspond to the n-heptane/nitrogen system
  at 28.6 bar. */

  rho1 = 626.7; rho2 = 17.51;
  mu1 = 1.e-3; mu2 = 1.e-5;
  Pref = 1e+6;

  /**
  We change the dimension of the domain as a function
  of the initial diameter of the droplet. */

  L0 = 20.*D0;
  origin (-10.*D0);

  /**
  We change the surface tension coefficient. and we
  decrease the tolerance of the Poisson solver. */

  f.sigma = 0.03;
  //TOLERANCE = 1.e-6;

  /**
  We run the simulation at different maximum
  levels of refinement. */

  for (maxlevel = 9; maxlevel <= 9; maxlevel++) {
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
  effective_radius0 = pow(3.*statsf(f).sum, 1./3.);

  /**
  We set the molecular weights of the chemial species
  involved in the simulation (by default inMW=1). */

  for (int jj=0; jj<NGS; jj++) {
    inMW[jj] = OpenSMOKE_MW (jj);
  }

  /**
  The proper Antoine equation function must be set
  to the attribute *antoine* of the liquid phase mass
  fraction fields. */

  scalar YL_C7 = YLList[0];
  YL_C7.antoine = antoine_heptane;

  scalar YL_C10 = YLList[1];
  YL_C10.antoine = antoine_hexadecane;

#ifdef CENTRIPETAL
  sfm.p = (coord){0.,0.};
  sfm.eps = 1.e-4;
#endif
}

/**
We use the same boundary conditions used by
[Pathak at al., 2018](#pathak2018steady). */

//event bcs (i = 0) {
//  scalar C7 = YGList[0];
//  scalar C10 = YGList[1];
//  scalar N2 = YGList[2];
//
//  C7[top] = dirichlet (0.);
//  C7[right] = dirichlet (0.);
//
//  C10[top] = dirichlet (0.);
//  C10[right] = dirichlet (0.);
//
//  N2[top] = dirichlet (1.);
//  N2[right] = dirichlet (1.);
//
//  TG[top] = dirichlet (TG0);
//  TG[right] = dirichlet (TG0);
//}

/**
We adapt the grid according to the mass fractions of the
mass fraction of n-heptane, the temperature, and the
velocity field. */

#if TREE
event adapt (i++) {
  scalar C7 = YList[0];
  scalar C10 = YList[1];
  adapt_wavelet_leave_interface ({C7,C10,T,u.x,u.y}, {f},
      (double[]){1.e-3,1.e-3,1.e-1,1.e-3,1.e-3}, maxlevel, minlevel, 1);
}
#endif

event acceleration (i++) {
  face vector av = a;
  foreach_face(x)
    av.x[] -= 9.81;
}

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
  double effective_radius = pow(3.*statsf(f).sum, 1./3.);
  double d_over_d02 = sq (effective_radius / effective_radius0);

  double TIntavg = avg_interface (TInt, f);
  double YIntavg = avg_interface (YInt_c7, f);
  double Tavg = avg_interface (T, f);
  double Yavg = avg_interface (Y_c7, f);

  fprintf (fp, "%g %g %g %g %g %g %g\n",
      t, effective_radius, d_over_d02, TIntavg, YIntavg, Tavg, Yavg);
  fflush (fp);
}

///**
//### Temperature and Mass Fraction Profiles
//
//We write on a file the temperature and mass fraction
//profiles at different time instants. */
//
//event profiles (t = {3.29e-6, 3.e-5, 1.05e-4, 1.5e-4}) {
//  char name[80];
//  sprintf (name, "Profiles-%d", maxlevel);
//
//  /**
//  We create an array with the temperature and mass
//  fraction profiles for each processor. */
//
//  scalar C7 = YList[0];
//
//  Array * arrtemps = array_new();
//  Array * arrmassf = array_new();
//  for (double x = 0.; x < L0; x += 0.5*L0/(1 << maxlevel)) {
//    double valt = interpolate (T, x, 0.);
//    double valm = interpolate (C7, x, 0.);
//    valt = (valt == nodata) ? 0. : valt;
//    valm = (valm == nodata) ? 0. : valm;
//    array_append (arrtemps, &valt, sizeof(double));
//    array_append (arrmassf, &valm, sizeof(double));
//  }
//  double * temps = (double *)arrtemps->p;
//  double * massf = (double *)arrmassf->p;
//
//  /**
//  We sum each element of the arrays in every processor. */
//
//  @if _MPI
//  int size = arrtemps->len/sizeof(double);
//  MPI_Allreduce (MPI_IN_PLACE, temps, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//  MPI_Allreduce (MPI_IN_PLACE, massf, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//  @endif
//
//  /**
//  The master node writes the profiles on a file. */
//
//  if (pid() == 0) {
//    static FILE * fpp = fopen (name, "w");
//    int count = 0;
//    for (double x = 0.; x < L0; x += 0.5*L0/(1 << maxlevel)) {
//      fprintf (fpp, "%g %g %g\n", x, temps[count], massf[count]);
//      count++;
//    }
//    fprintf (fpp, "\n\n");
//    fflush (fpp);
//  }
//  array_free (arrtemps);
//  array_free (arrmassf);
//}

/**
### Movie

We write the animation with the evolution of the
n-heptane mass fraction, the interface position
and the temperature field. */

//event movie (t += 2.e-6) {
event movie (t += 0.001) {
  clear();
  box();
  view (ty = -0.5, width=1200.);
  draw_vof ("f");
  squares ("NC7H16", min = 0., max = 0.5, linear = true);
  mirror ({1.,0.}) {
    draw_vof ("f");
    squares ("NC16H34", min = 0., max = 0.9, linear = true);
  }
  save ("movie.mp4");
}

#if DUMP
event snapshots (t += 1.e-5) {
  char name[80];
  sprintf (name, "snapshot-%g", t);
  dump (name);
}
#endif

#if TRACE > 1
event profiling (i += 20) {
  static FILE * fp = fopen ("profiling", "w");
  trace_print (fp, 1);
}
#endif

//event stop (t = 1.6e-4) {
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
