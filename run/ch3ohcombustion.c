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

//#include "ch3ohcombustion-input-multicomponent.h"

#define NGS 3
#define NLS 1

char* gas_species[NGS] = {"CH3OH", "N2", "O2"};
char* liq_species[NLS] = {"CH3OH"};
char* inert_species[1] = {"N2"};
double gas_start[NGS] = {0., 0.7670907862, 0.2329092138};
double liq_start[NLS] = {1.};
double inDmix1[NLS] = {0.};
double inDmix2[NGS] = {0., 0., 0.};
double inKeq[NLS] = {0.};

double lambda1 = 0.;
double lambda2 = 0.;
double dhev = 0.;
double cp1 = 0.;
double cp2 = 0.;

#define SOLVE_TEMPERATURE
#define USE_GSL 0
#define USE_ANTOINE_OPENSMOKE
#define FICK_CORRECTED
#define MOLAR_DIFFUSION
#define VARPROP
//#define MASS_DIFFUSION_ENTHALPY
//#define CHEMISTRY
//#define RADIATION
//#define GRAVIY

/**
## Simulation Setup

We use the centered solver with the evaporation source term
in the projection step. The extended velocity is obtained
from the doubled pressure-velocity coupling. We use the
evaporation model together with the multiomponent phase
change mechanism. */

#include "axi.h"
//#include "navier-stokes/velocity-jump.h"
#include "navier-stokes/centered-evaporation.h"
#include "navier-stokes/centered-doubled.h"
#include "opensmoke-properties.h"
#include "two-phase-varprop.h"
#ifdef GRAVIY
#include "pinning.h"
#endif
#include "tension.h"
#include "evaporation-varprop.h"
#include "multicomponent.h"
//#include "spark.h"
//#include "chemistry.h"
#include "view.h"

/**
### Data for multicomponent model

We define the data required by the multicomponent phase
change mechanism, including the solution of the temperature
field. The equilibrium constant *inKeq* is ignored when
*USE_ANTOINE* or *USE_CLAPEYRON* is used. */

double TL0 = 300.;
double TG0 = 500.;

/**
### Boundary conditions

Outflow boundary conditions are set at the top and right
sides of the domain. */

#ifdef VELOCITY_JUMP
u1.n[right] = neumann (0.);
u1.t[right] = neumann (0.);
u2.n[right] = neumann (0.);
u2.t[right] = neumann (0.);
p[right] = dirichlet (0.);
ps[right] = dirichlet (0.);
pg[right] = dirichlet (0.);

u1.n[top] = neumann (0.);
u1.t[top] = neumann (0.);
u2.n[top] = neumann (0.);
u2.t[top] = neumann (0.);
p[top] = dirichlet (0.);
ps[top] = dirichlet (0.);
pg[top] = dirichlet (0.);
#else
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
#endif

/**
### Simulation Data

We declare the maximum and minimum levels of refinement,
the initial radius and diameter, and the radius from the
numerical simulation. */

int maxlevel, minlevel = 2;
double D0 = 0.56e-3, effective_radius0, d_over_d02 = 1.;
double volumecorr = 0., volume0 = 0.;

int main (void) {
  //kinfolder = "skeletal/methanol";
  kinfolder = "evaporation/methanol-in-air";
  //kinfolder = "evaporation/ethanol-in-air";
  //kinfolder = "evaporation/n-heptane-in-air";
  //kinfolder = "evaporation/n-decane-in-air";

  /**
  We set the material properties of the fluids. The
  properties correspond to the n-heptane/nitrogen system
  at 28.6 bar. */

  mu1 = 1.e-3; mu2 = 1.e-5;
  rho1 = 0.; rho2 = 0.;
  Pref = 1*101325.;

  /**
  We change the dimension of the domain as a function
  of the initial diameter of the droplet. */

  double RR = 7.986462e+01;
  L0 = 0.5*RR*D0;

#ifdef GRAVIY
  double df = 0.1*D0;
  X0 = -0.5*L0;
  Y0 = 0.5*df;
  pinning.ap = 0.5*D0;
  pinning.ac = pinning.ap - 0.05*D0;
#endif

  /**
  We change the surface tension coefficient. and we
  decrease the tolerance of the Poisson solver. */

  f.sigma = 0.0227;
  TOLERANCE = 1.e-4;
  init_fields = true;

  /**
  We run the simulation at different maximum
  levels of refinement. */

  maxlevel = 10;
  init_grid (1 << (maxlevel - 2));
  run();
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

/**
We initialize the volume fraction field and we compute the
initial radius of the droplet. We don't use the value D0
because for small errors of initialization the squared
diameter decay would not start from 1. */

double mLiq0 = 0.;

event init (i = 0) {
  refine (circle (x, y, 2.*D0) > 0. && level < maxlevel);
  fraction (f, circle (x, y, 0.5*D0));
#ifdef GRAVIY
  effective_radius0 = pow(3./2.*statsf(f).sum, 1./3.);
#else
  effective_radius0 = pow(3.*statsf(f).sum, 1./3.);
#endif

  foreach (reduction(+:mLiq0))
    mLiq0 += rho1v[]*f[]*dv();

  /**
  We set the molecular weights of the chemial species
  involved in the simulation (by default inMW=1). */

  foreach_elem (YGList, jj) {
    inMW[jj] = OpenSMOKE_MW (jj);
    fprintf (stdout, "inMW[%d] = %g\n", jj, inMW[jj]);
  }

#ifdef SPARK
  spark.T = TG;
  spark.position = (coord){0., 0.7*D0};
  spark.diameter = 0.2*D0;
  //spark.time = 0.02;
  spark.time = 0.;
  spark.duration = 0.004;
  spark.temperature = 2700.;
#endif

//#ifdef RADIATION
//  divq_rad = optically_thin;
//#endif
}

//event set_spark (i++) {
//#if TREE
//  if (t <= spark.time)
//    refine (circle(x,y,D0+0.2*D0)*circle(x,y,D0-0.2*D0) < 0. && level < maxlevel);
//#endif
//  if (t >= spark.time && t <= (spark.time + spark.duration)) {
//    foreach() {
//      if (circle(x,y,0.7*D0+0.1*D0)*circle(x,y,0.7*D0-0.1*D0) < 0.)
//        TG[] = spark.baseline + (spark.temperature - spark.baseline)*(t - spark.time)/spark.duration;
//    }
//  }
//}

/**
We use the same boundary conditions used by
[Pathak at al., 2018](#pathak2018steady). */

//event bcs (i = 0) {
//  scalar CH3OH = YGList[OpenSMOKE_IndexOfSpecies ("CH3OH")];
//  scalar N2    = YGList[OpenSMOKE_IndexOfSpecies ("N2")];
//  scalar O2    = YGList[OpenSMOKE_IndexOfSpecies ("O2")];
//
//  CH3OH[top] = dirichlet (0.);
//  CH3OH[right] = dirichlet (0.);
//
//  N2[top] = dirichlet (0.7670907862);
//  N2[right] = dirichlet (0.7670907862);
//
//  O2[top] = dirichlet (0.2329092138);
//  O2[right] = dirichlet (0.2329092138);
//
//  TG[top] = dirichlet (TG0);
//  TG[right] = dirichlet (TG0);
//}

/**
We adapt the grid according to the mass fractions of the
mass fraction of n-heptane, the temperature, and the
velocity field. */

event adapt (i++) {
  scalar CH3OH = YList[OpenSMOKE_IndexOfSpecies ("CH3OH")];
  adapt_wavelet_leave_interface ({CH3OH,T,u.x,u.y}, {f},
      (double[]){1.e-1,1.e0,1.e-3,1.e-3}, maxlevel, minlevel, 1);
}

/**
We add the gravity contribution if the suspended droplet
configuration is considered. */

#ifdef GRAVIY
event acceleration (i++) {
  face vector av = a;
  foreach_face(x)
    av.x[] -= 9.81;
}
#endif

/**
## Post-Processing

The following lines of code are for post-processing purposes. */

/**
### Log file

Output variables to be monitored during the simulation. */

event logger (i++) {
  if (i == 0)
    fprintf (stdout, "i t min(u.x) max(u.x) min(u.y) max(u.y) min(T) max(T) min(F) max(F)\n");

  scalar FG = YGList[LSI[0]];

  double uxmin = statsf(u.x).min;
  double uxmax = statsf(u.x).max;
  double uymin = statsf(u.y).min;
  double uymax = statsf(u.y).max;
  double Tmin = statsf(T).min;
  double Tmax = statsf(T).max;
  double FGmin = statsf(FG).min;
  double FGmax = statsf(FG).max;

  fprintf (stdout, "%d %g %g %g %g %g %g %g %g %g\n",
      i, t, uxmin, uxmax, uymin, uymax, Tmin, Tmax, FGmin, FGmax);
}

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

#ifdef GRAVIY
  double effective_radius = pow(3./2.*statsf(f).sum, 1./3.);
#else
  double effective_radius = pow(3.*statsf(f).sum, 1./3.);
#endif
  d_over_d02 = sq (effective_radius / effective_radius0);

  double TIntavg = avg_interface (TInt, f, 0.1);
  double YIntavg = avg_interface (YInt_c7, f, 0.1);
  double mEvapavg = avg_interface (mEvapTot, f, 0.1);
  double Yavg = avg_interface (Y_c7, f, 0.1);

  double mLiq = 0.;
  foreach(reduction(+:mLiq))
    mLiq += rho1v[]*f[]*dv();

  fprintf (fp, "%g %g %g %g %g %g %g %g %g\n",
      t, t/sq(D0*1e3), effective_radius, d_over_d02, mLiq/mLiq0, TIntavg, YIntavg, mEvapavg, Yavg);
  fflush (fp);
}

/**
### Movie

We write the animation with the evolution of the
n-heptane mass fraction, the interface position
and the temperature field. */

event movie (t += 0.01; t <= 10) {
  clear();
  box();
  view (tx = -0.5, ty = -0.5);
  draw_vof ("f");
  squares ("T", min = TL0, max = statsf(T).max, linear = true);
  save ("movie.mp4");
}

event stop (i++) {
  if (d_over_d02 <= 0.05)
    return 1;
}

/**
### Snapshots

Output dump files for restore or post-processing. */

#if DUMP
event snapshots (t += 0.1) {
  char name[80];
  sprintf (name, "snapshots-%g", t);
  dump (name);
}
#endif

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
set xlabel "radius [m] x10^{6}"
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
set xlabel "radius [m] x10^{6}"
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
