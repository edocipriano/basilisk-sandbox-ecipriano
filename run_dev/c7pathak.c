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

/**
## Simulation Setup

We use the centered solver with the evaporation source term
in the projection step. The extended velocity is obtained
from the doubled pressure-velocity coupling. We use the
evaporation model together with the multiomponent phase
change mechanism. */

#include "axi.h"
#include "navier-stokes/low-mach.h"
#define FILTERED 1
#if VELOCITY_JUMP
# include "two-phase-clsvof.h"
#else
# include "two-phase.h"
#endif
#include "tension.h"
#include "evaporation.h"
#include "view.h"

/**
### Boundary conditions

Outflow boundary conditions are set at the top and right
sides of the domain. */

u.n[top] = neumann (0.);
u.t[top] = neumann (0.);
p[top] = dirichlet (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);

/**
### Simulation Data

We declare the maximum and minimum levels of refinement,
the initial radius and diameter, and the radius from the
numerical simulation. */

int maxlevel, minlevel = 2;
double D0 = 5.e-6, effective_radius0 = 1.;
FILE * fp, * fpp;
bool restored = false;

int main (void) {

  NGS = 2, NLS = 1;

  /**
  We set the material properties of the fluids. The
  properties correspond to the n-heptane/nitrogen system
  at 28.6 bar. */

  rho1 = 626.7; rho2 = 17.51;
  mu1 = 1.e-3; mu2 = 1.e-5;
  Dmix1 = 0., Dmix2 = 6.77e-7;
  lambda1 = 0.1121, lambda2 = 0.04428;
  cp1 = 2505., cp2 = 1053.;
  dhev = 3.23e5;
  Pref = 2860000.;

  TL0 = 363., TG0 = 563.;

  nv = 2;

  /**
  We change the dimension of the domain as a function
  of the initial diameter of the droplet. */

  L0 = 2.*D0;

  /**
  We change the surface tension coefficient. and we
  decrease the tolerance of the Poisson solver. */

  f.sigma = 0.01;

  /**
  We run the simulation at different maximum
  levels of refinement. */

  for (maxlevel = 6; maxlevel <= 6; maxlevel++) {
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
  if (!restore (file = "restart")) {
    fraction (f, circle (x, y, 0.5*D0));
#if VELOCITY_JUMP
    foreach()
      d[] = circle (x, y, 0.5*D0);
    vertex scalar phi[];
    foreach_vertex()
      phi[] = (d[] + d[-1] + d[0,-1] + d[-1,-1])/4.;
    fractions (phi, f);
#endif
    effective_radius0 = pow(3.*statsf(f).sum, 1./3.);
  }
  else
    restored = true;

  /**
  We set the molecular weights of the chemial species
  involved in the simulation (by default inMW=1). */

  ThermoState tsl, tsg;
  tsl.T = TL0, tsl.P = Pref, tsl.x = (double[]){1.};
  tsg.T = TG0, tsg.P = Pref, tsg.x = (double[]){0.,1.};

  phase_set_thermo_state (liq, &tsl);
  phase_set_thermo_state (gas, &tsg);

  phase_set_properties (liq, MWs = (double[]){100.2});
  phase_set_properties (gas, MWs = (double[]){100.2,29.});

  /**
  The proper Antoine equation function must be set
  to the attribute *antoine* of the liquid phase mass
  fraction fields. */

  scalar lfuel = liq->YList[0];
  lfuel.antoine = antoine_heptane;

  /**
  We use the same boundary conditions used by
  [Pathak at al., 2018](#pathak2018steady). */

  scalar C7 = gas->YList[0];
  scalar N2 = gas->YList[1];
  scalar TG = gas->T;

  C7[top] = dirichlet (0.);
  C7[right] = dirichlet (0.);

  N2[top] = dirichlet (1.);
  N2[right] = dirichlet (1.);

  TG[top] = dirichlet (TG0);
  TG[right] = dirichlet (TG0);

  /**
  Set post-processing files. */

  char name[80];
  sprintf (name, "OutputData-%d", maxlevel);

  char name2[80];
  sprintf (name2, "Profiles-%d", maxlevel);

  if (restored && access (name, F_OK) == 0) {
    fp = fopen (name, "a");
    fpp = fopen (name2, "a");
  }
  else {
    fp = fopen (name, "w");
    fpp = fopen (name2, "w");
  }
}

/**
We adapt the grid according to the mass fractions of the
mass fraction of n-heptane, the temperature, and the
velocity field. */

#if TREE
event adapt (i++) {
  adapt_wavelet_leave_interface ({Y,T,u.x,u.y}, {f},
      (double[]){1.e-2,1.e-1,1.e-2,1.e-2}, maxlevel, minlevel, 1);
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
  double effective_radius = pow(3.*statsf(f).sum, 1./3.);
  double d_over_d02 = sq (effective_radius / effective_radius0);

  double TIntavg = avg_interface (liq_int->T, f);
  double YIntavg = avg_interface (gas_int->YList[0], f);
  double Tavg = avg_interface (T, f);
  double Yavg = avg_interface (gas_int->YList[0], f);

  fprintf (fp, "%g %g %g %g %g %g %g %g\n",
      t, effective_radius0, effective_radius, d_over_d02, TIntavg, YIntavg, Tavg, Yavg);
  fflush (fp);
}

/**
### Temperature and Mass Fraction Profiles

We write on a file the temperature and mass fraction
profiles at different time instants. */

event profiles (t = {3.29e-6, 3.e-5, 1.05e-4, 1.5e-4}) {
  for (double x = 0.; x < L0; x += 0.5*L0/(1 << maxlevel)) {
    fprintf (fpp, "%g %g %g\n",
        x, interpolate (T, x, 0.), interpolate (Y, x, 0.));
  }
  fprintf (fpp, "\n\n");
  fflush (fpp);
}

/**
### Movie

We write the animation with the evolution of the
n-heptane mass fraction, the interface position
and the temperature field. */

event movie (t += 2.e-6; t <= 1.6e-4) {
  clear();
  box();
  view (ty = -0.5, width=1200.);
  draw_vof ("f");
  squares ("Y", min = 0., max = 1., linear = true);
  mirror ({1.,0.}) {
    draw_vof ("f");
    squares ("T", min = statsf(T).min, max = TG0, linear = true);
  }
  save ("movie.mp4");
}

#if TEST_RESTART
event stop (t = 8.e-5) {
  dump ("restart");
  return 1;
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
     "OutputData-6" u 1:($3/2.5e-6)**2 w l lw 2 t "LEVEL 6"
~~~

~~~gnuplot Evolution of the interface temperature
reset
set xlabel "t [s]"
set ylabel "Interface Temperature [K]"
set key bottom right
set grid

plot "../data/pathak-heptane-T563-temp.csv" w p ps 2 t "Pathank et al., 2018", \
     "OutputData-6" u 1:5 w l lw 2 t "LEVEL 6"
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
