/**
# Pressurization of a Closed Tank

When the domain is closed, the volume of the system cannot change, and the
expansion of the fluid must be balanced by an increase of the thermodynamic
pressure $P_0$, which is uniform in space for a low Mach number system. This
test case verifies the pressurization model implemented in
[low-mach.h](../src/navier-stokes/low-mach.h), using a configuration which has
an analytic solution.

A closed box is filled with a liquid layer at the bottom and with an ideal gas
ullage at the top. A uniform volumetric heat source $\dot{q}$ is applied to the
ullage only, while the thermal conductivity is set to zero in both phases, in
order to keep the temperature of the ullage uniform. The liquid is
incompressible, therefore the volume of the ullage does not change and the
system behaves as a constant-volume heating of an ideal gas:
$$
  \dfrac{dP_0}{dt} = \left(\gamma - 1\right)\dot{q}
$$
$$
  T_g(t) = T_g^0 + \dfrac{\gamma\dot{q}t}{\rho_g^0 c_{p,g}}
$$
Since the expansion of the gas is exactly balanced by the compression due to the
pressure rise, the velocity field must remain zero everywhere. This is a strict
test of the compressibility-weighted projection: if the volume variation was
redistributed uniformly over the domain, the incompressible liquid would be
forced to expand, generating a spurious velocity field.

![Evolution of the gas phase temperature](pressurization/movie.mp4)
*/

/**
## Simulation Setup

We use the low Mach number Navier--Stokes solver, which includes the volumetric
expansion terms in the projection step, together with the variable properties
module. The phase change is switched off using the fixed flux model with a null
vaporization rate, in order to focus exclusively on the pressurization.
*/

#include "navier-stokes/low-mach.h"
#define P_ERR 1.e-6
#include "two-phase-varprop.h"
#include "basilisk-properties.h"
#include "two-phase.h"
#include "tension.h"
#include "phasechange.h"
#include "fixedflux.h"
#include "view.h"

/**
### Model Data

The molecular weight of the gas phase, the volumetric heat source applied to the
ullage, and the initial liquid level. */

int maxlevel, minlevel = 2;
double q0 = 1.e5, level0 = 0.5, tend = 0.1;
double gamma0, rhog0, P0exact, Tg0exact;

#define MWG 29.

/**
The density of the gas phase follows the ideal gas law, which makes the ullage
compressible and sensitive to the variation of the thermodynamic pressure. */

double gasprop_density_idealgas (void * p) {
  ThermoState * ts = p;
  return ts->P*MWG/(R_GAS*1.e3*ts->T);
}

/**
### Boundary Conditions

No boundary condition is imposed: the default symmetry conditions of Basilisk
make the box closed and adiabatic, which is exactly the configuration we want to
reproduce. */

int main (void) {

  /**
  We use a single chemical species in each phase. */

  NGS = 1, NLS = 1;

  /**
  We set the material properties. The thermal conductivity and the diffusivity
  are null, in order to keep the ullage uniform and to obtain an analytic
  solution. The gas phase density is overwritten by the ideal gas law. */

  rho1 = 1000., rho2 = 1.;
  mu1 = 1.e-3, mu2 = 1.8e-5;
  Dmix1 = 0., Dmix2 = 0.;
  lambda1 = 0., lambda2 = 0.;
  cp1 = 4184., cp2 = 1004.5;
  dhev = 0.;

  Pref = 101325.;
  TG0 = 300., TL0 = 300.;

  /**
  The system is closed: the net expansion is converted into a variation of the
  thermodynamic pressure instead of leaving the domain. */

  closed = true;

  /**
  The composition of the two phases does not change, and the phase change is
  switched off setting a null vaporization rate. */

  pcm.isomassfrac = true;
  pcm.divergence = true;
  mEvapVal = 0.;

  f.sigma = 0.;

  L0 = 0.01;
  DT = 1.e-3;

  /**
  The heat capacity ratio of the ullage is used by the analytic solution. */

  gamma0 = cp2/(cp2 - R_GAS*1.e3/MWG);

  for (maxlevel = 4; maxlevel <= 6; maxlevel++) {
    init_grid (1 << maxlevel);
    run();
  }
}

/**
The simulation is stopped at `tend`. */

event stop (t = tend);

/**
We initialize a flat interface which splits the domain into a liquid layer and a
gas ullage. */

event init (i = 0) {
  fraction (f, level0*L0 - y);

  ThermoState tsl, tsg;
  tsl.T = TL0, tsl.P = Pref, tsl.x = (double[]){1.};
  tsg.T = TG0, tsg.P = Pref, tsg.x = (double[]){1.};

  phase_set_thermo_state (liq, &tsl);
  phase_set_thermo_state (gas, &tsg);

  phase_set_properties (liq, MWs = (double[]){18.});
  phase_set_properties (gas, MWs = (double[]){MWG});

  /**
  The gas phase is an ideal gas: its density, thermal expansion coefficient and
  isothermal compressibility are set to the corresponding analytic functions. */

  tp2.rhov = gasprop_density_idealgas;
  tp2.betaT = gasprop_thermal_expansion;
  tp2.chiT = gasprop_isothermal_compressibility;

  /**
  We store the initial density of the ullage, which does not change in time
  because the volume and the mass of the ullage are constant. */

  rhog0 = tp2.rhov (&tsg);
}

/**
### Heat Source

The volumetric heat source is applied to the gas phase only. It plays the role
of the heat released by a chemical reaction, therefore we use the `chemistry`
event, which is called after the source terms are reset and before the
`divergence` event. In this way, the source term is included both in the
velocity divergence and in the temperature equation, together with the
compression work added by the phase change model. The source is weighted on the
volume fraction of the gas phase, consistently with the `theta` coefficient used
by the diffusion equation of the phase. */

event chemistry (i++) {
  scalar STexpG = gas->STexp;
  foreach()
    STexpG[] += q0*cm[]*(1. - f[]);
}

/**
## Post-Processing

The following lines of code are for post-processing purposes.
*/

/**
### Output Files

We write the thermodynamic pressure, the pressurization rate, and the gas phase
temperature, together with the corresponding analytic values and with the
maximum velocity magnitude, which must remain close to zero. */

event output_data (i++) {
  char name[80];
  sprintf (name, "OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, "w");

  P0exact = Pref + (gamma0 - 1.)*q0*t;
  Tg0exact = TG0 + gamma0*q0*t/(rhog0*cp2);

  double umax = 0.;
  foreach (reduction(max:umax))
    umax = max (umax, norm (u));

  /**
  The temperature of the ullage is probed in the middle of the gas region,
  far from the interface. */

  scalar TG = gas->T;
  double Tgas = interpolate (TG, 0.5*L0, 0.5*(1. + level0)*L0);

  double relerr = fabs (P0 - P0exact)/P0exact;

  fprintf (fp, "%g %g %g %g %g %g %g %g %g\n", t, P0, P0exact, relerr,
      dP0dt, (gamma0 - 1.)*q0, Tgas, Tg0exact, umax),
    fflush (fp);
}

/**
### Logger

We output the thermodynamic pressure of the tank in time (for testing). */

event logger (t += 0.02) {
  fprintf (stderr, "%d %.2f %.6g\n", maxlevel, t, P0);
}

/**
### Movie

We write the animation with the evolution of the gas phase temperature. */

event movie (t += 0.002; t <= 0.1) {
  if (maxlevel == 6) {
    scalar TG = gas->T;
    scalar Tplot[];
    foreach()
      Tplot[] = TG[]*(1. - f[]) + TL0*f[];

    clear();
    box();
    view (tx = -0.5, ty = -0.5);
    draw_vof ("f", lw = 2.);
    squares ("Tplot", min = TG0, max = TG0 + 15., linear = true);
    save ("movie.mp4");
  }
}

/**
## Results

The thermodynamic pressure follows the analytic constant-volume heating
solution.

~~~gnuplot Evolution of the thermodynamic pressure
reset
set grid
set key top left
set xlabel "t [s]"
set ylabel "P_0 [Pa]"
set size square

plot "OutputData-5" u 1:3 every 5 w p ps 0.9 pt 6 lc rgb "black" t "Analytic", \
     "OutputData-4" u 1:2 w l lw 2 lc 1 t "LEVEL 4", \
     "OutputData-5" u 1:2 w l lw 2 lc 2 t "LEVEL 5", \
     "OutputData-6" u 1:2 w l lw 2 dt 2 lc 3 t "LEVEL 6"
~~~

The temperature of the ullage follows the analytic solution, which includes the
compression work performed by the increasing thermodynamic pressure.

~~~gnuplot Evolution of the ullage temperature
reset
set grid
set key top left
set xlabel "t [s]"
set ylabel "T_g [K]"
set size square

plot "OutputData-5" u 1:8 every 5 w p ps 0.9 pt 6 lc rgb "black" t "Analytic", \
     "OutputData-4" u 1:7 w l lw 2 lc 1 t "LEVEL 4", \
     "OutputData-5" u 1:7 w l lw 2 lc 2 t "LEVEL 5", \
     "OutputData-6" u 1:7 w l lw 2 dt 2 lc 3 t "LEVEL 6"
~~~

Since the expansion of the ullage is exactly balanced by the compression due to
the pressure rise, the velocity field remains at the machine zero for the first
time steps, and it keeps a residual value which is due to the explicit coupling
between the pressurization rate and the compression work. If the volume
variation was redistributed uniformly over the domain, instead of being weighted
on the local compressibility, the incompressible liquid would be forced to
expand and the spurious velocity would be about one order of magnitude larger,
without any pressurization of the tank.

~~~gnuplot Spurious velocity field
reset
set grid
set key top left
set xlabel "t [s]"
set ylabel "max |u| [m/s]"
set logscale y
set format y "10^{%T}"
set size square
set key bottom right

plot "OutputData-4" u 1:9 w l lw 2 lc 1 t "LEVEL 4", \
     "OutputData-5" u 1:9 w l lw 2 lc 2 t "LEVEL 5", \
     "OutputData-6" u 1:9 w l lw 2 dt 2 lc 3 t "LEVEL 6"
~~~

The pressurization rate matches the analytic value $(\gamma-1)\dot{q}$. The
explicit coupling between the rate and the compression work produces a startup
transient, which decays geometrically by a factor $(\gamma-1)/\gamma$ per time
step.

~~~gnuplot Convergence of the pressurization rate
reset
set grid
set key bottom right
set xlabel "time step"
set ylabel "dP_0/dt / (γ-1)q̇ [-]"
set yrange [0:1.2]
set size square

plot "OutputData-4" u 0:($5/$6) w l lw 2 lc 1 t "LEVEL 4", \
     "OutputData-5" u 0:($5/$6) w l lw 2 lc 2 t "LEVEL 5", \
     "OutputData-6" u 0:($5/$6) w l lw 2 dt 2 lc 3 t "LEVEL 6", \
     1 w l lc rgb "black" dt 3 notitle
~~~
*/
