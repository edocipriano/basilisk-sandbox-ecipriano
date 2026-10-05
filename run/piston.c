/**
# Bagnold Piston

Active pressurization of a tank due to the inflow of liquid. A channel of length
$l$ and height $h$ is filled with a liquid layer on the left and with an ideal
gas ullage on the right. Liquid is injected through the left boundary with a
mass flux $\dot{m}_{in}$, while the other boundaries are adiabatic slip walls.
The liquid is incompressible, therefore it pushes the interface to the right
like a piston, with velocity $u_\Gamma = \dot{m}_{in}/\rho_l$, compressing the
ullage.

The thermal conductivity is null in both phases, therefore the compression of
the gas is adiabatic and reversible, and the pressure and the temperature of the
ullage follow the isentropic relations:
$$
  P_0(t) = \left(\dfrac{l - x_{\Gamma,0}}{l - x_{\Gamma,0} - u_\Gamma t}
    \right)^\gamma P_0(0)
$$
$$
  T_g(t) = \left(\dfrac{P_0(0)}{P_0(t)}\right)^{\frac{1-\gamma}{\gamma}} T_0
$$
This test case verifies the pressurization model implemented in
[low-mach.h](../src/navier-stokes/low-mach.h) when the volume of the domain is
fixed but some boundaries are not walls: the inflow of liquid must be balanced
by the compression of the ullage, i.e. by an increase of the thermodynamic
pressure $P_0$.

![Evolution of the interface and the gas phase
temperature](piston/movie.mp4)(width="100%")
*/

/**
## Simulation Setup

The geometry is a uniform rectangular channel, which is obtained using a
multigrid with an aspect ratio $l/h = 10$. We use the low Mach number
Navier--Stokes solver together with the variable properties module. The phase
change is switched off using the fixed flux model with a null vaporization
rate. */

#include "grid/multigrid.h"
#include "navier-stokes/low-mach.h"
#define P_ERR 1.e-6
#include "two-phase-varprop.h"
#include "basilisk-properties.h"
#include "two-phase.h"
#include "phasechange.h"
#include "fixedflux.h"
#include "view.h"

/**
### Model Data

The data of the benchmark: a hydrogen-like ideal gas, the inflow mass flux, the
length and height of the channel. The initial position of the interface
$x_{\Gamma,0}$ is not given by the benchmark: the value $0.3\,l$ reproduces the
reference temperature evolution (from 20.4 K to about 24.3 K in 100 s). */

#define RGAS    5049.     // Specific gas constant [J/kg/K]
#define GAMMA   1.78      // Heat capacity ratio of the gas [-]
#define PIN     103.e3    // Initial thermodynamic pressure [Pa]
#define TIN     20.4      // Initial and inflow temperature [K]
#define MDOTIN  0.1       // Inflow mass flux [kg/m2/s]
#define RHORAT  70.       // Liquid to gas density ratio [-]
#define LENGTH  1.        // Length of the channel [m]
#define HEIGHT  0.1       // Height of the channel [m]
#define XGAMMA0 0.3       // Initial position of the interface [m]

int maxlevel;
double rhog0, uin;

/**
The density of the gas phase follows the ideal gas law, which makes the ullage
compressible and sensitive to the variation of the thermodynamic pressure. */

double gasprop_density_idealgas (void * p) {
  ThermoState * ts = p;
  return ts->P/(RGAS*ts->T);
}

/**
### Boundary Conditions

Liquid is injected from the left boundary. The other boundaries keep the
default symmetry conditions, which correspond to adiabatic slip walls. The
boundary conditions of the volume fraction and of the liquid temperature are set
in the `init` event, because the phase fields are not yet allocated here. */

u.n[left] = dirichlet (uin);
u.t[left] = dirichlet (0.);
f[left] = dirichlet (1.);
T[left] = dirichlet (TIN);

int main (void) {

  /**
  We use a single chemical species in each phase. */

  NGS = 1, NLS = 1;

  /**
  We set the material properties. The thermal conductivity and the diffusivity
  are null, in order to obtain an adiabatic compression of the ullage. The gas
  phase density is overwritten by the ideal gas law, and the heat capacity is
  obtained from the heat capacity ratio. The viscosities and the liquid heat
  capacity are not given by the benchmark, and we use values for hydrogen at
  20 K: they do not affect the solution, because the flow is a plug flow in a
  channel with slip walls, and the liquid temperature is uniform. */

  P0 = PIN;
  TG0 = TIN, TL0 = TIN;

  rhog0 = P0/(RGAS*TIN);
  rho1 = RHORAT*rhog0, rho2 = rhog0;
  mu1 = 1.3e-5, mu2 = 1.1e-6;
  Dmix1 = 0., Dmix2 = 0.;
  lambda1 = 0., lambda2 = 0.;
  cp1 = 9.7e3, cp2 = GAMMA*RGAS/(GAMMA - 1.);
  dhev = 0.;

  /**
  The liquid enters with velocity $u_\Gamma = \dot{m}_{in}/\rho_l$, which is
  also the velocity of the interface. */

  uin = MDOTIN/rho1;

  /**
  The volume of the domain is fixed: the net inflow is converted into a
  variation of the thermodynamic pressure. The volumetric flow rate entering
  from the left boundary (per unit depth) is not balanced by any outlet, and it
  is passed to the closed-system projection as a negative net outflow. */

  closed = true;
  closed_flux = -uin*HEIGHT;

  /**
  The composition of the two phases does not change, and the phase change is
  switched off setting a null vaporization rate. */

  pcm.isomassfrac = true;
  pcm.divergence = true;
  mEvapVal = 0.;

  /**
  The velocity divergence in the ullage is of the order of $u_\Gamma/l \approx
  10^{-3}$ s$^{-1}$, therefore the tolerance of the Poisson solver must be
  reduced with respect to the default value. */

  TOLERANCE = 1.e-6;
  DT = 0.1;

  L0 = LENGTH;
  int nx = round (LENGTH/HEIGHT);
  dimensions (nx = nx, ny = 1);

  for (maxlevel = 3; maxlevel <= 5; maxlevel++) {
    init_grid (nx*(1 << maxlevel));
    run();
  }
}

/**
We initialize a flat interface at $x_{\Gamma,0}$, with the liquid on the left
and the gas on the right. */

event init (i = 0) {
  fraction (f, XGAMMA0 - x);

  /**
  We initialize the velocity with the exact solution: the liquid moves with the
  inflow velocity, while the gas is compressed uniformly, and its velocity
  decreases linearly to zero at the right wall. Starting from rest, instead, the
  face velocity predicted at the inflow boundary differs from the imposed one
  during the first time steps, and the Poisson problem for `pf` is not
  compatible with `closed_flux`. */

  foreach()
    u.x[] = (x < XGAMMA0) ? uin : uin*(LENGTH - x)/(LENGTH - XGAMMA0);

  ThermoState tsl, tsg;
  tsl.T = TL0, tsl.P = P0, tsl.x = (double[]){1.};
  tsg.T = TG0, tsg.P = P0, tsg.x = (double[]){1.};

  phase_set_thermo_state (liq, &tsl);
  phase_set_thermo_state (gas, &tsg);

  phase_set_properties (liq, MWs = (double[]){R_GAS*1.e3/RGAS});
  phase_set_properties (gas, MWs = (double[]){R_GAS*1.e3/RGAS});

  /**
  The gas phase is an ideal gas: its density, thermal expansion coefficient and
  isothermal compressibility are set to the corresponding analytic functions. */

  tp2.rhov = gasprop_density_idealgas;
  tp2.betaT = gasprop_thermal_expansion;
  tp2.chiT = gasprop_isothermal_compressibility;
}

/**
### Thermodynamic Pressure

The thermodynamic pressure is integrated in time using the pressurization rate
computed by the projection step, which includes the volumetric flow rate
entering from the left boundary through `closed_flux`. */

event pressurization (i++) {
  P0 += dt*dP0dt;
}

/**
## Post-Processing

The following lines of code are for post-processing purposes.
*/

/**
### Output Files

We write the thermodynamic pressure and the volume-averaged temperature of the
ullage, together with the analytic values. We also write the position of the
interface, computed from the volume of liquid, and the mass of gas, which must
be conserved. */

event output_data (i++) {
  char name[80];
  sprintf (name, "OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, "w");

  /**
  The phase temperature is stored as a tracer, i.e. multiplied by the volume
  fraction of the phase. The mass of gas is computed from the current
  thermodynamic pressure and temperature, $\rho_g = P_0/(R_g T_g)$: the density
  field of the gas phase is updated at the beginning of the time step, while the
  volume fraction is already advected, and their product would underestimate
  the mass by $\Delta t\,d\ln\rho_g/dt$. */

  scalar TG = gas->T;
  double vliq = 0., vgas = 0., tgas = 0., mgas = 0.;
  foreach (reduction(+:vliq) reduction(+:vgas) reduction(+:tgas)
      reduction(+:mgas)) {
    vliq += f[]*dv();
    vgas += (1. - f[])*dv();
    tgas += TG[]*dv();
    if (1. - f[] > F_ERR && TG[] > 0.)
      mgas += P0/RGAS*sq(1. - f[])/TG[]*dv();
  }
  double Tgas = tgas/vgas;
  double xgamma = vliq/HEIGHT;

  static double mgas0 = 0.;
  if (i == 0)
    mgas0 = mgas;

  double P0exact = pow ((LENGTH - XGAMMA0)/(LENGTH - XGAMMA0 - uin*t),
      GAMMA)*PIN;
  double Tgexact = pow (PIN/P0exact, (1. - GAMMA)/GAMMA)*TIN;

  fprintf (fp, "%g %g %g %g %g %g %g %g\n", t, P0, P0exact, Tgas, Tgexact,
      xgamma, XGAMMA0 + uin*t, (mgas - mgas0)/mgas0),
    fflush (fp);
}

/**
### Logger

We output the thermodynamic pressure and the ullage temperature in time (for
testing). */

event logger (t += 10.) {
  double Tgas = interpolate (gas->T, 0.9*LENGTH, 0.5*HEIGHT);
  fprintf (stderr, "%d %g %.6g %.6g\n", maxlevel, t, P0, Tgas);
}

/**
### Movie

We write the animation with the evolution of the temperature field. */

event movie (t += 1; t <= 100) {
  if (maxlevel == 5) {
    clear();
    view (fov = 4.2, tx = -0.5, ty = -0.5*HEIGHT, width = 4000, height = 800);
    box (notics = true);
    draw_vof ("f", lw = 4.);
    squares ("T", min = TIN, max = 24.5, linear = true);
    save ("movie.mp4");
  }
}

/**
## Results

The thermodynamic pressure follows the isentropic compression of the ullage.

~~~gnuplot Evolution of the thermodynamic pressure
reset
set grid
set key top left
set xlabel "t [s]"
set ylabel "P_0 [kPa]"
set size square

plot "OutputData-4" u 1:($3/1e3) every 50 w p ps 0.9 pt 6 lc rgb "black" \
       t "Analytic", \
     "OutputData-3" u 1:($2/1e3) w l lw 2 lc 1 t "LEVEL 3", \
     "OutputData-4" u 1:($2/1e3) w l lw 2 lc 2 t "LEVEL 4", \
     "OutputData-5" u 1:($2/1e3) w l lw 2 dt 2 lc 3 t "LEVEL 5"
~~~

The volume-averaged temperature of the ullage is compared with the reference
solution of the benchmark.

~~~gnuplot Evolution of the ullage temperature
reset
set grid
set key top left
set xlabel "t [s]"
set ylabel "T_g [K]"
set size square

plot "OutputData-4" u 1:5 every 50 w p ps 0.9 pt 6 lc rgb "black" \
       t "Analytic", \
     "OutputData-3" u 1:4 w l lw 2 lc 1 t "LEVEL 3", \
     "OutputData-4" u 1:4 w l lw 2 lc 2 t "LEVEL 4", \
     "OutputData-5" u 1:4 w l lw 2 dt 2 lc 3 t "LEVEL 5"
~~~

The interface moves with the velocity of the inflow, since the liquid is
incompressible.

~~~gnuplot Position of the interface
reset
set grid
set key top left
set xlabel "t [s]"
set ylabel "x_Γ [m]"
set size square

plot "OutputData-4" u 1:7 every 50 w p ps 0.9 pt 6 lc rgb "black" \
       t "Analytic", \
     "OutputData-3" u 1:6 w l lw 2 lc 1 t "LEVEL 3", \
     "OutputData-4" u 1:6 w l lw 2 lc 2 t "LEVEL 4", \
     "OutputData-5" u 1:6 w l lw 2 dt 2 lc 3 t "LEVEL 5"
~~~

The mass of gas should be conserved: the ullage is compressed, not emptied. The
residual error is of order $\Delta$, and it grows only while the interface cell
is mostly filled with liquid ($f > 0.5$): the split VOF advection adds the term
$\Delta t\,c_c\nabla\cdot\mathbf{u}$ to the volume fraction, with $c_c = (f >
0.5)$, therefore the compression of the gas in the interface cell removes
liquid instead of gas.

~~~gnuplot Relative variation of the gas mass
reset
set grid
set key top left
set xlabel "t [s]"
set ylabel "(m_g - m_g^0)/m_g^0 [-]"
set size square

plot "OutputData-3" u 1:8 w l lw 2 lc 1 t "LEVEL 3", \
     "OutputData-4" u 1:8 w l lw 2 lc 2 t "LEVEL 4", \
     "OutputData-5" u 1:8 w l lw 2 dt 2 lc 3 t "LEVEL 5"
~~~
*/
