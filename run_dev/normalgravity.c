/**
# Evaporation of a Pure Droplet in Buoyancy-Driven Flows

This module reports the simulation setup for a generic pure droplet
at different gravity conditions. This configuration has been used in
several experimental works, by suspending the droplet in a static
configuration by means of a solid fiber. This module mimick the
experimental setup by combining a droplet suspension strategy with
the interface-resolved phase change model and the variable properties
formulation, which allows buoyancy-driven flows to be resolved directly.
These phenomena promote non-radial gas-phase convective fluxes, and
liquid phase internal recirculation, which break the spherical symmetry
of the evaporation process in microgravity conditions.
The simulation setup can be used with any gravity value, not only in
normal gravity. This feature allows small residual gravity values to be
simulated, making this simulation setup useful to drive the experimental
investigation also for droplets in microgravity.

This simulation includes variable thermodynamic and transport properties,
and the interface radiation contribution. We neglect the heat transfer
from the suspending solid fiber, which may affect experimental data if
the fiber is large. The simulation is performed by simulating a droplet
at the edge of a square domain, exploting the axial-symmetry, which is
acceptable for small Reynolds number.
*/

/**
## Default Simulation Data

The simulation setup is independent from the liquid under investigation,
and we may want to perform the same simulation studying the sensitivity
of the numerical results to the operative conditions (temperature,
pressure, ecc...) and to the liquid fuel or the ambient conditions.
To avoid code duplication, we set a bunch of default compiler variables,
which are overwritten in the `Makefile` to create different cases with
different operative conditions.

The default properties are:

* initial ambient temperature: `TEMPERATURE = 773 K`
* initial droplet temperature: `TEMPERATURE_DROPLET = 300 K`
* constant thermodynamic pressure: `PRESSURE = 10 atm`
* initial droplet diameter: `DIAMETER = 1 mm`
* emissivity of the liquid fuel: `RADIATION_INTERFACE = 0.93`
* name of the liquid fuel: `FUEL = n-heptane`
* name of the inert/ambient species: `INERT = nitrogen`
* path of the kinetics folder: `KINFOLDER = evaporation/n-heptane-in-nitrogen`
* path of the liquid properties folder: `LIQFOLDER = LiquidProperties`
* gravitational acceleration: `GRAVITY = -9.81 ` $m/s^2$
* ratio between solid fiber diameter and droplet diameter: `FIBER = 0.1`

these properties describe the evaporation of a n-heptane droplet in
nitrogen in normal gravity conditions, including the interface radiation.
These properties can be easily changed by overwriting the compilation
variables, for example by setting `-DTEMPERATURE=473 -DPRESSURE=1` to
change the ambient temperature and pressure. The interface radiation is
suppressed by setting a null value of emissivity: `-DRADIATION_INTERFACE=0`.

Eventually, one may consider the possibility of compiling a code which
reads the properties directly from an input file (using libconfig for example).

The following example considers a 1 mm droplet, without interface radiation,
in a 773 K environment. The simulation runs until 0.5 s. During this transient,
we observe the heating of the liquid phase and the formation of a downward
wake (P = 10 atm). The simulation is stopped before the complete consumption
of the droplet in order to reduce the simulation time on the basilisk server
as much as possible.

![Temperature field](normalgravity/temperature.mp4)
*/

/**
We solve both mass fractions and temperature fields, also in the interface
jump condition. The GSL library is used for root finding operations, whose
tolerance is controlled by the variable `FSOLVE_ABSTOL`. The thermodynamic
equlibrium is computed from Antoine's law, and we solve the momentum equation
in non-conservative form by setting the variable `NO_ADVECTION_DIV` to 1.

Additional compiler variables are used to activate terms in the governing
equations. In particular: `FICK_CORRECTED` is used to force the diffusive
fluxes to close to zero, even if the diffusivity of every chemical species
is different; `MOLAR_DIFFUSION` is used to correct Fick's law considering
that the diffusivity values are mole fractions-based; `MASS_DIFFUSION_ENTHALPY`
includes the species diffusion contribution in the temperature eqation. All
of them should be used.
*/

/**
## Simulation Setup

In this simulation, the Navier-Stokes equations can be solved both using the
centered solver with divergence source term, or using the centered solver with
velocity jump. The latter should be preferred since it is able to limit oscillations in the velocity field. The [two-phase.h](/src/two-phase.h)
solver is extended to variable physical properties by including a policy
for the calculation of such properties. In this case, we use correlations implemented in OpenSMOKE++ libraries.
The solution of the jump conditions and of the temperature and mass fraction
fields is performed by the multicomponent phase change model. We include the
recoil pressure contribution, and the non-reduced gravity contribution. Using
[gravity.h](/sandbox/ecipriano/src/gravity.h) simplifies the boundary condition
for pressure, analogously to the [reduced.h](/src/reduced.h) approach but
considering variable liquid and gas phase densities.
*/

#include "axi.h"
#if JUMP
# include "navier-stokes/velocity-jump.h"
#else
# include "navier-stokes/low-mach.h"
#endif
#define FILTERED 1
#define P_ERR 0.1
#include "two-phase-varprop.h"
#include "opensmoke/properties.h"
#include "pinning.h"
#if VELOCITY_JUMP
# include "two-phase-clsvof.h"
#else
# include "two-phase.h"
#endif
#include "tension.h"
#include "gravity.h"
//#include "recoil.h"
#include "evaporation.h"
#include "defaultvars.h"
#include "view.h"

/**
### Boundary conditions

We initialize the droplet at the lower edge of the domain, and we
exploit axial symmetry. We set no-slip boundary conditions on the side
in contact with the droplet (left), and outflow boundary
conditions on the other sides of the domain. */

u.n[top] = neumann (0.);
u.t[top] = neumann (0.);
p[top] = dirichlet (0.);

u.n[left] = neumann (0.);
u.t[left] = neumann (0.);
p[left] = dirichlet (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);

u.n[bottom] = dirichlet (0.);
u.t[bottom] = dirichlet (0.);
p[bottom] = neumann (0.);
uf.n[bottom] = 0.;
uf.t[bottom] = 0.;

/**
### Simulation Data

We declare the maximum and minimum levels of refinement,
the initial droplet diameter, and additional data for
post-processing. We also transport a scalar tracer `tr`
which can be useful to visualize liquid internal recirculation. */

int maxlevel, minlevel = 2;
double D0 = DIAMETER, effective_radius0, mLiq0;
double effective_radius = 0.5*DIAMETER, d_over_d02 = 1., tad = 0.;
double volumecorr = 0., trmin = 0., trmax = 0.;
bool restored = false;

scalar tr[];

int main (void) {

  /**
  We set the kinetics folder, which defines the species
  of the simulation, and it is used by OpenSMOKE++ for the
  calculation of the thermodynamic and transport properties.
  The path is relative to `OpenSMOKEppInterface/kinetics`.
  We do the same for the folder gathering the properties of
  the liquid species. */

#if TWO_PHASE_VARPROP
  kinetics (TOSTRING(KINFOLDER), &NGS);
  kinetics_liquid (TOSTRING(KINFOLDER), &NLS);
  properties_liquid (TOSTRING(LIQFOLDER));
#else
  NGS = 2, NLS = 1;
#endif

  /**
  We set additional simulation properties. Be careful, these are "dummy"
  properties which are overwritten by the variable-properties formulation
  already at the first iteration. To start the simulation without any
  problem of division by zero, the viscosity value should be non-null. */

  rho1 = 681.042; rho2 = 9.75415;
  mu1 = 0.00037446; mu2 = 2.02391e-05;
  lambda1 = 0.124069, lambda2 = 0.0295641;
  cp1 = 2244.92, cp2 = 1041.52;
  dhev = 364482;

  /**
  We set the thermodynamic pressure in SI units (Pa). */

  Pref = PRESSURE*101325.;
  TG0 = TEMPERATURE, TL0 = TEMPERATURE_DROPLET;

  nv = 2;

  /**
  We change the dimensions of the domain as a function of the initial
  diameter of the droplet. The domain is large with respect to the
  droplet diameter, in order to avoid the influence of the domain
  boundaries on the droplet evaporation dynamics. */

  double RR = ENVIRONMENT;
  L0 = 0.5*RR*D0;

  /**
  We set the gravity contribution, and we shift the origin along the `y`
  coordinate, according to the radius of the solid fiber. */

  G.x = GRAVITY;
  double df = FIBER*D0;
  X0 = -0.5*L0, Y0 = 0.5*df;

  /**
  We use a constant surface tension value, and we transport `tr` as
  a VOF tracer. */

  f.sigma = 0.03;
  f.tracers = {tr};

  /**
  We run the simulation at different maximum levels of refinement. */

  maxlevel = MAXLEVEL;

  /**
  The pinning point `ap` is adjusted according to the value of `Y0`.
  The pinning center `ac` is also adjusted depending on the grid. */

  pinning.ap = sqrt (sq (0.5*D0) - sq (Y0));
  pinning.ac = pinning.ap - 2.*L0/(1 << maxlevel);

  init_grid (1 << maxlevel);
  run();
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

/**
We initialize the volume fraction field according to `D0`,
after refining the region around the droplet to the maximum
level of refinement. */

event init (i = 0) {
  if (!restore (file = "restart")) {
    //refine (circle (x, y, 4.*D0) > 0. && level < maxlevel);
    fraction (f, circle (x, y, 0.5*D0));
#if VELOCITY_JUMP
    foreach()
      d[] = circle (x, y, 0.5*D0);
    vertex scalar phi[];
    foreach_vertex()
      phi[] = (d[] + d[-1] + d[0,-1] + d[-1,-1])/4.;
    fractions (phi, f);
#endif

    /**
    We compute initial variables useful for post-processing. */

    //volumecorr = 2.*pi*statsf(f).sum - (4./3.*pi*pow (0.5*D0, 3.));
    volumecorr = 0.;
    effective_radius0 = pow(3./4./pi*(2.*pi*statsf(f).sum - volumecorr), 1./3.);
    effective_radius = effective_radius0;

    scalar rhol = liq->rho;
    foreach (reduction(+:mLiq0))
      mLiq0 += rhol[]*f[]*dv();
  }
  else
    restored = true;

  /**
  We set the molecular weights of the chemial species
  involved in the simulation (by default inMW=1). */

  /**
  We set the properties of the system. */

  ThermoState tsl, tsg;
  tsl.T = TL0, tsl.P = Pref, tsl.x = (double[]){1.};
  tsg.T = TG0, tsg.P = Pref, tsg.x = (double[]){0.,1.};

  phase_set_thermo_state (liq, &tsl);
  phase_set_thermo_state (gas, &tsg);

  phase_set_properties (liq, MWs = (double[]){100.2});
  phase_set_properties (gas, MWs = (double[]){100.2,29.});

  scalar lfuel = liq->YList[0];
  lfuel.antoine = antoine_heptane;

#if TWO_PHASE_VARPROP
  foreach_species_in (gas)
    gas->MWs[i] = OpenSMOKE_MW (i);
  foreach_species_in (liq)
    liq->MWs[i] = OpenSMOKE_MW (LSI[i]);

  antoine = &opensmoke_antoine;
#endif

  /**
  On the edge in contact with the liquid droplet we set Neumann
  boundary conditions for the scalar fields, while we set Dirichlet
  boundary conditions for the other sides. */

  scalar fuel  = gas->YList[0];
  scalar inert = gas->YList[1];
  scalar TG = gas->T;

  fuel[top] = dirichlet (0.);
  fuel[left] = dirichlet (0.);
  fuel[right] = dirichlet (0.);

  inert[top] = dirichlet (1.);
  inert[left] = dirichlet (1.);
  inert[right] = dirichlet (1.);

  TG[top] = dirichlet (TG0);
  TG[left] = dirichlet (TG0);
  TG[right] = dirichlet (TG0);

  /**
  The scalar tracer is initialized with the value of the
  horizontal coordinate. */

  foreach()
    tr[] = x*f[];

  trmin = statsf(tr).min;
  trmax = statsf(tr).max;
}

/**
We adapt the grid according to the fuel mass fraction,
the temperature, and the velocity fields. */

#if TREE
event adapt (i++) {
  adapt_wavelet_leave_interface ({Y,T,u.x,u.y}, {f},
      (double[]){YTOL,TTOL,UTOL,UTOL}, maxlevel, minlevel, 1);
}
#endif

/**
## Post-Processing

The following lines of code are for post-processing purposes. */

/**
### Grashof Number

We compute the Grashof number to quantify the intensity of the
natural convective fluxes. */

struct Grashof {
  double rhob, rhos;
  double r, g, nu;
  double value;
};

struct Grashof Gr;

event grashof (i++) {
  if (i == 0) {
    ThermoState tsg;
    tsg.T = TG0;
    tsg.P = Pref;
    tsg.x = gas->ts0->x;

    Gr.rhob = tp2.rhov (&tsg);
    Gr.nu = tp2.muv (&tsg)/tp2.rhov (&tsg);
    effective_radius = effective_radius0;
  }
  Gr.r = effective_radius;
  Gr.g = fabs (GRAVITY);

  scalar YGIntFuel = gas_int->YList[0];
  scalar YGIntInert = gas_int->YList[1];

  double TIntAvg = avg_interface (gas_int->T, f, tol=0.1);
  double YIntAvgFuel = avg_interface (YGIntFuel, f, tol=0.1);
  double YIntAvgInert = avg_interface (YGIntInert, f, tol=0.1);

  double YIntAvg[] = {YIntAvgFuel, YIntAvgInert};
  double XIntAvg[NGS];

  correctfrac (YIntAvg, NGS);
  mass2molefrac (XIntAvg, YIntAvg, gas->MWs, NGS);

  ThermoState tsg;
  tsg.P = Pref;
  if (i == 0) {
    tsg.T = TL0;
    tsg.x = gas->ts0->x;
  }
  else {
    tsg.T = TIntAvg;
    tsg.x = XIntAvg;
  }
  Gr.rhos = tp2.rhov (&tsg);

  Gr.value = (Gr.rhos - Gr.rhob)*pow (Gr.r, 3.)*Gr.g/(Gr.rhob*sq(Gr.nu));
}

/**
### Output Files

We write on a file the squared diameter decay and the
dimensionless time. */

event output_data (i += 50) {
  char name[80];
  sprintf (name, "OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, "w");

  effective_radius = pow(3./4./pi*(2.*pi*statsf(f).sum - volumecorr), 1./3.);
  double d_over_d02_old = d_over_d02;
  double tad_old = tad;

  d_over_d02 = sq (effective_radius / effective_radius0);
  tad = t/sq(D0*1e3);

  /**
  The vaporization rate is computed according to the formula
  in Liu & Avedisian, 2011, pag. 777 bottom. */

  double kv = 0.;
  if (i > 1)
    kv = fabs ((d_over_d02 - d_over_d02_old)/(tad - tad_old));

  double mLiq = 0.;
  scalar rhol = liq->rho;
  foreach(reduction(+:mLiq))
    mLiq += rhol[]*f[]*dv();

  /**
  We compute and print useful average quantities such as the
  average interface temperature and mass fractions, and the
  average droplet temperature. */

  scalar YGIntFuel = gas_int->YList[0];
  double TIntAvg = avg_interface (gas_int->T, f, tol=0.1);
  double YIntAvg = avg_interface (YGIntFuel, f, tol=0.1);

  int counter = 0;
  double TDropAvg = 0.;
  scalar TL = liq->T;
  foreach(reduction(+:TDropAvg) reduction(+:counter)) {
    if (f[] > 1.-F_ERR) {
      counter++;
      TDropAvg += TL[];
    }
  }
  TDropAvg = (counter > 0.) ? TDropAvg/counter : 0.;

  fprintf (fp, "%g %g %g %g %g %g %g %g %g %g %g\n", t, tad, effective_radius0,
      effective_radius, d_over_d02, mLiq/mLiq0, kv, Gr.value, TIntAvg, YIntAvg,
      TDropAvg);
  fflush (fp);
}

/**
### Movie

We write the animation with the evolution of the
fuel mass fraction, the interface position
and the temperature field. */

#if MOVIE
event movie (t += MOVIE_EVERY) {
  clear();
  box();
  view (tx = 0.025, fov = 3.5, samples = 2,
      theta=0., phi=0., psi=-pi/2.);
  draw_vof ("f", lw = 1.5);
  squares ("T", min = TL0, max = TG0, linear = true);
  save ("temperature.mp4");

  clear();
  box();
  view (tx = 0.025, fov = 3.5, samples = 2,
      theta=0., phi=0., psi=-pi/2.);
  draw_vof ("f", lw = 1.5);
  squares ("Y", min = 0., max = 1., linear = true);
  save ("fuel.mp4");
}
#endif

/**
### Snapshots

Output dump files for restore or post-processing. */

#if DUMP
event snapshots (t += DUMP_EVERY) {
  if (i > 1) {
    char name[80];
    sprintf (name, "snapshots-%g", t);
    dump (name);
  }
}
#endif

/**
### Stopping Condition

We stop the simulation when the droplet is almost fully consumed. */

event stop (i++) {
  if (d_over_d02 <= MAX_DD02)
    return 1;
}

/**
We run the simulation for long time, which is not reached because
the stopping condition on the droplet diameter is reached first. */

#if BASILISK_SERVER
event end (t = 0.5);
#else
event end (t = 50.);
#endif

