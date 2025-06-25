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
## Phase Change Setup

We define the number of gas and liquid species in the domain,
and we initialize all the properties necessary for the multicomponent phase
change model. The properties are set to null values because they are
overwritten by the variable properties formulation, which computes all the
physical properties as a function of the thermodynamic state of the mixture. */

/**
We set the initial temperature of the liquid and of the gas phase. */

/**
We solve the temperature field, with variable properties, and we reduce the
tolerance for the calculation of the variable properties. The interfacial
temperature is not computed from the jump conditon, it is just set to the gas
phase temperature value. */

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
#if JUMP
# include "navier-stokes/velocity-jump.h"
#else
# include "navier-stokes/low-mach.h"
#endif
#define FILTERED 1
#define P_ERR 0.1
#include "two-phase-varprop.h"
#include "opensmoke/properties.h"
#if VELOCITY_JUMP
# include "two-phase-clsvof.h"
#else
# include "two-phase.h"
#endif
#include "tension.h"
//#include "recoil.h"
#include "evaporation.h"
#include "defaultvars.h"
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
numerical simulation, and additional data for post-processing. */

int maxlevel, minlevel = 2;
double D0 = DIAMETER, effective_radius0 = 1., d_over_d02 = 1., tad = 0.;
double volumecorr = 0., volume0 = 0., d_over_d02_stop = 1., mLiq0;
bool restored = false;

int main (void) {
  /**
  We set the kinetics folder, which defines the species
  of the simulation, and it is used by OpenSMOKE++ for the
  calculation of the thermodynamic and transport properties. */

#if TWO_PHASE_VARPROP
  kinetics (TOSTRING(KINFOLDER), &NGS);
  kinetics_liquid (TOSTRING(KINFOLDER), &NLS);
  properties_liquid (TOSTRING(LIQFOLDER));
#else
  NGS = 2, NLS = 1;
#endif

  /**
  We set additional data for the simulation. */

  rho1 = 681.042; rho2 = 4.4165;
  mu1 = 0.00037446; mu2 = 3.50286e-5;
  Dmix1 = 0., Dmix2 = 3.90227e-6;
  lambda1 = 0.124069, lambda2 = 0.0554676;
  cp1 = 2244.92, cp2 = 1115.04;
  dhev = 364482;
  Pref = PRESSURE*101325.;

  TG0 = TEMPERATURE;
  TL0 = TEMPERATURE_DROPLET;

  nv = 2;
  pcm.emissivity = EMISSIVITY;

  /**
  We change the dimension of the domain as a function
  of the initial diameter of the droplet. */

  double RR = 7.986462e+01;
  L0 = 0.5*RR*D0;

  /**
  We change the surface tension coefficient. */

  f.sigma = 0.03;

  /**
  We run the simulation at different maximum
  levels of refinement. */

  maxlevel = MAXLEVEL;
  init_grid (1 << maxlevel);
  run();
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

/**
We initialize the volume fraction field and we compute the
initial radius of the droplet. We don't use the value D0
because for small errors of initialization the squared
diameter decay would not start from 1. */

scalar qspark[];

event init (i = 0) {
  if (!restore (file = "restart")) {
#if TREE
    //refine (circle (x, y, 4.*D0) > 0. && level < maxlevel);
#endif
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

    effective_radius0 = pow (3.*statsf(f).sum, 1./3.);
    volume0 = 4./3.*pi*pow (effective_radius0, 3.);

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
  On the top and right sides we set Dirichlet boundary conditions
  for the temperature and mass fraction fields. */

  scalar fuel  = gas->YList[0];
  scalar inert = gas->YList[1];
  scalar TG    = gas->T;

  fuel[top] = dirichlet (0.);
  fuel[right] = dirichlet (0.);

  inert[top] = dirichlet (1.);
  inert[right] = dirichlet (1.);

  TG[top] = dirichlet (TG0);
  TG[right] = dirichlet (TG0);
}

/**
We adapt the grid according to the mass fractions of the
mass fraction of n-heptane, the temperature, and the
velocity field. */

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
### Output Files

We write on a file the squared diameter decay and the
dimensionless time. */

event output_data (i++) {
  char name[80];
  sprintf (name, "OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, "w");

  double effective_radius = pow (3.*statsf(f).sum, 1./3.);
  double d_over_d02_old = d_over_d02;
  double tad_old = tad;

  d_over_d02 = sq (effective_radius / effective_radius0);
  d_over_d02_stop = sq (effective_radius / (0.5*D0));
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
  We compute and print additional useful average quantities. */

  double TIntAvg = avg_interface (gas_int->T, f, tol=0.1);
  double YIntAvg = avg_interface (gas_int->YList[0], f, tol=0.1);

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

  fprintf (fp, "%g %g %g %g %g %g %g %g %g %g\n", t, tad, effective_radius0,
      effective_radius, d_over_d02, mLiq, kv, TIntAvg, YIntAvg, TDropAvg);
}

/**
### Movie

We write the animation with the evolution of the
n-heptane mass fraction, the interface position
and the temperature field. */

#if MOVIE
event movie (t += MOVIE_EVERY) {
  clear();
  box();
  view (fov = 3, samples = 2);
  draw_vof ("f", lw = 1.5);
  squares ("T", min = TL0, max = statsf(T).max, linear = true);
  save ("temperature.mp4");

  clear();
  box();
  view (fov = 3, samples = 2);
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
  char name[80];
  sprintf (name, "snapshots-%g", t);
  dump (name);
}
#endif

/**
### Stopping Condition

We stop the simulation when the droplet is almost fully consumed. */

event stop (i++) {
  if (d_over_d02_stop <= MAX_DD02)
    return 1;
}

/**
We run the simulation for long time, which is not reached because
the stopping condition on the droplet diameter is reached first. */

event end (t = 50.);

