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

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

#ifndef TEMPERATURE
# define TEMPERATURE 773.
#endif

#ifndef TEMPERATURE_DROPLET
# define TEMPERATURE_DROPLET 300.
#endif

#ifndef PRESSURE
# define PRESSURE 1.
#endif

#ifndef DIAMETER
# define DIAMETER 1.e-3
#endif

#ifndef FUEL
# define FUEL NC7H16
#endif

#ifndef INERT
# define INERT N2
#endif

#ifndef OXIDIZER
# define OXIDIZER O2
#endif

#ifndef KINFOLDER
# define KINFOLDER evaporation/n-heptane-in-nitrogen
#endif

#ifndef LIQFOLDER
# define LIQFOLDER LiquidProperties
#endif

#ifndef RADIATION_INTERFACE
# define RADIATION_INTERFACE 0.93
#endif

#ifndef GRAVITY
# define GRAVITY -9.81
#endif

#ifndef FIBER
# define FIBER 0.1
#endif

/**
## Phase Change Setup

We define the number of gas and liquid species in the domain,
and we initialize all the properties necessary for the multicomponent phase
change model. The properties are set to null values because they are
overwritten by the variable properties formulation, which computes all the
physical properties as a function of the thermodynamic state of the mixture. */

#ifndef NGS
# define NGS 2
#endif

#ifndef NLS
# define NLS 1
#endif

#ifndef MASSFRAC_INERT
# define MASSFRAC_INERT 1.
#endif

#ifndef MASSFRAC_OXIDIZER
# define MASSFRAC_OXIDIZER 0.
#endif

#ifndef SPARK_START
# define SPARK_START 0.
#endif

#ifndef SPARK_TIME
# define SPARK_TIME 0.008
#endif

#ifndef SPARK_VALUE
# define SPARK_VALUE 1e7
#endif

#if COMBUSTION
char* gas_species[NGS];
char* liq_species[NLS];
#else
char* gas_species[NGS] = {TOSTRING(FUEL), TOSTRING(INERT)};
char* liq_species[NLS] = {TOSTRING(FUEL)};
#endif
char* inert_species[1] = {TOSTRING(INERT)};
double gas_start[NGS] = {0., 1.};
double liq_start[NLS] = {1.};
double inDmix1[NLS] = {0.};
double inDmix2[NGS] = {0.};
double inKeq[NLS] = {0.};

double lambda1 = 0.;
double lambda2 = 0.;
double dhev = 0.;
double cp1 = 0.;
double cp2 = 0.;

/**
We set the initial temperature of the liquid and of the gas phase. */

double TL0 = TEMPERATURE_DROPLET;
double TG0 = TEMPERATURE;

/**
We solve the temperature field, with variable properties, and we reduce the
tolerance for the calculation of the variable properties. The interfacial
temperature is not computed from the jump conditon, it is just set to the gas
phase temperature value. */

#define SOLVE_TEMPERATURE
#define USE_GSL 0
#define FSOLVE_ABSTOL 1.e-3
#define USE_ANTOINE_OPENSMOKE
#define FICK_CORRECTED
#define MOLAR_DIFFUSION
#define MASS_DIFFUSION_ENTHALPY
#define NO_ADVECTION_DIV 1
#define FILTERED

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
# include "navier-stokes/centered-evaporation.h"
# include "navier-stokes/centered-doubled.h"
#endif
#include "opensmoke-properties.h"
#include "pinning.h"
#include "two-phase.h"
#include "tension.h"
#include "gravity.h"
#include "recoil.h"
#include "evaporation.h"
#include "multicomponent-varprop.h"
#if USE_SPARK
# include "spark.h"
#endif
#if COMBUSTION
# include "chemistry.h"
# include "flame.h"
#endif
#include "view.h"

/**
### Boundary conditions

Outflow boundary conditions are set at the top and right
sides of the domain. */

#if JUMP
u1.n[top] = neumann (0.);
u1.t[top] = neumann (0.);
u2.n[top] = neumann (0.);
u2.t[top] = neumann (0.);
p[top] = dirichlet (0.);
ps[top] = dirichlet (0.);
pg[top] = dirichlet (0.);

u1.n[left] = neumann (0.);
u1.t[left] = neumann (0.);
u2.n[left] = neumann (0.);
u2.t[left] = neumann (0.);
p[left] = dirichlet (0.);
ps[left] = dirichlet (0.);
pg[left] = dirichlet (0.);

u1.n[right] = neumann (0.);
u1.t[right] = neumann (0.);
u2.n[right] = neumann (0.);
u2.t[right] = neumann (0.);
p[right] = dirichlet (0.);
ps[right] = dirichlet (0.);
pg[right] = dirichlet (0.);

u1.n[bottom] = dirichlet (0.);
u1.t[bottom] = dirichlet (0.);
u2.n[bottom] = dirichlet (0.);
u2.t[bottom] = dirichlet (0.);
p[bottom] = neumann (0.);
ps[bottom] = neumann (0.);
pg[bottom] = neumann (0.);
uf1.n[bottom] = 0.;
uf1.t[bottom] = 0.;
uf2.n[bottom] = 0.;
uf2.t[bottom] = 0.;
#else
u.n[top] = neumann (0.);
u.t[top] = neumann (0.);
p[top] = dirichlet (0.);
uext.n[top] = neumann (0.);
uext.t[top] = neumann (0.);
pext[top] = dirichlet (0.);

u.n[left] = neumann (0.);
u.t[left] = neumann (0.);
p[left] = dirichlet (0.);
uext.n[left] = neumann (0.);
uext.t[left] = neumann (0.);
pext[left] = dirichlet (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);
uext.n[right] = neumann (0.);
uext.t[right] = neumann (0.);
pext[right] = dirichlet (0.);

u.n[bottom] = dirichlet (0.);
u.t[bottom] = dirichlet (0.);
p[bottom] = neumann (0.);
uext.n[bottom] = dirichlet (0.);
uext.t[bottom] = dirichlet (0.);
pext[bottom] = neumann (0.);
uf.n[bottom] = 0.;
uf.t[bottom] = 0.;
ufext.n[bottom] = 0.;
ufext.t[bottom] = 0.;
#endif

/**
### Simulation Data

We declare the maximum and minimum levels of refinement,
the initial radius and diameter, and the radius from the
numerical simulation, and additional data for post-processing. */

int maxlevel, minlevel = 2;
double D0 = DIAMETER, effective_radius0;
double effective_radius = 0.5*DIAMETER, d_over_d02 = 1., tad = 0.;
double volumecorr = 0., trmin = 0., trmax = 0.;

scalar tr[];

int main (void) {
  /**
  We set the kinetics folder, which defines the species
  of the simulation, and it is used by OpenSMOKE++ for the
  calculation of the thermodynamic and transport properties. */

  kinfolder = TOSTRING(KINFOLDER);
  liqfolder = TOSTRING(LIQFOLDER);

  /**
  We set additional data for the simulation. */

  rho1 = 1.; rho2 = 1.;
  mu1 = 1.; mu2 = 1.;   // can't use 0 or [two-phase.h] won't create muf
  Pref = PRESSURE*101325.;

  /**
  We change the dimension of the domain as a function
  of the initial diameter of the droplet. */

  //double RR = 7.986462e+01;
  double RR = 160.;
  L0 = 0.5*RR*D0;

  G.x = GRAVITY;
  double df = FIBER*D0;
  X0 = -0.15*L0, Y0 = 0.5*df;

  /**
  We change the surface tension coefficient. */

  f.sigma = 0.03;
  f.tracers = {tr};

  /**
  We run the simulation at different maximum
  levels of refinement. */

  for (maxlevel = 11; maxlevel <= 11; maxlevel++) {

    pinning.ap = sqrt (sq (0.5*D0) - sq (Y0));
    pinning.ac = pinning.ap - 2.*L0/(1 << maxlevel);

    init_grid (1 << 9);
    run();
  }
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

/**
We initialize the volume fraction field and we compute the
initial radius of the droplet. We don't use the value D0
because for small errors of initialization the squared
diameter decay would not start from 1. */

scalar qspark[];

event init (i = 0) {
  refine (circle (x, y, 4.*D0) > 0. && level < maxlevel);
  fraction (f, circle (x, y, 0.5*D0));

  /**
  We compute initial variables useful for post-processing. */

  //volumecorr = 2.*pi*statsf(f).sum - (4./3.*pi*pow (0.5*D0, 3.));
  volumecorr = 0.;
  effective_radius0 = pow(3./4./pi*(2.*pi*statsf(f).sum - volumecorr), 1./3.);
  effective_radius = effective_radius0;

  foreach (reduction(+:mLiq0))
    mLiq0 += rho1v[]*f[]*dv();

  /**
  We set the molecular weights of the chemial species
  involved in the simulation (by default inMW=1). */

  foreach_elem (YGList, jj)
    inMW[jj] = OpenSMOKE_MW (jj);

  /**
  We introduce a scalar tracer to quantify the internal
  mixing in liquid phase. */

  foreach()
    tr[] = x*f[];

  trmin = statsf(tr).min;
  trmax = statsf(tr).max;

#ifdef RADIATION
  divq_rad = opensmoke_optically_thin;
#endif

#if USE_SPARK
  spark.T = qspark;
  //spark.position = (coord){0., 0.8*D0};
  //spark.position = (coord){0., 1.*D0};
  spark.position = (coord){0., 1.2*D0};
  spark.diameter = 0.2*D0;
  spark.time = SPARK_START;
  spark.duration = SPARK_TIME;
  spark.temperature = SPARK_VALUE;
#endif
}

#if COMBUSTION
event defaults (i = 0) {

  /**
  Initialize OpenSMOKE pointers. */

  char kinfolder_root[120];
  sprintf (kinfolder_root, "%s/kinetics/%s/kinetics",
      getenv ("OPENSMOKE_INTERFACE"), kinfolder);

  char liqfolder_root[120];
  sprintf (liqfolder_root, "%s/kinetics/LiquidProperties/%s",
      getenv ("OPENSMOKE_INTERFACE"), liqfolder);

  OpenSMOKE_Init();
  OpenSMOKE_ReadKinetics (kinfolder_root);
  OpenSMOKE_ReadLiquidKinetics (kinfolder_root);
  OpenSMOKE_ReadLiquidProperties (liqfolder_root);

  /**
  Reset species name using OpenSMOKE data. */

  for (int jj=0; jj<NGS; jj++) {
    gas_species[jj] = strdup (OpenSMOKE_NamesOfSpecies (jj));
    fprintf (stdout, "gas_species[%d] = %s\n", jj, gas_species[jj]);
  }

  for (int jj=0; jj<NLS; jj++) {
    const char * species = OpenSMOKE_NamesOfLiquidSpecies (jj);
    int len = strlen (species);
    char corrname[len+1];
    strcpy (corrname, species);
    corrname[3 <= len ? len-3 : 0] = '\0';
    liq_species[jj] = strdup (corrname);
    fprintf (stdout, "liq_species[%d] = %s\n", jj, liq_species[jj]);
  }

  /**
  Reset initial composition according to the new species. */

  for (int jj=0; jj<NGS; jj++)
    gas_start[jj] = 0.;

  gas_start[OpenSMOKE_IndexOfSpecies (TOSTRING(OXIDIZER))] = MASSFRAC_OXIDIZER;
  gas_start[OpenSMOKE_IndexOfSpecies (TOSTRING(INERT))] = MASSFRAC_INERT;

  /**
  Clean OpenSMOKE pointers to avoid conflicts with successive
  events. */

  OpenSMOKE_Clean();
}

event cleanup (t = end) {
  for (int jj=0; jj<NGS; jj++) {
    free (gas_species[jj]);
    gas_species[jj] = NULL;
  }

  for (int jj=0; jj<NLS; jj++) {
    free (liq_species[jj]);
    liq_species[jj] = NULL;
  }
}
#endif

/**
We use the same boundary conditions used by
[Pathak at al., 2018](#pathak2018steady). */

//event bcs (i = 0) {
//  scalar fuel  = YGList[OpenSMOKE_IndexOfSpecies (TOSTRING(FUEL))];
//  scalar inert = YGList[OpenSMOKE_IndexOfSpecies (TOSTRING(INERT))];
//
//  fuel[top] = dirichlet (0.);
//  fuel[left] = dirichlet (0.);
//  fuel[right] = dirichlet (0.);
//
//  inert[top] = dirichlet (MASSFRAC_INERT);
//  inert[left] = dirichlet (MASSFRAC_INERT);
//  inert[right] = dirichlet (MASSFRAC_INERT);
//
//#if COMBUSTION
//  scalar oxidizer = YGList[OpenSMOKE_IndexOfSpecies (TOSTRING(OXIDIZER))];
//
//  oxidizer[top] = dirichlet (MASSFRAC_OXIDIZER);
//  oxidizer[left] = dirichlet (MASSFRAC_OXIDIZER);
//  oxidizer[right] = dirichlet (MASSFRAC_OXIDIZER);
//#endif
//
//  TG[top] = dirichlet (TG0);
//  TG[left] = dirichlet (TG0);
//  TG[right] = dirichlet (TG0);
//}

/**
We adapt the grid according to the mass fractions of the
mass fraction of n-heptane, the temperature, and the
velocity field. */

#if TREE
event adapt (i++) {
  scalar fuel = YList[OpenSMOKE_IndexOfSpecies (TOSTRING(FUEL))];
  adapt_wavelet_leave_interface ({fuel,T,u.x,u.y}, {f},
      (double[]){1.e-1,1.e0,1.e-2,1.e-2}, maxlevel, minlevel, 1);
  unrefine (x >= (0.5*L0 + X0));
}
#endif

/**
## Post-Processing

The following lines of code are for post-processing purposes. */

/**
### Grashof Number

We compute the Grashof number to quantify the importance of the
natural convective fluxes. */

struct Grashof {
  double rhob, rhos;
  double r, g, nu;
  double value;
};

struct Grashof Gr;

event grashof (i++) {
#if COMBUSTION
  Gr.value = 1.;
#else
  if (i == 0) {
    ThermoState tsg;
    tsg.T = TG0;
    tsg.P = Pref;
    tsg.x = gas_start;

    Gr.rhob = tp2.rhov (&tsg);
    Gr.nu = tp2.muv (&tsg)/tp2.rhov (&tsg);
    effective_radius = effective_radius0;
  }
  Gr.r = effective_radius;
  Gr.g = fabs (GRAVITY);

  scalar YGIntFuel = YGIntList[0];
  scalar YGIntInert = YGIntList[1];

  double TIntAvg = avg_interface (TInt, f, tol=0.1);
  double YIntAvgFuel = avg_interface (YGIntFuel, f, tol=0.1);
  double YIntAvgInert = avg_interface (YGIntInert, f, tol=0.1);

  double YIntAvg[] = {YIntAvgFuel, YIntAvgInert};
  double XIntAvg[NGS];

  correctfrac (YIntAvg, NGS);
  mass2molefrac (XIntAvg, YIntAvg, inMW, NGS);

  ThermoState tsg;
  tsg.P = Pref;
  if (i == 0) {
    tsg.T = TL0;
    tsg.x = gas_start;
  }
  else {
    tsg.T = TIntAvg;
    tsg.x = XIntAvg;
  }
  Gr.rhos = tp2.rhov (&tsg);

  Gr.value = (Gr.rhos - Gr.rhob)*pow (Gr.r, 3.)*Gr.g/(Gr.rhob*sq(Gr.nu));
#endif
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
  foreach(reduction(+:mLiq))
    mLiq += rho1v[]*f[]*dv();

  /**
  We compute and print additional useful average quantities. */

  scalar YGIntFuel = YGIntList[0];
  double TIntAvg = avg_interface (TInt, f, tol=0.1);
  double YIntAvg = avg_interface (YGIntFuel, f, tol=0.1);

  int counter = 0;
  double TDropAvg = 0.;
  foreach(reduction(+:TDropAvg) reduction(+:counter)) {
    if (f[] > 1.-F_ERR) {
      counter++;
      TDropAvg += TL[];
    }
  }
  TDropAvg = (counter > 0.) ? TDropAvg/counter : 0.;

  fprintf (fp, "%g %g %g %g %g %g %g %g %g %g\n", t, tad, effective_radius,
      d_over_d02, mLiq/mLiq0, kv, Gr.value, TIntAvg, YIntAvg, TDropAvg);
  fflush (fp);
}

/**
### Movie

We write the animation with the evolution of the
n-heptane mass fraction, the interface position
and the temperature field. */

#if MOVIE
event movie (t += 0.001) {
# if COMBUSTION
  clear();
  box();
  view (tx = -0.025, fov = 5.5, samples = 2);
  draw_vof ("f", lw = 1.5);
  squares ("T", min = TL0, max = statsf(T).max, linear = true);
  isoline ("zmix - zsto", lw = 1.5, lc = {1.,1.,1.});
  save ("temperature.mp4");

  clear();
  box();
  view (tx = -0.025, fov = 5.5, samples = 2);
  draw_vof ("f", lw = 1.5);
  squares (TOSTRING(FUEL), min = 0., max = 1., linear = true);
  isoline ("zmix - zsto", lw = 1.5, lc = {1.,1.,1.});
  save ("fuel.mp4");
# else
  clear();
  box();
  view (tx = 0.025, fov = 3.5, samples = 2);
  draw_vof ("f", lw = 1.5);
  squares ("T", min = TL0, max = TG0, linear = true);
  save ("temperature.mp4");

  clear();
  box();
  view (tx = 0.025, fov = 3.5, samples = 2);
  draw_vof ("f", lw = 1.5);
  squares (TOSTRING(FUEL), min = 0., max = 1., linear = true);
  save ("fuel.mp4");
# endif
}
#endif

/**
### Snapshots

Output dump files for restore or post-processing. */

#if DUMP
event snapshots (t += 0.005) {
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
  if (d_over_d02 <= 0.05)
    return 1;
#if USE_SPARK
  if (t >= (SPARK_START + SPARK_TIME) && statsf(T).max < 1200.) {
    fprintf (ferr, "WARNING: Combustion did not start.\n");
    fflush (ferr);
    return 1;
  }
#endif
}

event end (t = 50.);

