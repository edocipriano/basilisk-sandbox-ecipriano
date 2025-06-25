/**
# Combustion of an axi droplet in microgravity or in any gravity conditions

This simulation file helps setting up numerical simulations of droplet
combustion in axial symmetry, in any gravity conditions, and considering
different kinetic schemes and chemical species.
*/

#include "axi.h"
#include "navier-stokes/low-mach.h"
#define FILTERED 1
#define P_ERR 0.1
#include "two-phase-varprop.h"
#include "opensmoke/properties.h"
#include "opensmoke/chemistry.h"
#include "pinning.h"
#if VELOCITY_JUMP
# include "two-phase-clsvof.h"
#else
# include "two-phase.h"
#endif
#include "tension.h"
#include "gravity.h"
#include "evaporation.h"
#include "spark.h"
#include "opensmoke/flame.h"
#include "defaultvars.h"
#include "view.h"

/**
## Simulation setups

We define an enum with different simulation setups: in microgravity and in
normal gravity conditions. We use this approach because the first configuration
place the droplet on the lower-left corner of the domain, in order to exploit
the axis of symmetry and to avoid non-physical lateral displacement. In normal
gravity, instead, we place the droplet along the edge of the bottom boundary,
and we include the droplet suspension strategy. This enum helps managing
differences among the two configurations.
*/

enum Setups { microgravity, gravity };
enum Setups setup = microgravity;

/**
## Input data

We declare a structure containing the input parameters for the simulations. We
are resolving not only the fluid dynamics but also the temperature field and
many chemical species (we don't know how many a-priori). Additionally, we need
to include a kinetic scheme, and the rules for calculating the material
properties in both the liquid and the gas phases. Since too many parameter have
to be provided, we gather them in a structure of input data and we read all the
parameters from an input file.
*/

struct InputData {
  char * kinetics_folder;
  char * liquid_properties_folder;
  int maxlevel;
  double Ytol;
  double Ttol;
  double Utol;
  double TL0;
  double TG0;
  double P0;
  double D0;
  double G;
  double sigma;
  double emissivity;
  double rfiber_to_rdrop;
  double rgas_to_rdrop;
  bool combustion;
  bool fick_corrected;
  bool molar_diffusion;
  bool mass_diffusion_enthalpy;
  bool divergence;
  char * liq_start;
  char * gas_start;
  double spark_diameter;
  double spark_start;
  double spark_time;
  double spark_value;
  double dump_every;
  double movie_every;
} inputdata = {
  TOSTRING(KINFOLDER),
  TOSTRING(LIQFOLDER),
  MAXLEVEL,
  YTOL,
  TTOL,
  UTOL,
  TEMPERATURE_DROPLET,
  TEMPERATURE,
  PRESSURE,
  DIAMETER,
  GRAVITY,
  SIGMA,
  EMISSIVITY,
  FIBER,
  ENVIRONMENT,
  COMBUSTION,
  true,         // fick_corrected
  true,         // molar_diffusion
  true,         // mass_diffusion_enthalpy
  true,         // divergence
  TOSTRING(LIQUID),
  TOSTRING(GAS),
  SPARK_DIAMETER,
  SPARK_START,
  SPARK_TIME,
  SPARK_VALUE,
  DUMP_EVERY,
  MOVIE_EVERY,
};

/**
## Parsing input file

We write a function which parse an input file with the simulation data, in case
one wants to avoid using compiler macros. We exploit libconfig. */

#if USE_LIBCONFIG

#include <libconfig.h>
#pragma autolink -lconfig

#define config_get_string(name, variable)                                     \
  config_lookup_string (&cfg, name, (const char **)&variable);

#define config_get_int(name, variable)                                        \
  config_lookup_int (&cfg, name, &variable);

#define config_get_bool(name, variable)                                       \
  config_lookup_bool (&cfg, name, (int *)&variable);

#define config_get_float(name, variable)                                      \
  config_lookup_float (&cfg, name, &variable);

#define config_get(type, name, variable)                                      \
  {                                                                           \
    int success = config_get_##type(name, variable);                          \
    if (success == CONFIG_FALSE)                                              \
      fprintf (stdout, "config: variable %s is set to default\n", #name),     \
        fflush (stdout);                                                      \
  }

void parse_input_file (char * filename, struct InputData * inputdata) {
  config_t cfg;
  config_init (&cfg);

  if (!config_read_file (&cfg, filename)) {
    fprintf (ferr, "error:%s:%d: %s:%d - %s\n",
        __FILE__, __LINE__,
        config_error_file(&cfg),
        config_error_line(&cfg),
        config_error_text(&cfg));
    config_destroy (&cfg);
    abort();
  }

  // Get strings
  config_get (string, "kinetics.kinfolder", inputdata->kinetics_folder);
  config_get (string, "kinetics.liqfolder", inputdata->liquid_properties_folder);
  config_get (string, "liquid.composition", inputdata->liq_start);
  config_get (string, "gas.composition", inputdata->gas_start);

  // Get integer values
  config_get (int, "domain.maxlevel", inputdata->maxlevel);

  // Get booleans
  config_get (bool, "physics.combustion", inputdata->combustion);
  config_get (bool, "physics.fick_corrected", inputdata->fick_corrected);
  config_get (bool, "physics.molar_diffusion", inputdata->molar_diffusion);
  config_get (bool, "physics.mass_diffusion_enthalpy", inputdata->mass_diffusion_enthalpy);
  config_get (bool, "physics.divergence", inputdata->divergence);

  // Get floating point values
  config_get (float, "domain.Ytol", inputdata->Ytol);
  config_get (float, "domain.Ttol", inputdata->Ttol);
  config_get (float, "domain.Utol", inputdata->Utol);
  config_get (float, "gas.temperature", inputdata->TG0);
  config_get (float, "liquid.temperature", inputdata->TL0);
  config_get (float, "gas.pressure", inputdata->P0);
  config_get (float, "domain.diameter", inputdata->D0);
  config_get (float, "domain.gravity", inputdata->G);
  config_get (float, "domain.fiber", inputdata->rfiber_to_rdrop);
  config_get (float, "domain.environment", inputdata->rgas_to_rdrop);
  config_get (float, "liquid.sigma", inputdata->sigma);
  config_get (float, "liquid.emissivity", inputdata->emissivity);
  config_get (float, "spark.diameter", inputdata->spark_diameter);
  config_get (float, "spark.start", inputdata->spark_start);
  config_get (float, "spark.time", inputdata->spark_time);
  config_get (float, "spark.value", inputdata->spark_value);
  config_get (float, "postprocessing.dump_every", inputdata->dump_every);
  config_get (float, "postprocessing.movie_every", inputdata->movie_every);
}

#endif

/**
## Useful variables

We declare other useful variabiles for this simulations.
*/

int maxlevel, minlevel = 2;
double D0, R0, R, M0, DD02 = 1., tad = 0.;
bool restored = false;
scalar qspark[];

/**
## Simulation setups

Here we have the main function and the intial event which set the properties
needed by the navier-stokes solver and by the phasechange model. */

int main (int argc, char ** argv) {

  /**
  We read input data from a file, if we prefer to avoid compiler macros. */

#if USE_LIBCONFIG
  char * inputfile = (argc > 1) ? argv[1] : "burningdroplet.input";
  parse_input_file (inputfile, &inputdata);
#endif

  /**
  We set the kinetics folders, for the gas and for the liquid kinetics. The
  liquid properties are initialized as well. */

  kinetics (inputdata.kinetics_folder, &NGS);
  kinetics_liquid (inputdata.kinetics_folder, &NLS);
  properties_liquid (inputdata.liquid_properties_folder);

  /**
  We specify the species names, in order to facilitate operations on the species
  when lage kinetic schemes are used (in this context large means more than 4
  species). */

  gas_species = new_species_names (NGS);
  liq_species = new_species_names_liquid (NLS);

  /**
  We set default simulation properties which will be overwritten by the variable
  material properties formulations. We can skip this part but WE NEED TO DEFINE
  NON-NULL VISCOSITY, because the [two-phase-generic.h](src/two-phase-generic.h)
  module creates the variable viscosity field only if mu is non-zero. */

  mu1 = mu2 = 1.;

  /**
  We set the initial thermodynamic state of the two-phase system (SI units). */

  Pref = inputdata.P0*101325.;
  TG0 = inputdata.TG0;
  TL0 = inputdata.TL0;

  /**
  We use a two-velocities model, and we set additional options for the phase
  change model `pcm`. */

  nv = 2;
  pcm.divergence = inputdata.divergence;
  pcm.chemistry = inputdata.combustion;
  pcm.fick_corrected = inputdata.fick_corrected;
  pcm.molar_diffusion = inputdata.molar_diffusion;
  pcm.mass_diffusion_enthalpy = inputdata.mass_diffusion_enthalpy;
  pcm.emissivity = inputdata.emissivity;

  /**
  We change the dimensions of the domain according with the simulation setup
  used. */

  double RR = inputdata.rgas_to_rdrop;
  D0 = inputdata.D0;
  L0 = 0.5*RR*D0;
  maxlevel = inputdata.maxlevel;

  /**
  Depending on the gravity value, we change the simulation setup. */

  G.x = inputdata.G;
  setup = G.x ? gravity : microgravity;

  if (setup == gravity) {
    double df = inputdata.rfiber_to_rdrop*D0;
    X0 = -0.5*L0, Y0 = 0.5*df;

    pinning.ap = sqrt (sq (0.5*D0) - sq (Y0));
    pinning.ac = pinning.ap - 2.*L0/(1 << maxlevel);
  }
  else {
    pinning.ap = 1e10;
    pinning.ac = 1e10;
    pinning_warnings = false;
  }

  /**
  We set the surface tension and additional tracers if needed. */

  f.sigma = inputdata.sigma;

  /**
  Setting up the grid and the pinning point if the setup requires droplet
  suspension. Finally, we run the simulation. */

  //init_grid (1 << min (maxlevel, 9));
  init_grid (1 << maxlevel);
  run();

  /**
  We de-allocate memory used for the vectors of species names. */

  free_species_names (NGS, gas_species), gas_species = NULL;
  free_species_names (NLS, liq_species), liq_species = NULL;
}

macro number circle (double x, double y, double R,
    double x0 = 0., double y0 = 0.)
{
  return (sq(R) - sq(x - x0) - sq(y - y0));
}

event init (i = 0) {
  if (!restore (file = "restart")) {
#if TREE
    refine (circle (x, y, 4.*D0) > 0. && level < maxlevel);
#endif
    fraction (f, circle (x, y, 0.5*D0));

    /**
    We compute initial variables, useful for post-processing. */

    if (setup == gravity)
      R = R0 = pow(3./4./pi*(2.*pi*statsf(f).sum), 1./3.);
    else
      R = R0 = pow (3.*statsf(f).sum, 1./3.);

    scalar rhol = liq->rho;
    M0 = 0.;
    foreach (reduction(+:M0))
      M0 += rhol[]*f[]*dv();
  }

  /**
  We set the initial thermodynamic state of the phases, and we overwrite some of
  their properties. The function `phase_set_thermo_state` sets the initial
  conditions of temperature, pressure, and mass fraction fields only if they are
  uniform. Therefore, if `restore` was called, the function should not overwrite
  the fields. However, it is still important to call it in order to store the
  initial thermo state */

  ThermoState tsl, tsg;
  tsl.T = TL0, tsl.P = Pref, tsl.x = NULL;
  tsg.T = TG0, tsg.P = Pref, tsg.x = NULL;

  phase_set_thermo_state (liq, &tsl);
  phase_set_thermo_state (gas, &tsg);

  phase_set_composition_from_string (liq, inputdata.liq_start, sep = "_");
  phase_set_composition_from_string (gas, inputdata.gas_start, sep = "_");

  /**
  The only property that we need to set, and that remains constant throughout
  the simulation, is the vector with the molecular weight of each species. */

  double MWLs[NLS], MWGs[NGS];
  molecular_weights (NGS, MWGs);
  foreach_species_in (liq)
    MWLs[i] = MWGs[LSI[i]];

  phase_set_properties (liq, MWs = MWLs);
  phase_set_properties (gas, MWs = MWGs);

#if TWO_PHASE_VARPROP
  antoine = &opensmoke_antoine;
#endif

  /**
  We set the spark properties. */

#if SPARK
  spark.T = qspark;
  spark.position = (setup == gravity) ? (coord){0., 1.2*D0}
                                      : (coord){0.75*D0, 0.75*D0};
  spark.diameter = inputdata.spark_diameter*D0;
  spark.time = inputdata.spark_start;
  spark.duration = inputdata.spark_time;
  spark.temperature = inputdata.spark_value;
  spark.phase = gas;
#endif
}

/**
## Boundary conditions

We set the boundary conditions in the `defaults` event because, in
[navier-stokes/centered-new.h](navier-stokes/centered-new.h) the boundary
conditions for the additional velocity fields are copied from u and p.
*/

event defaults (i = 0) {
  u.n[top] = neumann (0.);
  u.t[top] = neumann (0.);
  p[top] = dirichlet (0.);

  u.n[right] = neumann (0.);
  u.t[right] = neumann (0.);
  p[right] = dirichlet (0.);

  if (setup == gravity) {
    u.n[left] = neumann (0.);
    u.t[left] = neumann (0.);
    p[left] = dirichlet (0.);

    u.n[bottom] = dirichlet (0.);
    u.t[bottom] = dirichlet (0.);
    p[bottom] = neumann (0.);
    uf.n[bottom] = 0.;
    uf.t[bottom] = 0.;
  }
}

/**
## Mesh adaptation

We adapt the grid according to the fuel mass fraction, the temperature, and the
velocity field. */

#if TREE
event adapt (i++) {
  double Ytol = inputdata.Ytol;
  double Ttol = inputdata.Ttol;
  double Utol = inputdata.Utol;

  adapt_wavelet_leave_interface ({Y,T,u.x,u.y}, {f},
      (double[]){Ytol,Ttol,Utol,Utol}, maxlevel, minlevel, 1);

  if (setup == gravity)
    unrefine (x >= (0.45*L0));
}
#endif

/**
## Output Files

We write on a file interesting average properties in time, such as the droplet
diameter, surface, mass, interface temperature and mass fractions. */

event output_data (i += 50) {
  if (setup == gravity)
    R = pow(3./4./pi*(2.*pi*statsf(f).sum), 1./3.);
  else
    R = pow (3.*statsf(f).sum, 1./3.);

  DD02 = sq (R/R0);
  tad = t/sq (D0*1e3);

  double TIntAvg = avg_interface (gas_int->T, f, tol = 0.1);

  int counter = 0.;
  double TDropAvg = 0.;
  scalar TL = liq->T;
  foreach(reduction(+:TDropAvg) reduction(+:counter)) {
    if (f[] > 1.-F_ERR) {
      counter++;
      TDropAvg += TL[];
    }
  }
  TDropAvg = (counter > 0.) ? TDropAvg/counter : 0.;

  char name[80];
  sprintf (name, "OutputData-%d", maxlevel);

  static FILE * fp = fopen (name, "w");

  fprintf (fp, "%g %g %g %g %g %g %g", t, tad, R0, R, DD02,
      TIntAvg, TDropAvg);

  foreach_species_in (liq_int) {
    scalar YGInt = gas_int->YList[LSI[i]];
    double YGIntAvg = avg_interface (YGInt, f, tol = 0.1);
    fprintf (fp, " %g", YGIntAvg);
  }
  fprintf (fp, "\n"), fflush (fp);
}

/**
## Snapshots

We write periodic sanpshots of the simulation in order to restart the simulation
if needed, or for post-processing. */

event snapshots (t += inputdata.dump_every) {
  if (i > 1) {
    char name[80];
    sprintf (name, "snapshots-%g", t);
    dump (name);
  }
}

/**
## Movie

We write the animation with the evolution of the fuel mass fraction, the
temperature field, the interface position and the flame front. */

#if 0
event movie (t += inputdata.movie_every) {
  if (setup == microgravity) {
    clear();
    view (fov = 2);
    draw_vof ("f", lw = 1.5);
    squares ("Y", min = 0., max = 1., linear = true);
    isoline ("zmix - zsto", lw = 1.5, lc = {1.,1.,1.});
    mirror ({0,1}) {
      draw_vof ("f", lw = 1.5);
      squares ("YGO2", min = 0., max = 0.21, linear = true);
      isoline ("zmix - zsto", lw = 1.5, lc = {1.,1.,1.});
    }
    mirror ({1,0}) {
      draw_vof ("f", lw = 1.5);
      squares ("T", min = TL0, max = 2100., linear = true);
      isoline ("zmix - zsto", lw = 1.5, lc = {1.,1.,1.});
      mirror ({0,1}) {
        cells();
        draw_vof ("f", lw = 1.5);
        isoline ("zmix - zsto", lw = 1.5, lc = {1.,1.,1.});
      }
    }
    save ("movie.mp4");
  }
}
#endif

/**
## Stopping conditions

We assume that in 50 s the droplet is fully consumed. Alternatively, we stop the
simulation when the droplet is almost fully consumed, i.e. when the square
diameter decay is small: $(D/D_0)^2 < 0.005$. */

event stop (DD02 < MAX_DD02) {
  return 1;
}

#if BASILISK_SANDBOX
event end (t = 0.02);
#else
event end (t = 50);
#endif

