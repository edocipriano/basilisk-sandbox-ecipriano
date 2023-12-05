/**
# Variable Properties

Simulations involving two-phase flows with variable material
properties can be performed using this module, which defines
structures and functions that help to setup such cases.
*/

#define VARPROP

/**
## Fields Allocations

We allocate fields required by this module. All properties
are intialized as constant fields, and they are initialized
just by the solver that needs them. If a thermal solver is
used, there is no need to initialize a non-constant
diffusivity of the chemical species, for example.
*/

(const) scalar rho1v = zeroc, rho2v = zeroc;
(const) scalar mu1v = zeroc, mu2v = zeroc;
(const) scalar lambda1v = zeroc, lambda2v = zeroc;
(const) scalar cp1v = zeroc, cp2v = zeroc;

/**
## Thermodynamic State

We define a structure with variables that define the thermodynamic
state of the mixture: temperature *T*, pressure *P*, composition *x*.
*/

typedef struct {
  double T, P;
  double * x;
} ThermoState;

/**
## Thermodynamic Properties

We define a struct that contains a bunch of function pointers
for each property of interest for the gas phase.
*/

typedef struct {
  // Mixture properties
  double (* rhov) (void *);
  double (* muv) (void *);
  double (* lambdav) (void *);
  double (* cpv) (void *);
  // Species properties
  double (* pvap) (void *, int);
  double (* dhev) (void *, int);
  double (* diff) (void *, int);
  double (* cps)  (void *, int);
} ThermoProps;

/**
## Useful functions

We define functions that are useful for variable properties
simulations.
*/

/**
### *check_termostate()*: check that the thermodynamic state is
reasonable. */

int check_thermostate (ThermoState * ts, int NS) {
  double sum = 0.;
  for (int jj=0; jj<NS; jj++)
    sum += ts->x[jj];

  int T_ok = (ts->T > 180. && ts->T < 4000.) ? true : false;
  int P_ok = (ts->P > 1e3 && ts->P < 1e7) ? true : false;
  int X_ok = (sum > 1.-1.e-3 && sum < 1.+1.e-3) ? true : false;

  return T_ok*P_ok*X_ok;
}

/**
### *print_thermostate()*: print the thermodynamic state of the mixture.
*/

void print_thermostate (ThermoState * ts, int NS, FILE * fp = stdout) {
  fprintf (fp, "Temperature = %g - Pressure = %g\n", ts->T, ts->P);
  for (int jj=0; jj<NS; jj++)
    fprintf (fp, "  Composition[%d] = %g\n", jj, ts->x[jj]);
  fprintf (fp, "\n");
}

/**
### *gasprop_thermal_expansion()*: Thermal expansion coefficient of an ideal gas
*/

double gasprop_thermal_expansion (ThermoState * ts) {
  return ts->T > 0. ? 1./ts->T : 0.;
}

/**
### *liqprop_thermal_expansion()*: Thermal expansion coefficient of a liquid
*/

double liqprop_thermal_expansion (ThermoProps * tp, ThermoState * ts) {
  double epsT = 1.e-3;
  double Ttop = ts->T + epsT, Tbot = ts->T - epsT;
  ThermoState tstop, tsbot;
  tstop.T = Ttop, tstop.P = ts->P, tstop.x = ts->x;
  tsbot.T = Tbot, tsbot.P = ts->P, tsbot.x = ts->x;
  double rhotop = tp->rhov (&tstop), rhobot = tp->rhov (&tsbot);
  double rhoval = tp->rhov (ts);
  return (rhoval > 0.) ? -1./rhoval*(rhotop - rhobot)/(2.*epsT) : 0.;
}

