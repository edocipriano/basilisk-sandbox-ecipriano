/**
# Cantera Properties

We compute the material properties of a mixture using the Cantera library.
*/

#include "cantera/cantera.h"
#include "variable-properties.h"

#define VARIABLE_PROPERTIES 1

/**
## Properties Functions

Functions for the update of the density, given the thermodynamic
state.
*/

/**
### *cantera_gasprop_density()*: gas phase density according to the ideal gas low
*/

double cantera_gasprop_density (void * p) {
  ThermoState * ts = (ThermoState *)p;
  thermo_setTemperature (thermo, ts->T);
  thermo_setPressure (thermo, ts->P);
  size_t ns = thermo_nSpecies (thermo);
  thermo_setMoleFractions (thermo, ns, ts->x, 1);
  //return thermo_density (thermo);
  double MWmix = thermo_meanMolecularWeight (thermo);
  return ts->P*MWmix / (R_GAS*1e3*ts->T);
}

/**
### *cantera_gasprop_viscosity()*: gas phase dynamic viscosity
*/

double cantera_gasprop_viscosity (void * p) {
  ThermoState * ts = (ThermoState *)p;
  thermo_setTemperature (thermo, ts->T);
  thermo_setPressure (thermo, ts->P);
  size_t ns = thermo_nSpecies (thermo);
  thermo_setMoleFractions (thermo, ns, ts->x, 1);
  return trans_viscosity (tran);
}

/**
### *cantera_gasprop_thermalconductivity()*: gas phase thermal conductivity
*/

double cantera_gasprop_thermalconductivity (void * p) {
  ThermoState * ts = (ThermoState *)p;
  thermo_setTemperature (thermo, ts->T);
  thermo_setPressure (thermo, ts->P);
  size_t ns = thermo_nSpecies (thermo);
  thermo_setMoleFractions (thermo, ns, ts->x, 1);
  return trans_thermalConductivity (thermo);
}

/**
### *cantera_gasprop_heatcapacity()*: gas phase specific heat capacity
*/

double cantera_gasprop_heatcapacity (void * p) {
  ThermoState * ts = (ThermoState *)p;
  thermo_setTemperature (thermo, ts->T);
  thermo_setPressure (thermo, ts->P);
  size_t ns = thermo_nSpecies (thermo);
  thermo_setMoleFractions (thermo, ns, ts->x, 1);
  return thermo_cp_mass (thermo);
}

/**
### *cantera_gasprop_heatcapacity_species()*: gas phase species heat capacity
*/

void cantera_gasprop_heatcapacity_species (void * p, double * r) {
  ThermoState * ts = (ThermoState *)p;
  size_t ns = thermo_nSpecies (thermo);
  thermo_setTemperature (thermo, ts->T);
  thermo_setPressure (thermo, ts->P);
  thermo_setMoleFractions (thermo, ns, ts->x, 1);

  thermo_getPartialMolarCp (thermo, ns, r);

  double MW[ns];
  thermo_getMolecularWeights (thermo, ns, MW);

  for (int i = 0; i < ns; i++)
    r[i] /= MW[i];
}

/**
### *cantera_gasprop_diff()*: diffusion coefficient of a species in gas phase
*/

void cantera_gasprop_diff (void * p, double * r) {
  ThermoState * ts = (ThermoState *)p;
  size_t ns = thermo_nSpecies (thermo);
  thermo_setTemperature (thermo, ts->T);
  thermo_setPressure (thermo, ts->P);
  double xm[ns];
  for (int i = 0; i < ns; i++)
    xm[i] = max (ts->x[i], 1e-10);
  thermo_setMoleFractions (thermo, ns, xm, 1);
  trans_getMixDiffCoeffs (tran, ns, r);
}

#if 0
/**
### *cantera_liqprop_density_addvol()*: liquid phase mixture density with the additive volume method
*/

double cantera_liqprop_density_addvol (void * p) {
  extern double rho1;
  return rho1;
}

/**
### *cantera_liqprop_viscosity()*: liquid phase mixture dynamic viscosity
*/

double cantera_liqprop_viscosity (void * p) {
  extern double mu1;
  return mu1;
}

/**
### *cantera_liqprop_thermalconductivity()*: liquid phase mixture thermal conductivity
*/

double cantera_liqprop_thermalconductivity (void * p) {
  extern double lambda1;
  return lambda1;
}

/**
### *cantera_liqprop_heatcapacity()*: liquid phase mixture specific heat capacity
*/

double cantera_liqprop_heatcapacity (void * p) {
  extern double cp1;
  return cp1;
}

/**
### *cantera_liqprop_heatcapacity_species()*: liquid phase species heat capacity
*/

void cantera_liqprop_heatcapacity_species (void * p, double * r) {
  extern double cp1;
  for (int i = 0; i < thermo_nSpecies (thermo_liq); i++)
    r[i] = cp1;
}

/**
### *cantera_liqprop_dhev()*: vapor pressure of the chemical species
*/

void cantera_liqprop_dhev (void * p, double * r) {
  extern double dhev;
  for (int i = 0; i < thermo_nSpecies (thermo_liq); i++)
    r[i] = dhev;
}

/**
### *cantera_liqprop_sigma()*: surface tension of the chemical species
*/

void cantera_liqprop_sigma (void * p, double * r) {
  for (int i = 0; i < thermo_nSpecies (thermo_liq); i++)
      r[i] = 0.03;
}

/**
### *cantera_liqprop_diff()*: diffusion coefficient of a species in liquid phase
*/

void cantera_liqprop_diff (void * p, double * r) {
  extern double Dmix1;
  for (int i = 0; i < thermo_nSpecies (thermo_liq); i++)
    r[i] = Dmix1;
}
#endif

/**
### *cantera_antoine()*: implementation of the antoine function using opensmoke
*/

double cantera_antoine (double T, double P, int i) {
  return 0.;
}

/**
### *cantera_gasprop_thermal_expansion()*: gas thermal expansion coefficient
*/

double cantera_gasprop_thermal_expansion (const void * p, void * s) {
  ThermoState * ts = (ThermoState *)s;
  return (ts->T > 0.) ? 1./ts->T : 0.;
}

/**
### *cantera_gasprop_species_expansion()*: gas species expansion coefficient
*/

void cantera_gasprop_species_expansion (const void * p, void * s, double * r) {
  ThermoState * ts = (ThermoState *)s;
  thermo_setTemperature (thermo, ts->T);
  thermo_setPressure (thermo, ts->P);
  size_t ns = thermo_nSpecies (thermo);
  thermo_setMoleFractions (thermo, ns, ts->x, 1);
  double MWmix = thermo_meanMolecularWeight (thermo);

  double MW[ns];
  thermo_getMolecularWeights (thermo, ns, MW);

  for (int i = 0; i < ns; i++)
    r[i] = MWmix / MW[i];
}

/**
### *cantera_liqprop_thermal_expansion()*: liq thermal expansion coefficient
*/

double cantera_liqprop_thermal_expansion (const void * p, void * s) {
  ThermoProps * tp = (ThermoProps *)p;
  ThermoState * ts = (ThermoState *)s;

  double epsT = 1.e-3;
  double Ttop = ts->T + epsT, Tbot = ts->T - epsT;
  ThermoState tstop, tsbot;
  tstop.T = Ttop, tstop.P = ts->P, tstop.x = ts->x;
  tsbot.T = Tbot, tsbot.P = ts->P, tsbot.x = ts->x;
  double rhotop = tp->rhov (&tstop), rhobot = tp->rhov (&tsbot);
  double rhoval = tp->rhov (ts);
  return (rhoval > 0.) ? -1./rhoval*(rhotop - rhobot)/(2.*epsT) : 0.;
}

/**
### *cantera_liqprop_species_expansion()*: liq species expansion coefficient
*/

void cantera_liqprop_species_expansion (const void * p, void * s, double * r)
{
  for (unsigned int i = 0; i < thermo_nSpecies (thermo_liq); i++)
    r[i] = 0.;
}

/**
## Thermodynamic Properties

We create the instance of two structures with the
thermodynamic properties, *tp1* for the liquid phase
and *tp2* for the gas phase. The same nomenclature is used
for the thermodynamic states.
*/

ThermoProps tp1, tp2;

/**
## Initialization

We set the thermodynamic properties function pointers
to the specific opensmoke functions declared above.
*/

event defaults (i = 0) {

  /**
  We set the thermodynamic properties functions to the
  correct opensmoke functions that compute material
  properties. */

#if 1
  tp2.rhov    = cantera_gasprop_density;
  tp2.muv     = cantera_gasprop_viscosity;
  tp2.lambdav = cantera_gasprop_thermalconductivity;
  tp2.cpv     = cantera_gasprop_heatcapacity;
  tp2.diff    = cantera_gasprop_diff;
  tp2.cps     = cantera_gasprop_heatcapacity_species;
  tp2.betaT   = cantera_gasprop_thermal_expansion;
  tp2.betaY   = cantera_gasprop_species_expansion;

  //tp2.rhov    = cantera_gasprop_density;
  //tp2.muv     = cantera_gasprop_viscosity;
  //tp2.lambdav = cantera_gasprop_thermalconductivity;
  //tp2.cpv     = cantera_gasprop_heatcapacity;
  //tp2.diff    = cantera_gasprop_diff;
  //tp2.cps     = cantera_gasprop_heatcapacity_species;
  //tp2.betaT   = cantera_gasprop_thermal_expansion;
  //tp2.betaY   = cantera_gasprop_species_expansion;
#else
  tp2.rhov    = NULL;
  tp2.muv     = NULL;
  tp2.lambdav = NULL;
  tp2.cpv     = NULL;
  tp2.diff    = NULL;
  tp2.cps     = NULL;
  tp2.betaT   = NULL;
  tp2.betaY   = NULL;
#endif

#if 0
  tp1.rhov    = cantera_liqprop_density_addvol;
  tp1.muv     = cantera_liqprop_viscosity;
  tp1.lambdav = cantera_liqprop_thermalconductivity;
  tp1.cpv     = cantera_liqprop_heatcapacity;
  tp1.dhev    = cantera_liqprop_dhev;
  tp1.diff    = cantera_liqprop_diff;
  tp1.cps     = cantera_liqprop_heatcapacity_species;
  tp1.sigmas  = NULL;
  tp1.betaT   = cantera_liqprop_thermal_expansion;
  tp1.betaY   = cantera_liqprop_species_expansion;
#else
  tp1.rhov    = NULL;
  tp1.muv     = NULL;
  tp1.lambdav = NULL;
  tp1.cpv     = NULL;
  tp1.dhev    = NULL;
  tp1.diff    = NULL;
  tp1.cps     = NULL;
  tp1.sigmas  = NULL;
  tp1.betaT   = NULL;
  tp1.betaY   = NULL;
#endif
}

