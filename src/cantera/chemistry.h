#include "cantera/cantera.h"
#include "common-chemistry.h"

void cantera_ode_solver (ode_function batch,
    unsigned int NEQ, double dt, double * y0, void * args)
{
  size_t ns = thermo_nSpecies (thermo);
  bool isothermal = (ns == NEQ) ? true : false;

  UserDataODE * data = (UserDataODE *)args;
  double rho = data->rho;
  double cp = data->cp;
  double P = data->P;
  double T = isothermal ? data->T : y0[ns];
  double * sources = data->sources;

  double ymass[ns];
  for (int i = 0; i < ns; i++)
    ymass[i] = y0[i];

  thermo_setTemperature (thermo, T);
  thermo_setPressure (thermo, P);
  thermo_setMassFractions (thermo, ns, ymass, 1);

  int reactor = reactor_new ("IdealGasConstPressureReactor", soln, "batch");
  reactor_setEnergy (reactor, !isothermal);
  int net = reactornet_new();
  reactornet_addreactor (net, reactor);

  reactornet_advance (net, dt);

  // Recover results
  for (int i = 0; i < ns; i++) {
    ymass[i] = reactor_massFraction (reactor, i);
    y0[i] = ymass[i];
  }
  if (!isothermal) {
    T = reactor_temperature (reactor);
    y0[ns] = T;
  }

  // Computes derivatives (RHS) for the expansion term
  thermo_setTemperature (thermo, T);
  thermo_setPressure (thermo, P);
  thermo_setMassFractions (thermo, ns, ymass, 1);

  double wdot[ns], hm[ns];
  thermo_getPartialMolarEnthalpies (thermo, ns, hm);
  kin_getNetProductionRates (kin, ns, wdot);

  rho = thermo_density (thermo);
  cp = thermo_cp_mass (thermo);

  sources[ns] = 0.;
  for (int i = 0; i < ns; i++) {
    sources[i] = wdot[i]/rho;
    sources[ns] -= wdot[i]*hm[i];
  }
  sources[ns] /= (rho*cp);

  // Re-assign values to the UserDataODE, it may be useful for binning
  data->rho = rho;
  data->cp = cp;

  ct_clearReactors();
}

event defaults (i = 0) {
  stiff_ode_solver = cantera_ode_solver;
}
