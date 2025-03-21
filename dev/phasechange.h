#include "common-evaporation.h"
#include "phase.h"

scalar mEvapTot[];
extern scalar * mEvapList;

int NLS = 0, NGS = 0;
double lambda1 = 1., lambda2 = 1., dhev = 1., cp1 = 1., cp2 = 1.;
double TG0 = 300., TL0 = 300., TIntVal = 300., Pref = 101325.;
double Dmix1 = 0., Dmix2 = 1., YIntVal = 0., YG0 = 0., YL0 = 1.;

Phase * liq, * gas;

scalar fu[], fu1[], * fulist = {fu,fu1};  // fixme: fu1 should be dynamic
scalar T[], Y[];

#include "regression.h"

enum advection_policy {
  NO_ADVECTION, ADVECTION_VELOCITY, SOURCE_TERM, PLANE_SHIFTING
};

enum shifting_policy {
  NO_SHIFTING, SHIFT_TO_LIQUID, SHIFT_TO_GAS
};

enum velocity_policy {
  WITHOUT_EXPANSION, WITH_EXPANSION
};

struct PhaseChangeModel {
  enum advection_policy advection;
  enum shifting_policy shifting;
  enum velocity_policy velocity;
  bool boiling;
  bool byrhogas;
  bool consistent;
  bool isothermal;
  bool isomassfrac;
} pcm = {
  ADVECTION_VELOCITY,
  SHIFT_TO_LIQUID,
  WITHOUT_EXPANSION,
  false,
  false,
  false,
  false,
  false,
};

void intexp_explicit (scalar intexp, scalar f, scalar mEvapTot) {
  scalar rhol = liq->rho, rhog = gas->rho;
  foreach_interfacial_plic (f, F_ERR)
    intexp[] = (rhol[] > 0. && rhog[] > 0.) ?
      mEvapTot[]*(1./rhog[] - 1./rhol[])*dirac : 0.;
}

event defaults (i = 0) {
  liq = new_phase ("L", NLS, false);
  gas = new_phase ("G", NGS, true);

  liq->isothermal = pcm.isothermal;
  gas->isothermal = pcm.isothermal;
  liq->isomassfrac = pcm.isomassfrac;
  gas->isomassfrac = pcm.isomassfrac;

  ThermoState tsl, tsg;
  tsl.T = TL0, tsl.P = Pref, tsl.x = (double[]){1};
  tsg.T = TG0, tsg.P = Pref, tsg.x = (double[]){1};

  phase_set_thermo_state (liq, &tsl);
  phase_set_thermo_state (gas, &tsg);

  double * Dl = malloc (liq->n*sizeof (double));
  double * Dg = malloc (gas->n*sizeof (double));
  double * dhevsl = malloc (liq->n*sizeof (double));
  double * dhevsg = malloc (gas->n*sizeof (double));
  double * cpsl = malloc (liq->n*sizeof (double));
  double * cpsg = malloc (gas->n*sizeof (double));

  foreach_species_in (liq) {
    Dl[i] = Dmix1;
    dhevsl[i] = dhev;
    cpsl[i] = cp1;
  }

  foreach_species_in (gas) {
    Dg[i] = Dmix2;
    dhevsg[i] = dhev;
    cpsg[i] = cp2;
  }

  phase_set_properties (liq,
      rho = rho1, mu = mu1,
      lambda = lambda1, cp = cp1,
      dhev = dhev, dhevs = dhevsl,
      D = Dl, cps = cpsl);

  phase_set_properties (gas,
      rho = rho2, mu = mu2,
      lambda = lambda2, cp = cp2,
      dhev = dhev, dhevs = dhevsg,
      D = Dg, cps = cpsg);

  free (Dl);
  free (Dg);
  free (dhevsl);
  free (dhevsg);
  free (cpsl);
  free (cpsg);

  phase_scalars_to_tracers (liq, f);
  phase_scalars_to_tracers (gas, f);

  phase_set_tracers (liq);
  phase_set_tracers (gas);

  if (nv == 1) {
    scalar fug = fulist[0];
    fug.tracers = list_concat (fug.tracers, liq->tracers);
    fug.tracers = list_concat (fug.tracers, gas->tracers);
  }
  else if (nv == 2) {
    scalar ful = fulist[1], fug = fulist[0];
    ful.tracers = list_concat (ful.tracers, liq->tracers);
    fug.tracers = list_concat (fug.tracers, gas->tracers);
  }

  for (scalar c in fulist) {
    c.refine = c.prolongation = fraction_refine;
    c.dirty = true;
    scalar * tracers = c.tracers;
    for (scalar t in tracers) {
      t.restriction = restriction_volume_average;
      t.refine = t.prolongation = vof_concentration_refine;
      t.dirty = true;
      t.c = c;
    }
  }
}

event init (i = 0) {
#if VELOCITY_JUMP
  _boiling = pcm.boiling;
#endif
}

event cleanup (t = end) {
  for (scalar fu in fulist)
    free (fu.tracers), fu.tracers = NULL;
  free (fulist), fulist = NULL;

  delete_phase (liq);
  delete_phase (gas);
}

event reset_sources (i++) {
  phase_reset_sources (liq);
  phase_reset_sources (gas);

  foreach() {
    for (scalar mEvap in mEvapList)
      mEvap[] = 0.;
    for (scalar fu in fulist)
      fu[] = f[];
    for (scalar intexp in intexplist)
      intexp[] = 0.;
    for (scalar drhodt in drhodtlist)
      drhodt[] = 0.;
  }
}

event end_init (i = 0);

event reset_sources (i++);

event chemistry (i++);

event phasechange (i++) {
  foreach() {
    mEvapTot[] = 0.;
    for (scalar mEvap in mEvapList)
      mEvapTot[] += mEvap[];
  }
#if VELOCITY_JUMP
  scalar rhol = liq->rho, rhog = gas->rho;
  foreach()
    jump[] = (rhol[] > 0. && rhog[] > 0.) ?
      -mEvapTot[]*(1./rhol[] - 1./rhog[]) : 0.;

  scalar fext[], dext[];
  foreach() {
    fext[] = f[];
    dext[] = -d[];
  }

  extern scalar d;
  constant_extrapolation (jump, dext, 0.5, 10, c=fext, nl=0,
            nointerface=true, inverse=false, tol=1e-4);
  constant_extrapolation (jump, dext, 0.5, 10, c=fext, nl=0,
            nointerface=true, inverse=true, tol=1e-4);
#else
  scalar intexp = (pcm.boiling && nv > 1) ? intexplist[1] : intexplist[0];
  intexp_explicit (intexp, f, mEvapTot);

  switch (pcm.shifting) {
    case NO_SHIFTING: break;
    case SHIFT_TO_LIQUID: shift_field (intexp, f, 1); break;
    case SHIFT_TO_GAS: shift_field (intexp, f, 0); break;
  }
#endif
}

face vector ufsave[];

event vof (i++) {
  // Transport due to the phase change
  scalar rhoh = pcm.byrhogas ? gas->rho : liq->rho;
  switch (pcm.advection) {
    case NO_ADVECTION:
      break;
    case ADVECTION_VELOCITY:
      vof_advection_phasechange (f, mEvapTot, rhoh, i);
      break;
    case SOURCE_TERM:
      vof_expl_sources (f, mEvapTot, rhoh, dt);
      break;
    case PLANE_SHIFTING:
      vof_plane_shifting (f, mEvapTot, rhoh, dt);
  }

  // Transport tracers with uf by default (i.e. nv == 1)
  scalar fug = fulist[0];
  vof_advection ({fug}, i);

  // Transport tracers with liquid velocity
  if (nv == 2) {
    scalar ful = fulist[1];
    face vector ufl = uflist[1];
    foreach_face() {
      ufsave.x[] = uf.x[];
      uf.x[] = ufl.x[];
    }
    vof_advection ({ful}, i);
  }
}

event vof_sources (i++) {
  if (nv == 2)
    foreach_face()
      uf.x[] = ufsave.x[];

  foreach() {
    f[] = clamp (f[], 0., 1.);
    f[] = (f[] < F_ERR) ? 0. : (f[] > 1.-F_ERR) ? 1. : f[];
  }
}

event tracer_advection (i++);

event tracer_diffusion (i++) {
  if (nv == 1) {
    phase_tracers_to_scalars (liq, fu, tol = F_ERR);
    phase_tracers_to_scalars (gas, fu, tol = F_ERR);
  }
  else if (nv == 2) {
    scalar ful = fulist[1], fug = fulist[0];
    phase_tracers_to_scalars (liq, ful, tol = F_ERR);
    phase_tracers_to_scalars (gas, fug, tol = F_ERR);
  }

#if TWO_PHASE_VARPROP
  phase_diffusion (gas, f, varcoeff = true);
  phase_diffusion (liq, f, varcoeff = true);
#else
  phase_diffusion (gas, f, varcoeff = false);
  phase_diffusion (liq, f, varcoeff = false);
#endif

  phase_scalars_to_tracers (liq, f);
  phase_scalars_to_tracers (gas, f);

  scalar TL = liq->T, TG = gas->T;
  foreach()
    T[] = TL[] + TG[];
}

#if TWO_PHASE_VARPROP
event properties (i++) {
  scalar rhol = liq->rho, rhog = gas->rho;
  scalar mul = liq->mu, mug = gas->mu;
  foreach() {
    rho1v[] = rhol[];
    rho2v[] = rhog[];
    mu1v[] = mul[];
    mu2v[] = mug[];
  }
}
#endif

