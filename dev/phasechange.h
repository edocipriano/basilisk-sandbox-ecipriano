#include "common-evaporation.h"
#include "phase.h"

scalar mEvapTot[];
extern scalar * mEvapList;

int NLS = 0, NGS = 0;
double lambda1 = 1., lambda2 = 1., dhev = 1., cp1 = 1., cp2 = 1.;
double TG0 = 300., TL0 = 300., TIntVal = 300., Pref = 101325.;
double Dmix1 = 0., Dmix2 = 1., YIntVal = 0., YG0 = 0., YL0 = 1.;
double MW1 = 1., MW2 = 1.;

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

enum diffusion_policy {
  EXPLICIT_ONLY, EXPLICIT_IMPLICIT
};

struct PhaseChangeModel {
  enum advection_policy advection;
  enum shifting_policy shifting;
  enum velocity_policy velocity;
  enum diffusion_policy diffusion;
  bool boiling;
  bool byrhogas;
  bool expansion;
  bool consistent;
  bool isothermal;
  bool isomassfrac;
  bool isothermal_interface;
  bool fick_corrected;
  bool molar_diffusion;
  bool mass_diffusion_enthalpy;
  bool divergence_free;
  double emissivity;
} pcm = {
  ADVECTION_VELOCITY,
  SHIFT_TO_LIQUID,
  WITHOUT_EXPANSION,
  EXPLICIT_IMPLICIT,
  false,
  false,
  true,
  false,
  false,
  false,
  false,
  false,
  false,
  false,
  false,
  0.,
};

void intexp_explicit (scalar intexp, scalar f, scalar mEvapTot) {
  scalar rhol = liq->rho, rhog = gas->rho;
  foreach_interfacial_plic (f, F_ERR)
    intexp[] = (rhol[] > 0. && rhog[] > 0.) ?
      mEvapTot[]*(1./rhog[] - 1./rhol[])*dirac : 0.;
}

static scalar * f_tracers = NULL;

event defaults (i = 0) {
  liq = new_phase ("L", NLS, false);
  gas = new_phase ("G", NGS, true);

  liq->isothermal = pcm.isothermal;
  gas->isothermal = pcm.isothermal;
  liq->isomassfrac = pcm.isomassfrac;
  gas->isomassfrac = pcm.isomassfrac;

  double * xl = malloc (liq->n*sizeof (double));
  double * xg = malloc (gas->n*sizeof (double));
  double * Dl = malloc (liq->n*sizeof (double));
  double * Dg = malloc (gas->n*sizeof (double));
  double * dhevsl = malloc (liq->n*sizeof (double));
  double * dhevsg = malloc (gas->n*sizeof (double));
  double * cpsl = malloc (liq->n*sizeof (double));
  double * cpsg = malloc (gas->n*sizeof (double));
  double * MWl = malloc (liq->n*sizeof (double));
  double * MWg = malloc (gas->n*sizeof (double));

  foreach_species_in (liq)
    xl[i] = 1.;
  correctfrac (xl, liq->n);

  foreach_species_in (gas)
    xg[i] = 1.;
  correctfrac (xg, gas->n);

  if (NLS == 1)
    xl[0] = 1.;

  if (NGS == 1)
    xg[0] = 1.;
  else if (NGS == 2)
    xg[0] = YG0, xg[1] = 1. - YG0;

  ThermoState tsl, tsg;
  tsl.T = TL0, tsl.P = Pref, tsl.x = xl;
  tsg.T = TG0, tsg.P = Pref, tsg.x = xg;

  phase_set_thermo_state (liq, &tsl);
  phase_set_thermo_state (gas, &tsg);

  foreach_species_in (liq) {
    Dl[i] = Dmix1;
    dhevsl[i] = dhev;
    cpsl[i] = cp1;
    MWl[i] = MW1;
  }

  foreach_species_in (gas) {
    Dg[i] = Dmix2;
    dhevsg[i] = dhev;
    cpsg[i] = cp2;
    MWg[i] = MW2;
  }

  phase_set_properties (liq,
      rho = rho1, mu = mu1,
      lambda = lambda1, cp = cp1,
      dhev = dhev, dhevs = dhevsl,
      D = Dl, cps = cpsl, MWs = MWl);

  phase_set_properties (gas,
      rho = rho2, mu = mu2,
      lambda = lambda2, cp = cp2,
      dhev = dhev, dhevs = dhevsg,
      D = Dg, cps = cpsg, MWs = MWg);

  free (xl);
  free (xg);
  free (Dl);
  free (Dg);
  free (dhevsl);
  free (dhevsg);
  free (cpsl);
  free (cpsg);
  free (MWl);
  free (MWg);

  phase_set_tracers (liq);
  phase_set_tracers (gas);

  if (pcm.consistent) {
    f_tracers = f.tracers;
    f.tracers = list_concat (f.tracers, liq->tracers);
    f.tracers = list_concat (f.tracers, gas->tracers);
  }
  else {
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
  }

#if TREE
  // Only f.tracers are set to the correct refine functions in
  // [vof.h](src/vof.h). Other helper fractions with their tracers must be
  // set here. Do not remove.
  for (scalar c in fulist) {
    c.refine = c.prolongation = fraction_refine;
    c.dirty = true;
    scalar * tracers = c.tracers;
    for (scalar t in tracers) {
      t.restriction = restriction_volume_average;
      t.refine = t.prolongation = vof_concentration_refine;
      t.dirty = true;
      t.c = c;

      t.depends = list_add (t.depends, c);
    }
  }
#endif
}

event init (i = 0) {
#if VELOCITY_JUMP
  _boiling = pcm.boiling;
#endif

#if TWO_PHASE_VARPROP
  no_advection_div = true;
  pcm.divergence_free = false;
  pcm.fick_corrected = true;
  pcm.molar_diffusion = true;
  pcm.mass_diffusion_enthalpy = true;
#endif

  if (phase_is_uniform (liq))
    phase_scalars_to_tracers (liq, f);
  if (phase_is_uniform (gas))
    phase_scalars_to_tracers (gas, f);
}

event cleanup (t = end) {
  for (scalar fu in fulist)
    free (fu.tracers), fu.tracers = NULL;

  if (pcm.consistent)
    free (f.tracers), f.tracers = f_tracers;

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

// fixme: I should be able to compute MW and Moles for the molar diffusion even
// if the model is used with constant properties
#if TWO_PHASE_VARPROP
event phase_properties (i++) {
  phase_tracers_to_scalars (liq, f, tol = F_ERR);
  phase_tracers_to_scalars (gas, f, tol = F_ERR);

  phase_update_mw_moles (liq, f, tol = P_ERR, extend = true);
  phase_update_mw_moles (gas, f, tol = P_ERR, extend = true);

  phase_update_properties (liq, &tp1, f, P_ERR);
  phase_update_properties (gas, &tp2, f, P_ERR);

  phase_extend_properties (liq, f, P_ERR);
  phase_extend_properties (gas, f, P_ERR);

  phase_scalars_to_tracers (liq, f);
  phase_scalars_to_tracers (gas, f);
}
#endif

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
  if (pcm.expansion) {
    intexp_explicit (intexp, f, mEvapTot);

    switch (pcm.shifting) {
      case NO_SHIFTING: break;
      case SHIFT_TO_LIQUID: shift_field (intexp, f, 1); break;
      case SHIFT_TO_GAS: shift_field (intexp, f, 0); break;
    }
  }
#endif
}

#if TWO_PHASE_VARPROP
event divergence (i++) {
  if (!pcm.divergence_free) {
    phase_tracers_to_scalars (liq, f, tol = F_ERR);
    phase_tracers_to_scalars (gas, f, tol = F_ERR);

    phase_update_divergence (liq, f, pcm.fick_corrected, pcm.molar_diffusion);
    phase_update_divergence (gas, f, pcm.fick_corrected, pcm.molar_diffusion);

    scalar divu1 = liq->divu, divu2 = gas->divu;
    foreach()
      drhodt[] = divu1[] + divu2[];

    if (nv > 1) {
      scalar drhodt1 = drhodtlist[1];
      foreach()
        drhodt1[] = divu1[];
    }

    phase_scalars_to_tracers (liq, f);
    phase_scalars_to_tracers (gas, f);
  }
}
#endif

face vector ufsave[];

event vof (i++) {
  // Transport due to the phase change
  // fixme: add consistency also for source and plane shifting
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

  if (pcm.consistent) {
    phase_tracers_to_scalars (liq, f, tol = F_ERR);
    phase_tracers_to_scalars (gas, f, tol = F_ERR);
  }
  else {
    if (nv == 1) {
      phase_tracers_to_scalars (liq, fu, tol = F_ERR);
      phase_tracers_to_scalars (gas, fu, tol = F_ERR);
    }
    else if (nv == 2) {
      scalar ful = fulist[1], fug = fulist[0];
      phase_tracers_to_scalars (liq, ful, tol = F_ERR);
      phase_tracers_to_scalars (gas, fug, tol = F_ERR);
    }
  }

  phase_update_mw_moles (liq, f, tol = P_ERR, extend = true);
  phase_update_mw_moles (gas, f, tol = P_ERR, extend = true);

  // Let's perform the velocity correction step
  phase_diffusion_velocity (liq, f, pcm.fick_corrected, pcm.molar_diffusion);
  phase_diffusion_velocity (gas, f, pcm.fick_corrected, pcm.molar_diffusion);

#if TWO_PHASE_VARPROP
  phase_diffusion (gas, f, varcoeff = true);
  phase_diffusion (liq, f, varcoeff = true);
#else
  phase_diffusion (gas, f, varcoeff = false);
  phase_diffusion (liq, f, varcoeff = false);
#endif

  phase_normalize_fractions (liq);
  phase_normalize_fractions (gas);

  phase_scalars_to_tracers (liq, f);
  phase_scalars_to_tracers (gas, f);

  scalar TL = liq->T, TG = gas->T;
  foreach()
    T[] = TL[] + TG[];

  // fixme: make it general choosing the right gas species
  foreach() {
    double Ysum = 0.;
    foreach_species_in (liq) {
      scalar YL = liq->YList[i], YG = gas->YList[i];
      Ysum += YL[] + YG[];
    }
    Y[] = Ysum;
  }
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

