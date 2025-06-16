#include "fractions.h"
#include "phase.h"
#include "tracer.h"

int NS;
double rhov = 1., muv = 1., lambda = 1., cp = 1.;
double T0 = 300., Pref = 101325., Dmix = 1., MW = 1.;
bool isothermal = false, isomassfrac = false;

Phase * phase;
scalar * tracers = NULL;
char ** species = NULL;

scalar rhot[];
face vector mut[];

struct CombustionModel {
  bool normalize;
  bool fick_corrected;
  bool molar_diffusion;
  bool heat_species_diffusion;
  bool chemistry;
  bool divergence;
  bool no_advection_div;
} combustion = {
  true,
  true,
  true,
  true,
  true,
  true,
  true,
};

attribute {
  bool inverse;
}

event defaults (i = 0) {
  phase = new_phase ("", NS, false, species);

  phase->isothermal = isothermal;
  phase->isomassfrac = isomassfrac;

  double * x = malloc (phase->n*sizeof (double));
  double * Ds = malloc (phase->n*sizeof (double));
  double * cps = malloc (phase->n*sizeof (double));
  double * MWs = malloc (phase->n*sizeof (double));

  foreach_species_in (phase)
    x[i] = 1.;
  correctfrac (x, phase->n);

  ThermoState ts;
  ts.T = T0, ts.P = Pref, ts.x = x;
  phase_set_thermo_state (phase, &ts);

  foreach_species_in (phase) {
    Ds[i] = Dmix;
    cps[i] = cp;
    MWs[i] = MW;
  }

  phase_set_properties (phase,
      rho = rhov, mu = muv,
      lambda = lambda, cp = cp,
      D = Ds, cps = cps, MWs = MWs);

  phase_set_tracers (phase);
  tracers = phase->tracers;

#if TREE
 for (scalar s in tracers) {
#if EMBED
   s.refine = s.prolongation = refine_embed_linear;
#else
   s.refine  = refine_linear;
#endif
   s.restriction = restriction_volume_average;
   s.dirty = true; // boundary conditions need to be updated
 }
#endif

  free (x);
  free (Ds);
  free (cps);
  free (MWs);

  // Default properties
  rho = rhot;
  mu = mut;
}

event init (i = 0) {
  assert (nv == 1);

  no_advection_div = combustion.no_advection_div;

  //phase_update_properties (phase, &tp2);
}

event cleanup (t = end) {
  delete_phase (phase);
}

event reset_sources (i++) {
  phase_reset_sources (phase);

  foreach()
    for (scalar drhodt in drhodtlist)
      drhodt[] = 0.;
}

event phase_properties (i++) {
  phase_update_mw_moles (phase);
  phase_update_properties (phase, &tp2);

  scalar rhop = phase->rho;
  foreach()
    rhot[] = rhop[]*cm[];

  scalar mup = phase->mu;
  foreach_face()
    mut.x[] = face_value (mup, 0)*fm.x[];
}

#if CHEMISTRY
event chemistry (i++) {
  if (combustion.chemistry) {
    ode_function batch = batch_isothermal_constantpressure;
    unsigned int NEQ = phase->n;
    if (!phase->isothermal) {
      batch = batch_nonisothermal_constantpressure;
      NEQ++;
    }
    phase_chemistry_direct (phase, dt, batch, NEQ);
  }

  if (combustion.heat_species_diffusion)
    phase_add_heat_species_diffusion (phase,
        molar_diffusion = combustion.molar_diffusion);
}
#endif

event divergence (i++) {
  if (combustion.divergence) {
    phase_update_mw_moles (phase);
    phase_update_divergence (phase,
        fick_corrected = combustion.fick_corrected,
        molar_diffusion = combustion.molar_diffusion);
    //if (i > 1) phase_update_divergence_density (phase, u);

    scalar divp = phase->divu;
    foreach()
      drhodt[] = divp[];
  }
}

event tracer_advection (i++);

event tracer_diffusion (i++) {
  phase_update_mw_moles (phase);

  phase_diffusion_velocity (phase,
      fick_corrected = combustion.fick_corrected,
      molar_diffusion = combustion.molar_diffusion);

  phase_diffusion (phase, varcoeff = true);

  if (combustion.normalize)
    phase_normalize_fractions (phase);
}

