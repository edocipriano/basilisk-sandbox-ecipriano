/**
# Phase Change

Building blocks of the phase change solver. Each function performs one
operation of the time step, and the phase change models
([evaporation.h](evaporation.h), [boiling.h](boiling.h),
[fixedflux.h](fixedflux.h)) assemble them in the standard Basilisk events `vof`
and `tracer_diffusion`, in the order that suits their interface conditions. */

#include "common-evaporation.h"
#include "phase.h"

scalar mEvapTot[];
extern scalar * mEvapList;

int NLS = 0, NGS = 0;
double lambda1 = 1., lambda2 = 1., dhev = 1., cp1 = 1., cp2 = 1.;
double TG0 = 300., TL0 = 300., TIntVal = 300.;
double Dmix1 = 0., Dmix2 = 1., YIntVal = 0., YG0 = 0., YL0 = 1.;
double MW1 = 1., MW2 = 1.;

Phase * liq, * gas;
char ** liq_species = NULL, ** gas_species = NULL;

// volume fractions transporting the gas (fu) and liquid (fu1) for nv = 2
scalar fu[], fu1[], * fulist = {fu,fu1};
scalar T[], Y[];

#include "regression.h"

#ifndef vof_gradients
# define vof_gradients true
#endif

static scalar * f_tracers = NULL;

/**
## Parameters

The structure `pcm` collects the options of the phase change model: the
transport of the interface due to the phase change (`advection`, `byrhogas`),
the localization of the expansion term (`expansion`, `shifting`), the
treatment of the diffusion (`fick_corrected`, `molar_diffusion`,
`mass_diffusion_enthalpy`, `normalize`), and the physics
included in the model (`boiling`, `isothermal`, `isomassfrac`,
`isothermal_interface`, `divergence`, `chemistry`, `emissivity`). The models
set their own values in `init`, which can be overridden by the simulation in
`end_init`. */

enum advection_policy {
  NO_ADVECTION, ADVECTION_VELOCITY, SOURCE_TERM, PLANE_SHIFTING
};

enum shifting_policy {
  NO_SHIFTING, SHIFT_TO_LIQUID, SHIFT_TO_GAS
};

struct PhaseChangeModel {
  enum advection_policy advection;
  enum shifting_policy shifting;
  bool boiling;
  bool byrhogas;
  bool expansion;
  bool isothermal;
  bool isomassfrac;
  bool isothermal_interface;
  bool fick_corrected;
  bool molar_diffusion;
  bool mass_diffusion_enthalpy;
  bool divergence;
  bool chemistry;
  bool normalize;
  double emissivity;
} pcm = {
  ADVECTION_VELOCITY,         // advection
  SHIFT_TO_LIQUID,            // shifting
  false,                      // boiling
  false,                      // byrhogas
  true,                       // expansion
  false,                      // isothermal
  false,                      // isomassfrac
  false,                      // isothermal_interface
#if TWO_PHASE_VARPROP
  true,                       // fick_corrected
  true,                       // molar_diffusion
  true,                       // mass_diffusion_enthalpy
  true,                       // divergence
#else
  false,                      // fick_corrected
  false,                      // molar_diffusion
  false,                      // mass_diffusion_enthalpy
  false,                      // divergence
#endif
  false,                      // chemistry
  true,                       // normalize
  0.,                         // emissivity
};

/**
## Functions

### *intexp_explicit()*: expansion term due to the phase change */

void intexp_explicit (scalar intexp, scalar f, scalar mEvapTot) {
  scalar rhol = liq->rho, rhog = gas->rho;
  foreach_interfacial_plic (f, F_ERR)
    intexp[] = (rhol[] > 0. && rhog[] > 0.) ?
      mEvapTot[]*(1./rhog[] - 1./rhol[])*dirac : 0.;
}

/**
### *phasechange_reset()*: set to zero vaporization rates and expansion terms */

void phasechange_reset (void) {
  foreach() {
    for (scalar mEvap in mEvapList)
      mEvap[] = 0.;
    for (scalar intexp in intexplist)
      intexp[] = 0.;
    for (scalar drhodt in drhodtlist)
      drhodt[] = 0.;
  }
}

/**
### *phasechange_reset_sources()*: set to zero the phase source terms

The source terms of the phases (`STexp`, `SYexp`, ...) are reset after the
diffusion, so that the simulation can add its own sources before the interface
moves (`vof` event) or after (`tracer_diffusion` event). */

void phasechange_reset_sources (void) {
  phase_reset_sources (liq);
  phase_reset_sources (gas);
}

/**
### *phasechange_properties()*: update the material properties of the phases */

// fixme: I should be able to compute MW and Moles for the molar diffusion even
// if the model is used with constant properties
void phasechange_properties (void) {
#if TWO_PHASE_VARPROP
  phase_tracers_to_scalars (liq, f, tol = F_ERR);
  phase_tracers_to_scalars (gas, f, tol = F_ERR);

  /**
  In a closed system the thermodynamic pressure changes in time, and the
  material properties must be updated using its current value. */

#if LOW_MACH
  if (closed)
    foreach() {
      foreach_scalar_in (liq)
        P[] = P0;
      foreach_scalar_in (gas)
        P[] = P0;
    }
#endif

  phase_update_mw_moles (liq, f, tol = P_ERR, extend = true);
  phase_update_mw_moles (gas, f, tol = P_ERR, extend = true);

  phase_update_properties (liq, &tp1, f, P_ERR);
  phase_update_properties (gas, &tp2, f, P_ERR);

  phase_extend_properties (liq, f, P_ERR);
  phase_extend_properties (gas, f, P_ERR);

#if LOW_MACH
  /**
  The work performed by the variation of the thermodynamic pressure is added to
  the temperature equation, see `phase_add_compression_work()`. The rate of the
  current time step is known only after the projection, therefore we use the one
  computed by the projection of the previous time step.

  The source is added here, together with the material properties, because
  `betaT` is evaluated by `phase_update_properties()` using the temperature of
  this function: the product $\beta_T T$ must be formed from the same
  thermodynamic state, and it reduces to one for an ideal gas. The source must
  also be available before *phasechange_divergence()*, in order to be included in
  the material derivative of the temperature computed by
  `phase_update_divergence()`: without this contribution the closed-system
  projection underestimates the pressurization rate by a factor $\gamma$. */

  if (closed) {
    phase_add_compression_work (liq, dP0dt, f, F_ERR);
    phase_add_compression_work (gas, dP0dt, f, F_ERR);
  }
#endif

  phase_scalars_to_tracers (liq, f);
  phase_scalars_to_tracers (gas, f);
#endif
}

/**
### *phasechange_chemistry()*: gas phase reactions */

void phasechange_chemistry (void) {
#if CHEMISTRY
  if (pcm.chemistry) {
    ode_function batch = batch_isothermal_constantpressure;
    unsigned int NEQ = gas->n;
    if (!gas->isothermal) {
      batch = batch_nonisothermal_constantpressure;
      NEQ++;
    }
    phase_chemistry_direct (gas, dt, batch, NEQ, f, tol = 1-F_ERR);
  }
#endif
}

/**
### *phasechange_rate()*: vaporization rates and interface source terms

This function is defined by the model, and it computes the vaporization rates
`mEvapList` and the interface source terms of the phases. */

void phasechange_rate (void);

/**
### *phasechange_expansion()*: total vaporization rate and expansion term */

void phasechange_expansion (void) {
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

  scalar fext[];
  foreach()
    fext[] = f[];
  vof_to_ls (fext, ls, imax = 5);

  constant_extrapolation (jump, ls, 0.5, 10, c=fext, nl=0,
            nointerface=true, inverse=false, tol=1e-4);
  constant_extrapolation (jump, ls, 0.5, 10, c=fext, nl=0,
            nointerface=true, inverse=true, tol=1e-4);

  foreach_face()
    jumpf.x[] = face_value (jump, 0);
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

/**
### *phasechange_divergence()*: velocity divergence due to density changes */

void phasechange_divergence (void) {
#if TWO_PHASE_VARPROP
  if (pcm.divergence) {
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
#endif
}

/**
### *phasechange_regression()*: phase change step of the interface

The phase change moves the interface with the total vaporization rate
`mEvapTot`. */

void phasechange_regression (void) {
  // fixme: add consistency also for source and plane shifting
  scalar rhoh = pcm.byrhogas ? gas->rho : liq->rho;
  switch (pcm.advection) {
    case NO_ADVECTION:
      break;
    case ADVECTION_VELOCITY:
      vof_advection_phasechange (f, mEvapTot, rhoh, iter);
      break;
    case SOURCE_TERM:
      vof_expl_sources (f, mEvapTot, rhoh, dt);
      break;
    case PLANE_SHIFTING:
      vof_plane_shifting (f, mEvapTot, rhoh, dt);
  }
}

/**
### *phasechange_transport()*: transport of the tracers

The volume fractions `fulist`, which carry the phase tracers, are transported
by the fluid velocity. With `nv == 2`, the gas tracers are transported with
the gas velocity by *fu*, and the liquid tracers with the liquid velocity by
*fu1*, while *f* is transported with the liquid velocity by
[vof.h](/src/vof.h), which is set here and restored by
*phasechange_advection_end()*. With `nv == 1`, the phase tracers are
transported with *f*. */

face vector ufsave[];

void phasechange_transport (void) {
  // Volume fractions used to transport the tracers
  foreach()
    for (scalar fu in fulist)
      fu[] = f[];

  if (nv == 2) {
    scalar fug = fulist[0], ful = fulist[1];
    vof_advection ({fug}, iter);

    face vector ufl = uflist[1];
    foreach_face() {
      ufsave.x[] = uf.x[];
      uf.x[] = ufl.x[];
    }
    boundary ((scalar *){uf});  // do not remove
    vof_advection ({ful}, iter);
  }
}

/**
### *phasechange_advection_end()*: restore the velocity and clamp the interface */

void phasechange_advection_end (void) {
  if (nv == 2)
    foreach_face()
      uf.x[] = ufsave.x[];
    boundary ((scalar *){uf});  // do not remove

  foreach() {
    f[] = clamp (f[], 0., 1.);
    f[] = (f[] < F_ERR) ? 0. : (f[] > 1.-F_ERR) ? 1. : f[];
  }
}

/**
### *phasechange_diffusion()*: diffusion of temperature and species

The phase tracers are diffused, and the one-field temperature `T` and mass
fraction `Y` are reconstructed. */

void phasechange_diffusion (void) {

  if (nv == 1) {
    phase_tracers_to_scalars (liq, f, tol = F_ERR);
    phase_tracers_to_scalars (gas, f, tol = F_ERR);
  }
  else {
    scalar ful = fulist[1], fug = fulist[0];
    phase_tracers_to_scalars (liq, ful, tol = F_ERR);
    phase_tracers_to_scalars (gas, fug, tol = F_ERR);
  }

  phase_update_mw_moles (liq, f, tol = P_ERR, extend = true);
  phase_update_mw_moles (gas, f, tol = P_ERR, extend = true);

  // [DIFF]
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

  if (pcm.normalize) {
    phase_normalize_fractions (liq);
    phase_normalize_fractions (gas);
  }

  phase_scalars_to_tracers (liq, f);
  phase_scalars_to_tracers (gas, f);

  scalar TL = liq->T, TG = gas->T;
  foreach() {
    T[] = 0.;
    T[] += liq->isothermal ? f[]*TL[] : TL[];
    T[] += gas->isothermal ? (1. - f[])*TG[] : TG[];
  }

  /**
  The one-field mass fraction sums the species of both phases with the same
  index, which are the transferring species if they come first in both
  lists. Only the first *min (NLS, NGS)* species exist in both phases.
  fixme: match the species by name. */

  int nmin = min (liq->n, gas->n);
  foreach() {
    double Ysum = 0.;
    for (int i = 0; i < nmin; i++) {
      scalar YL = liq->YList[i], YG = gas->YList[i];
      Ysum += liq->isomassfrac ? f[]*YL[] : YL[];
      Ysum += gas->isomassfrac ? (1. - f[])*YG[] : YG[];
    }
    Y[] = Ysum;
  }
}

/**
### *phasechange_mixture()*: mixture properties for the Navier--Stokes equations */

void phasechange_mixture (void) {
#if TWO_PHASE_VARPROP
  scalar rhol = liq->rho, rhog = gas->rho;
  scalar mul = liq->mu, mug = gas->mu;
#if LOW_MACH
  scalar chiTl = liq->chiT, chiTg = gas->chiT;
#endif
  foreach() {
    rho1v[] = rhol[];
    rho2v[] = rhog[];
    mu1v[] = mul[];
    mu2v[] = mug[];

    /**
    The isothermal compressibility of the two phases is published into the
    field used by the closed-system projection. It is a volumetric quantity,
    therefore the two phases are combined with a volume-weighted average. */

#if LOW_MACH
    chiT[] = aavg (f[], chiTl[], chiTg[]);
#endif
  }
#endif
}

/**
## Events

The events `vof` and `tracer_diffusion` set the order of integration, through
the macros *phasechange_vof()* and *phasechange_tracer_diffusion()*. By default,
the phase change step precedes the transport of the interface and of the
tracers (phase change, VOF, advection, diffusion, Navier--Stokes). A model with
a different order overloads the two macros (see
[macro.h](/src/ast/macro.h#overloading)).

The simulation can add its own source terms to the phases in the `vof` event,
before the interface moves, or in the `tracer_diffusion` event, after the
interface moves.

We create the phases, set their properties, and attach the phase tracers to
the volume fractions used to transport them. */

event defaults (i = 0) {
  liq = new_phase ("L", NLS, false, liq_species);
  gas = new_phase ("G", NGS, true, gas_species);

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
  tsl.T = TL0, tsl.P = P0, tsl.x = xl;
  tsg.T = TG0, tsg.P = P0, tsg.x = xg;

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

  /**
  With `nv == 1`, the phase tracers are transported with *f*. With `nv == 2`,
  the gas and liquid tracers are transported with the gas and liquid
  velocities by *fu* and *fu1*, respectively. */

  if (nv == 1) {
    f_tracers = f.tracers;
    // These passages are necessary to avoid memory leaks
    scalar * f_tracers2 = list_concat (f.tracers, liq->tracers);
    f.tracers = list_concat (f_tracers2, gas->tracers);
    free (f_tracers2), f_tracers2 = NULL;
  }
  else {
    scalar ful = fulist[1], fug = fulist[0];
    ful.tracers = list_concat (ful.tracers, liq->tracers);
    fug.tracers = list_concat (fug.tracers, gas->tracers);
  }

#if TREE
  // Only f.tracers are set to the correct refine functions in
  // [vof.h](src/vof.h). Other helper fractions with their tracers must be
  // set here. Do not remove.
  for (scalar c in fulist) {
    c.refine = fraction_refine;
    set_prolongation (c, fraction_refine);
    scalar * tracers = c.tracers;
    for (scalar t in tracers) {
      set_restriction (t, restriction_volume_average);
      t.refine = vof_concentration_refine;
      set_prolongation (t, vof_concentration_refine);
      t.c = c;

      t.depends = list_add (t.depends, c);
    }
  }
#endif
}

/**
We set the options of the solver, and the initial material properties. */

event init (i = 0) {
#if VELOCITY_JUMP
  _boiling = pcm.boiling;
#endif

  if (pcm.divergence)
    no_advection_div = true;

  if (phase_is_uniform (liq))
    phase_scalars_to_tracers (liq, f);
  if (phase_is_uniform (gas))
    phase_scalars_to_tracers (gas, f);

#if TWO_PHASE_VARPROP
  phasechange_properties();

  // The compression work is added again in the first time step
  phasechange_reset_sources();
#endif
}

/**
We free the tracer lists and the phases. */

event cleanup (t = end) {
  for (scalar fu in fulist)
    free (fu.tracers), fu.tracers = NULL;

  if (nv == 1)
    free (f.tracers), f.tracers = f_tracers;

  delete_phase (liq);
  delete_phase (gas);
}

/**
The simulation can override the `pcm` parameters set by the model. */

event end_init (i = 0);

/**
Phase change step, and transport of the interface and of the tracers. */

auto macro phasechange_vof() {
  phasechange_reset();
  phasechange_properties();
  phasechange_chemistry();
  phasechange_rate();
  phasechange_expansion();
  phasechange_divergence();
  phasechange_transport();
  phasechange_regression();
}

event vof (i++) {
  phasechange_vof();
}

/**
We complete the transport of the interface, before the advection of the other
tracers. */

event tracer_advection (i++) {
  phasechange_advection_end();
}

/**
Diffusion of temperature and species, which uses and resets the source
terms. */

auto macro phasechange_tracer_diffusion() {
  phasechange_diffusion();
  phasechange_reset_sources();
}

event tracer_diffusion (i++) {
  phasechange_tracer_diffusion();
}

/**
We set the mixture properties used by the Navier--Stokes solver. */

#if TWO_PHASE_VARPROP
event properties (i++) {
  phasechange_mixture();
}
#endif
