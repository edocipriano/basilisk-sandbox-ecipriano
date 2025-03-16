scalar mEvapTot[];
extern scalar * mEvapList;

#include "common-evaporation.h"
#include "phase.h"
#include "regression.h"

Phase * liq, * gas;

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
} pcm = {
  ADVECTION_VELOCITY,
  SHIFT_TO_LIQUID,
  WITHOUT_EXPANSION,
  false,
  false,
};

void intexp_explicit (scalar intexp, scalar f, scalar mEvapTot) {
  foreach()
    intexp[] = 0.;

  foreach_interfacial_plic (f, F_ERR)
    intexp[] = (rho1 > 0. && rho2 > 0.) ?
      mEvapTot[]*(1./rho2 - 1./rho1)*dirac : 0.;
}

event init (i = 0) {
#if VELOCITY_JUMP
  _boiling = pcm.boiling;
#endif
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
  foreach()
    jump[] = -mEvapTot[]*(1./rho1 - 1./rho2);

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
  switch (pcm.advection) {
    case NO_ADVECTION:
      break;
    case ADVECTION_VELOCITY:
      vof_advection_phasechange (f, mEvapTot, i, byrhogas = pcm.byrhogas);
      break;
    case SOURCE_TERM:
      vof_expl_sources (f, mEvapTot, dt, byrhogas = pcm.byrhogas);
      break;
    case PLANE_SHIFTING:
      vof_plane_shifting (f, mEvapTot, dt, byrhogas = pcm.byrhogas);
  }

  if (nv > 1) {
    face vector uf1 = uflist[1];
    foreach_face() {
      ufsave.x[] = uf.x[];
      uf.x[] = uf1.x[];
    }
  }
}

event vof_sources (i++) {
  if (nv > 1)
    foreach_face()
      uf.x[] = ufsave.x[];

  foreach()
    f[] = clamp (f[], 0., 1.);
}

event tracer_advection (i++);

