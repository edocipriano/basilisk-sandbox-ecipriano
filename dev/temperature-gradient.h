#include "intgrad.h"

scalar mEvap[], * mEvapList = {mEvap}, T[];
scalar fu[], * fulist = NULL;

int NLS = 0, NGS = 0;
double lambda1 = 1., lambda2 = 1., dhev = 1., cp1 = 1., cp2 = 1.;
double TIntVal = 300., TG0 = 1., TL0 = 1., Pref = 101325.;

//Phase * liq, * gas;

event defaults (i = 0) {
  liq = new_phase ("L", NLS, false);
  gas = new_phase ("G", NGS, true);

  ThermoState tsl, tsg;
  tsl.T = TL0, tsl.P = Pref, tsl.x = (double[]){1};
  tsg.T = TG0, tsg.P = Pref, tsg.x = (double[]){1};

  phase_set_thermo_state (liq, &tsl);
  phase_set_thermo_state (gas, &tsg);

  phase_set_properties (liq,
      rho = rho1, mu = mu1,
      lambda = lambda1, cp = cp1,
      dhev = dhev);

  phase_set_properties (gas,
      rho = rho2, mu = mu2,
      lambda = lambda2, cp = cp2,
      dhev = dhev);

  phase_scalars_to_tracers (liq, f);
  phase_scalars_to_tracers (gas, f);

  phase_set_tracers (liq);
  phase_set_tracers (gas);

  fulist = list_add (fulist, fu);
  for (int i = 1; i < nv; i++) {
    scalar fu = new scalar;
    char name[80];
    sprintf (name, "fu%d", i);
    free (fu.name);
    fu.name = strdup (name);
    fulist = list_add (fulist, fu);
  }

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

  for (scalar fu in fulist) {
    for (int d = 0; d < nboundary; d++) {
      fu.boundary[d] = f.boundary[d];
      fu.boundary_homogeneous[d] = f.boundary_homogeneous[d];
    }
  }
}

event cleanup (t = end) {
  for (int i = 1; i < nv; i++) {
    scalar fu = fulist[i];
    delete ({fu});
  }
  for (scalar fu in fulist)
    free (fu.tracers), fu.tracers = NULL;
  free (fulist), fulist = NULL;

  delete_phase (liq);
  delete_phase (gas);
}

event reset_sources (i++) {
  phase_reset_sources (liq);
  phase_reset_sources (gas);

  scalar TL = liq->T, TG = gas->T;
  foreach() {
   mEvap[] = 0.;
    T[] = TL[] + TG[];
    for (scalar fu in fulist)
      fu[] = f[];
  }
}

// fixme: put these operations in a function, such that the used can easily move
// it from the phase change event
event phasechange (i++) {
  phase_tracers_to_scalars (liq, f, tol = F_ERR);
  phase_tracers_to_scalars (gas, f, tol = F_ERR);

  scalar fl[], fg[];
  foreach() {
    fl[] = f[];
    fg[] = 1. - f[];
  }

  face vector fsl[], fsg[];
  face_fraction (fl, fsl);
  face_fraction (fg, fsg);

  //// fixme as:
  //// face vector fsl[], fsg[];
  //// face_fraction (f, fsl, inverse = false);
  //// face_fraction (f, fsg, inverse = true);

  foreach_interfacial (f, F_ERR) {
    scalar TL = liq->T, TG = gas->T;
    double ltrgrad = ebmgrad (point, TL, fl, fg, fsl, fsg, false, TIntVal, false);
    double gtrgrad = ebmgrad (point, TG, fl, fg, fsl, fsg, true, TIntVal, false);

    scalar lambdal = liq->lambda, lambdag = gas->lambda;
    scalar deltahev = liq->dhev;

    mEvap[] += (deltahev[] > 0.) ? lambdal[]*ltrgrad/deltahev[] : 0.;
    mEvap[] += (deltahev[] > 0.) ? lambdag[]*gtrgrad/deltahev[] : 0.;
  }

  foreach_interfacial_plic (f, F_ERR) {
    scalar TL = liq->T, TG = gas->T;
    double ltrgrad = ebmgrad (point, TL, fl, fg, fsl, fsg, false, TIntVal, false);
    double gtrgrad = ebmgrad (point, TG, fl, fg, fsl, fsg, true, TIntVal, false);

    scalar slT = liq->STexp, sgT = gas->STexp;
    scalar lambdal = liq->lambda, lambdag = gas->lambda;

    slT[] += lambdal[]*ltrgrad*dirac;
    sgT[] += lambdag[]*gtrgrad*dirac;
  }

  phase_scalars_to_tracers (liq, f);
  phase_scalars_to_tracers (gas, f);
}

event vof (i++) {
  // Transport tracers with uf by default (i.e. nv == 1)
  scalar fug = fulist[0];
  vof_advection ({fug}, i);

  // Transport tracers with liquid velocity
  if (nv == 2) {
    scalar ful = fulist[1];
    face vector ufsave[], ufl = uflist[1];
    foreach_face() {
      ufsave.x[] = uf.x[];
      uf.x[] = ufl.x[];
    }
    vof_advection ({ful}, i);
    foreach_face()
      uf.x[] = ufsave.x[];
  }
}

event tracer_advection (i++);

event tracer_diffusion (i++) {
  foreach() {
    f[] = clamp (f[], 0., 1.);
    if (f[] < F_ERR)
      f[] = 0.;
    if (f[] > 1. - F_ERR)
      f[] = 1.;
  }

  if (nv == 1) {
    phase_tracers_to_scalars (liq, fu, tol = F_ERR);
    phase_tracers_to_scalars (gas, fu, tol = F_ERR);
  }
  else if (nv == 2) {
    scalar ful = fulist[1], fug = fulist[0];
    phase_tracers_to_scalars (liq, ful, tol = F_ERR);
    phase_tracers_to_scalars (gas, fug, tol = F_ERR);
  }

  phase_diffusion (gas, f, varcoeff = false);
  phase_diffusion (liq, f, varcoeff = false);

  phase_scalars_to_tracers (liq, f);
  phase_scalars_to_tracers (gas, f);

  scalar TL = liq->T, TG = gas->T;
  foreach()
    T[] = TL[] + TG[];
}

