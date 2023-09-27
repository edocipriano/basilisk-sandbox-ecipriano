/**
# Temperature Gradient Phase Change Model

This phase change model is suitable for boiling conditions.
The vaporization rate is computed from an energy balance at
the interface (in every interfacial cell), assuming that
the interface is at saturation temperature and that the
heat conduction to the interface is used for the phase change
phenomena. Therefore, the vaporization rate is computed from
the following balance:

$$
\dot{m}_{Evap} \Delta h_{ev} =
    \mathbf{\dot{q}}_l \cdot \mathbf{n_\Gamma}
  + \mathbf{\dot{q}}_g \cdot \mathbf{n_\Gamma} =
  - \lambda_l \left.\dfrac{\partial T_l}{\partial \mathbf{n}_\Gamma}\right\vert_l
  - \lambda_g \left.\dfrac{\partial T_g}{\partial \mathbf{n}_\Gamma}\right\vert_g
$$

where $\lambda$ is the thermal conductivity, $\Delta h_{ev}$
is the enthalpy of evaporation and $\mathbf{n_\Gamma}$ is the
normal to the interface.
After the calculation of the vaporization rate, this model
solves the transport equation for the temperature field $T$
using a two-field formulation and including the presence of
the phase change, both in the diffusion equation and in the
transport of the temperature field which keeps into account
the Stefan convection.
*/

#include "intgrad.h"
#include "fracface.h"
#include "diffusion.h"

/**
## Memory Allocations

This phase change model defines a list of vaporization
rates *mEvapList* which is populated by a single scalar
*mEvap* because a single contribution to the vaporization
rate is considered.
*/

#ifndef PHASECHANGE
(const) scalar mEvap = zeroc;
scalar * mEvapList;

/**
*fL* and *fG* store the value of volume fractions for the
calculation of the interface gradients, while the vectors
*fsL* and *fsG* contain the face fraction fields computed
using the [fracface.h](/sandbox/ecipriano/src/fracface.h)
module. */

scalar fL[], fG[];
face vector fsL[], fsG[];

face vector s[];
# define PHASECHANGE
#endif

/**
## Field Allocations

* *T* one-field temperature
* *TL* liquid-phase temperature field
* *TG* gas-phase temperature field
* *TInt* value of temperature at the interface
*/

scalar T[], TL[], TG[], TInt[];

#ifdef VARPROP
scalar rhovold[];
scalar frhocp1[], frhocp2[];
scalar frhocp1old[], frhocp2old[];
#endif

/**
## User Data

Using this phase change model, the user should define the
following variables (SI units):

* *lambda1* Thermal conductivity in liquid phase
* *lambda2* Thermal conductivity in gas phase
* *dhev* Latent heat of vaporization
* *cp1* Specific heat capacity in liquid phase
* *cp2* Specific heat capacity in gas phase
* *TIntVal* Interface temperature value
*/

extern double lambda1, lambda2, dhev, cp1, cp2;
extern double TIntVal;
double Pref = 101325.;

/**
Allocating diffusivity fields *lambdaf*, volume correction
*theta* and explicit source terms *sgT* and *sTS* for
the diffusion equation. */

face vector lambda1f[], lambda2f[];
scalar thetacorr1[], thetacorr2[];
scalar sgT[], slT[], sgTimp[], slTimp[];

/**
## Useful Functions

We define the functions that update the fields
for the material properties.
*/

static void update_properties_constant (void) {
  double x[] = {1.};
  ts1.T = 300., ts1.P = Pref, ts1.x = x;
  ts2.T = 300., ts2.P = Pref, ts2.x = x;

  foreach() {
    rho1v[] = rho1;
    rho2v[] = rho2;
    mu1v[] = mu1;
    mu2v[] = mu2;
    cp1v[] = cp1;
    cp2v[] = cp2;
    lambda1v[] = lambda1;
    lambda2v[] = lambda2;

    frhocp1[] = f[]*rho1v[]*cp1v[];
    frhocp2[] = (1. - f[])*rho2v[]*cp2v[];

    rhovold[] = rhov[];
    rhov[] = aavg (f[], rho1v[], rho2v[]);
  }
  boundary({lambda1v, lambda2v});
}

static void update_properties (void) {
  double x[] = {1.};
  ts1.T = 300., ts1.P = Pref, ts1.x = x;
  ts2.T = 300., ts2.P = Pref, ts2.x = x;

  foreach() {
    //ts1.T = T[], ts2.T = T[];
    double * There1 = &ts1.T;
    double * There2 = &ts2.T;
    *There1 = T[];
    *There2 = T[];

    //rho1v[] = rho1;
    rho1v[] = tp1.rhov (&ts1);
    //mu1v[] = tp1.muv (&ts1);
    //cp1v[] = tp1.cpv (&ts1);
    //lambda1v[] = tp1.lambdav (&ts1);

    rho2v[] = tp2.rhov (&ts2);
    //mu2v[] = tp2.muv (&ts2);
    //cp2v[] = tp2.cpv (&ts2);
    //lambda2v[] = tp2.lambdav (&ts2);

    frhocp1[] = f[]*rho1v[]*cp1;
    frhocp2[] = (1. - f[])*rho2v[]*cp2;
  }
}

/**
## Defaults

In the defaults event we setup the tracer lists for the
advection of the temperature fields. */

event defaults (i = 0)
{
  TL.inverse = false;
  TG.inverse = true;

  /**
  On adaptive meshes, tracers need to use linear interpolation (rather
  than the default bilinear interpolation) to ensure conservation when
  refining cells. */

#if TREE
#if EMBED
      TL.refine = TL.prolongation = refine_embed_linear;
      TG.refine = TG.prolongation = refine_embed_linear;
#else
      TL.refine  = refine_linear;
      TG.refine  = refine_linear;
#endif
      TL.restriction = restriction_volume_average;
      TL.dirty = true; // boundary conditions need to be updated
      TG.restriction = restriction_volume_average;
      TG.dirty = true; // boundary conditions need to be updated
#endif
}

/**
## Init

In the init event, we avoid dumping all the fields that we
don't need to visualize. */

event init (i = 0)
{
  sgT.nodump = true;
  slT.nodump = true;
  fL.nodump  = true;
  fG.nodump  = true;
  thetacorr1.nodump = true;
  thetacorr2.nodump = true;

  if (is_constant (rho2v)) {
    scalar * l = list_copy (all);
    rho2v = new scalar;
    free (all);
    all = list_concat ({rho2v}, l);
    free (l);
  }
  if (is_constant (mu2v)) {
    scalar * l = list_copy (all);
    mu2v = new scalar;
    free (all);
    all = list_concat ({mu2v}, l);
    free (l);
  }
  if (is_constant (cp2v)) {
    scalar * l = list_copy (all);
    cp2v = new scalar;
    free (all);
    all = list_concat ({cp2v}, l);
    free (l);
  }
  if (is_constant (lambda2v)) {
    scalar * l = list_copy (all);
    lambda2v = new scalar;
    free (all);
    all = list_concat ({lambda2v}, l);
    free (l);
  }
  if (is_constant (rho1v)) {
    scalar * l = list_copy (all);
    rho1v = new scalar;
    free (all);
    all = list_concat ({rho1v}, l);
    free (l);
  }
  if (is_constant (mu1v)) {
    scalar * l = list_copy (all);
    mu1v = new scalar;
    free (all);
    all = list_concat ({mu1v}, l);
    free (l);
  }
  if (is_constant (cp1v)) {
    scalar * l = list_copy (all);
    cp1v = new scalar;
    free (all);
    all = list_concat ({cp1v}, l);
    free (l);
  }
  if (is_constant (lambda1v)) {
    scalar * l = list_copy (all);
    lambda1v = new scalar;
    free (all);
    all = list_concat ({lambda1v}, l);
    free (l);
  }
  update_properties_constant();
}

/**
## Finalise

We deallocate the various lists from the memory. */

event cleanup (t = end);

/**
## Phase Change

In the *phasechange* event, the vaporization rate is computed
and the diffusion step for the mass fraction field (in liquid
and gas phase) is solved. */

event phasechange (i++)
{
  //update_properties();

  /**
  We compute the lagrangian derivative which will
  be introduced in the continuity equation in order
  to correct the divergence according to the temperature
  variations. */

  scalar lambdav[];
  foreach()
    lambdav[] = f[]*lambda1 + (1. - f[])*lambda2;
    //lambdav[] = f[]*lambda1v[] + (1. - f[])*lambda2v[];

  face vector lambdagT[];
  foreach_face() {
    double lambdavf = 0.5*(lambdav[] + lambdav[-1]);
    lambdagT.x[] = lambdavf*face_gradient_x (T, 0); /// !<<
  }

  // Compute lagrangian derivative
  foreach() {
    double laplT = 0.;
    foreach_dimension()
      laplT += (lambdagT.x[1] - lambdagT.x[]);
    laplT /= Delta;
    //ts1.T = T[], ts2.T = T[];
    double * There1 = &ts1.T;
    double * There2 = &ts2.T;
    *There1 = TL[];
    *There2 = T[];

    //double beta2exp = liqprop_thermal_expansion (&tp1, &ts1);
    double beta2exp = -1./T[];
    double drhodt1 = -beta2exp/(rho1v[]*cp1)*laplT;
    double drhodt2 = -1./(rho2v[]*cp2*T[])*laplT;

    drhodt[] = aavg (f[], drhodt1, drhodt2);
  }
}

/**
## Tracer Advection

We let the volume fractions *fu* and *fuext* to
advect the fields YL and YG, as implemented in
the tracer_advection event of [evaporation.h](evaporation.h)
*/

event tracer_advection (i++) {
  foreach() {
    //TL[] = (f[] > F_ERR) ? TL[]/f[] : 0.;
    //TG[] = (1. - f[] > F_ERR) ? TG[]/(1. - f[]) : 0.;
    //TL[] *= frhocp1[];
    //TG[] *= frhocp2[];
    //frhocp1old[] = frhocp1[];
    //frhocp2old[] = frhocp2[];
    TL[] *= rho1*cp1;
    TG[] *= rho2*cp2;
    frhocp1[] = f[]*rho1*cp1;
    frhocp2[] = (1. - f[])*rho2*cp2;
  }
}

/**
## Tracer Diffusion

We solve the diffusion equations for the temperature fields
accounting for the phase change contributions. */

event tracer_diffusion (i++)
{
  /**
  We remove the fractions of f and mass fractions
  lower than F_ERR and we reconstruct the non-volume
  averaged form of the mass fraction fields, in order
  to improve the discretization of the face gradients
  in the diffusion equation. */

  foreach() {
    f[] = clamp (f[], 0., 1.);
    fL[] = f[]; fG[] = 1. - f[];
    TL[] = (f[] > F_ERR) ? TL[]/frhocp1[] : 0.;
    TG[] = (1. - f[] > F_ERR) ? TG[]/frhocp2[] : 0.;
  }

  /**
  We compute the value of volume fraction *f* on the
  cell-faces using a geometric approach (necessary
  for interface gradients and diffusion equations). */

  face_fraction (fL, fsL);
  face_fraction (fG, fsG);

  /**
  The calculation of the interface gradients is used
  also for the calculation of the source terms for the
  diffusion equation of the temperature fields. */

  foreach() {
    sgT[] = 0., slT[] = 0.;
    if (f[] > F_ERR && f[] < 1.-F_ERR) {
      coord n = facet_normal (point, fL, fsL), p;
      double alpha = plane_alpha (fL[], n);
      double area = plane_area_center (n, alpha, &p);
      normalize (&n);

      double ltrgrad = ebmgrad (point, TL, fL, fG, fsL, fsG, false, TIntVal, false);
      double gtrgrad = ebmgrad (point, TG, fL, fG, fsL, fsG, true, TIntVal, false);

      double lheatflux = lambda1v[]*ltrgrad;
      double gheatflux = lambda2v[]*gtrgrad;

#ifdef AXI
      slT[] = lheatflux*area*(y + p.y*Delta)/(Delta*y)*cm[];
      sgT[] = gheatflux*area*(y + p.y*Delta)/(Delta*y)*cm[];
#else
      slT[] = lheatflux*area/Delta*cm[];
      sgT[] = gheatflux*area/Delta*cm[];
#endif
    }
  }

  /**
  We solve the diffusion equations, confined by means of
  the face fraction fields *fsL* and *fsG*. */

  foreach_face() {
    //double lambda1vf = 0.5*(lambda1v[] + lambda1v[-1]);
    //double lambda2vf = 0.5*(lambda2v[] + lambda2v[-1]);
    //lambda1f.x[] = lambda1vf*fsL.x[]*fm.x[];
    //lambda2f.x[] = lambda2vf*fsG.x[]*fm.x[];
    lambda1f.x[] = lambda1*fsL.x[]*fm.x[];
    lambda2f.x[] = lambda2*fsG.x[]*fm.x[];
  }

  foreach() {
    thetacorr1[] = cm[]*max (frhocp1[], F_ERR);
    thetacorr2[] = cm[]*max (frhocp2[], F_ERR);
  }

  // Qui c'è da aggiungere beta
  foreach() {
    slTimp[] = -(frhocp1[] - frhocp1old[])/dt;
    sgTimp[] = -(frhocp2[] - frhocp2old[])/dt;
  }

  diffusion (TG, dt, D=lambda2f, r=sgT, theta=thetacorr2);
  diffusion (TL, dt, D=lambda1f, r=slT, theta=thetacorr1);

  foreach() {
    TL[] *= f[];
    TG[] *= (1. - f[]);
    T[] = TL[] + TG[];
  }
}

