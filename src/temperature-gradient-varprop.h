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
//scalar mEvap[];
(const) scalar mEvap = zeroc;
scalar * mEvapList;

/**
*fL* and *fG* store the value of volume fractions for the
calculation of the interface gradients, while the vectors
*fsL* and *fsG* contain the face fraction fields computed
using the [fracface.h](/sandbox/ecipriano/src/fracface.h)
module. */

scalar fL[], fG[], f0[];
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
scalar frhocp1[], frhocp2[];
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

/**
Allocating diffusivity fields *lambdaf*, volume correction
*theta* and explicit source terms *sgT* and *sTS* for
the diffusion equation. */

face vector lambda1f[], lambda2f[];
scalar thetacorr1[], thetacorr2[];
scalar sgT[], slT[], sgTimp[], slTimp[];

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

  /**
  With variable properties, we initialize the propertiy
  fields in order to be non-constant. */
}

/**
## Init

In the init event, we avoid dumping all the fields that we
don't need to visualize. */

event init (i = 0,last)
{
  sgT.nodump = true;
  slT.nodump = true;
  fL.nodump  = true;
  fG.nodump  = true;
  thetacorr1.nodump = true;
  thetacorr2.nodump = true;

#ifdef VARPROP
  if (is_constant (rho2v)) {
    scalar * l = list_copy (all);
    rho2v = new scalar;
    free (all);
    all = list_concat ({rho2v}, l);
    free (l);
  }
  if (is_constant (cp2v)) {
    scalar * l = list_copy (all);
    cp2v = new scalar;
    free (all);
    all = list_concat ({cp2v}, l);
    free (l);
  }
  if (is_constant (lambda2v.x)) {
    scalar * l = list_copy (all);
    lambda2v = new face vector;
    free (all);
    all = list_concat ((scalar *){lambda2v}, l);
    free (l);
  }
#endif
}

/**
## Finalise

We deallocate the various lists from the memory. */

event cleanup (t = end)
{
  delete (fu.tracers), free (fu.tracers), fu.tracers = NULL;
  delete (fuext.tracers), free (fuext.tracers), fuext.tracers = NULL;
}

/**
## Phase Change

In the *phasechange* event, the vaporization rate is computed
and the diffusion step for the mass fraction field (in liquid
and gas phase) is solved. */

event phasechange (i++)
{

#ifdef VARPROP
  {
    double x[] = {1.};
    ts2.T = 300., ts2.P = 101325., ts2.x = x;
    rho2 = tp2.rhov (&ts2);

    // Update cell properties
    foreach() {
      ts2.T = T[];
      rho2v[] = tp2.rhov (&ts2);
      cp2v[] = tp2.cpv (&ts2);

      frhocp2[] = (1. - f[])*rho2v[]*cp2v[];
    }

    // Update face properties
    foreach_face() {
      ts2.T = 0.5*(T[] + T[-1]);
      lambda2v.x[] = tp2.lambdav (&ts2);
    }

    // Compute lagrangian derivative
    foreach() {
      double laplT = 0.;
      foreach_dimension()
        laplT += lambda2*face_gradient_x (T, 0);
      laplT /= Delta;
      drhodt[] = -1./(rho2v[]*cp2*T[])*laplT*(1. - f[]);
    }
  }
#endif

  /**
  First we compute the value of the non volume-averaged
  temperature fields. This procedure allows a better
  calculation of the gradients close to the interface. */

  foreach() {
    f[] = clamp (f[], 0., 1.);
    f[] = (f[] > F_ERR) ? f[] : 0.;
    f0[] = f[];
    fL[] = f[]; fG[] = 1. - f[];
    TL[] = f[] > F_ERR ? TL[]/f[] : 0.;
    TG[] = ((1. - f[]) > F_ERR) ? TG[]/(1. - f[]) : 0.;
  }
  //boundary({fL,fG,TL,TG,f0});

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

#ifdef VARPROP
      {
        double x[] = {1.};
        ts2.T = T[], ts2.P = 101325., ts2.x = x;
      }
      double lheatflux = lambda1*ltrgrad;
      double gheatflux = tp2.lambdav (&ts2)*gtrgrad;
#else
      double lheatflux = lambda1*ltrgrad;
      double gheatflux = lambda2*gtrgrad;
#endif

#ifdef AXI
      slT[] = lheatflux/rho1/cp1*area*(y + p.y*Delta)/(Delta*y)*cm[];
#ifdef VARPROP
      sgT[] = gheatflux*area*(y + p.y*Delta)/(Delta*y)*cm[];
#else
      sgT[] = gheatflux/rho2/cp2*area*(y + p.y*Delta)/(Delta*y)*cm[];
#endif
#else
      slT[] = lheatflux/rho1/cp1*area/Delta*cm[];
#ifdef VARPROP
      sgT[] = gheatflux*area/Delta*cm[];
#else
      sgT[] = gheatflux/rho2/cp2*area/Delta*cm[];
#endif
#endif
    }
  }

  /**
  We restore the tracer form of the liquid and gas-phase
  temperature fields. */

  foreach() {
    TL[] *= f[]*(f[] > F_ERR);
    TG[] *= (1. - f[])*((1. - f[]) > F_ERR);
    T[]  = TL[] + TG[];
  }
}

/**
## Tracer Advection

We let the volume fractions *fu* and *fuext* to
advect the fields YL and YG, as implemented in
the tracer_advection event of [evaporation.h](evaporation.h)
*/

event tracer_advection (i++) {
#ifdef VARPROP
  foreach() {
    TG[] = (1. - f[]) > F_ERR ? TG[]/(1. - f[]) : 0.;
    TG[] *= frhocp2[];
  }
#endif
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
    f[] = (f[] > F_ERR) ? f[] : 0.;
    fuext[] = f[];
    fu[] = f[];

#ifdef CONSISTENTPHASE1
    TL[] = fuext[] > F_ERR ? TL[]/fuext[] : 0.;
#else
    TL[] = fu[] > F_ERR ? TL[]/fu[] : 0.;
#endif
#ifdef CONSISTENTPHASE2
    //TG[] = ((1. - fuext[]) > F_ERR) ? TG[]/(1. - fuext[]) : 0.;
    TG[] = ((1. - fuext[]) > F_ERR) ? TG[]/frhocp2[] : 0.;
#else
    //TG[] = ((1. - fu[]) > F_ERR) ? TG[]/(1. - fu[]) : 0.;
    TG[] = ((1. - fu[]) > F_ERR) ? TG[]/frhocp2[] : 0.;
#endif
    fL[] = f[]; fG[] = 1. - f[];
  }

  /**
  We compute the value of volume fraction *f* on the
  cell-faces using a geometric approach (necessary
  for interface gradients and diffusion equations). */

  face_fraction (fL, fsL);
  face_fraction (fG, fsG);

  /**
  We solve the diffusion equations, confined by means of
  the face fraction fields *fsL* and *fsG*. */

  foreach_face() {
    lambda1f.x[] = lambda1/rho1/cp1*fsL.x[]*fm.x[];
#ifdef VARPROP
    lambda2f.x[] = lambda2v.x[]*fsG.x[]*fm.x[];
#else
    lambda2f.x[] = lambda2/rho2/cp2*fsG.x[]*fm.x[];
#endif
  }

  foreach() {
    thetacorr1[] = cm[]*max(fL[], F_ERR);
#ifdef VARPROP
    thetacorr2[] = cm[]*max(fG[], F_ERR)*rho2v[]*cp2v[];
#else
    thetacorr2[] = cm[]*max(fG[], F_ERR);
#endif
  }

  diffusion (TG, dt, D=lambda2f, r=sgT, theta=thetacorr2);
  diffusion (TL, dt, D=lambda1f, r=slT, theta=thetacorr1);

  foreach() {
    TL[] *= fL[];
    TG[] *= (1. - fL[]);
  }

  /**
  We reconstruct the one-field temperature field summing
  the two fields $T_L$ and $T_G$ in tracer form. */

  foreach() {
    TL[] = (fL[] > F_ERR) ? TL[] : 0.;
    TG[] = (fG[] > F_ERR) ? TG[] : 0.;
    T[] = TL[] + TG[];
  }
}

