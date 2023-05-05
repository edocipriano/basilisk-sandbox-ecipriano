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

#ifndef PHASECHANGE
scalar mEvap[];
scalar * mEvapList = {mEvap};

scalar * gas_tracers = NULL;
scalar * liq_tracers = NULL;
scalar * ftracersorig = NULL;

scalar fL[], fG[], f0[];
face vector fsL[], fsG[];

face vector s[];
# define PHASECHANGE
#endif

scalar T[], TL[], TG[], TInt[];

extern double lambda1, lambda2, dhev, cp1, cp2;
extern double TIntVal;

/**
Allocating diffusivity fields *lambdaf*, volume correction
*theta* and explicit source terms *sgT* and *sTS* for
the diffusion equation. */

face vector lambda1f[], lambda2f[];
scalar thetacorr1[], thetacorr2[];
scalar sgT[], slT[], sgTimp[], slTimp[];

event defaults (i = 0,last)
{
  TL.inverse = false;
  TG.inverse = true;

  gas_tracers = NULL;
  liq_tracers = NULL;

  gas_tracers = {TG};
  liq_tracers = {TL};

  ftracersorig = f.tracers;

#ifdef CONSISTENTPHASE1
  //f.tracers = list_concat (f.tracers, liq_tracers);
  fuext.tracers = list_concat (fuext.tracers, liq_tracers);
#else
  //f.tracers = list_concat (f.tracers, liq_tracers);
  fu.tracers = list_concat (fu.tracers, liq_tracers);
#endif
#ifdef CONSISTENTPHASE2
  //f.tracers = list_concat (f.tracers, gas_tracers);
  fuext.tracers = list_concat (fuext.tracers, gas_tracers);
#else
  //f.tracers = list_concat (f.tracers, gas_tracers);
  fu.tracers = list_concat (fu.tracers, gas_tracers);
#endif
}

event init (i = 0,last)
{
  sgT.nodump = true;
  slT.nodump = true;
  fL.nodump  = true;
  fG.nodump  = true;
  thetacorr1.nodump = true;
  thetacorr2.nodump = true;
}

event finalise (t = end)
{
  gas_tracers    = NULL;
  liq_tracers    = NULL;
  fu.tracers     = NULL;
  fuext.tracers  = NULL;
  f.tracers      = ftracersorig;
}

/**
In the *phasechange* event, the vaporization rate is computed
and the diffusion step for the temperature field (in liquid
and gas phase) is solved. */

event phasechange (i++)
{
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
  boundary({fL,fG,TL,TG,f0});

  /**
  We compute the value of volume fraction *f* on the
  cell-faces using a geometric approach (necessary
  for interface gradients and diffusion equations). */

  face_fraction (fL, fsL);
  face_fraction (fG, fsG);

  /**
  We compute the vaporization rate from the interface
  jump condition, obtaining the vaporization rate per
  unit of interface surface, stored in *mEvap*. */

  foreach() {
    mEvap[] = 0.; TInt[] = 0.;
    if (f[] > F_ERR && f[] < 1.-F_ERR) {
      TInt[] = TIntVal;

      double ltrgrad = ebmgrad (point, TL, fL, fG, fsL, fsG, false, TIntVal, false);
      double gtrgrad = ebmgrad (point, TG, fL, fG, fsL, fsG, true,  TIntVal, false);

      mEvap[] += lambda1*ltrgrad/dhev;
      mEvap[] += lambda2*gtrgrad/dhev;

#ifdef SOLVE_LIQONLY
      mEvap[] = lambda1*ltrgrad/dhev;
#endif

#ifdef SOLVE_GASONLY
      mEvap[] = lambda2*gtrgrad/dhev;
#endif
    }
  }
  boundary({mEvap});

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

      double lheatflux = lambda1*ltrgrad;
      double gheatflux = lambda2*gtrgrad;
      //double lheatflux = mEvapTot[]*dhev;
      //double gheatflux = mEvapTot[]*dhev;

#ifdef AXI
      slT[] = lheatflux/rho1/cp1*area*(y + p.y*Delta)/(Delta*y)*cm[];
      sgT[] = gheatflux/rho2/cp2*area*(y + p.y*Delta)/(Delta*y)*cm[];
      //slT[] = mEvap[]*dhev/rho1/cp1*area*(y + p.y*Delta)/(Delta*y)*cm[];
      //sgT[] = mEvap[]*dhev/rho2/cp2*area*(y + p.y*Delta)/(Delta*y)*cm[];
#else
      slT[] = lheatflux/rho1/cp1*area/Delta*cm[];
      sgT[] = gheatflux/rho2/cp2*area/Delta*cm[];
      //slT[] = mEvap[]*dhev/rho1/cp1*area/Delta*cm[]*f[];
      //sgT[] = mEvap[]*dhev/rho2/cp2*area/Delta*cm[]*(1. - f[]);
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
  boundary({TL,TG,T});
}

event tracer_advection (i++)
{
  /**
  We let the volume fractions *fu* and *fuext* to
  advect the fields TL and TG, as implemented in
  the tracer_advection event of [evaporation.h](evaporation.h) */
}

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
    //fL[] = f[]; fG[] = 1. - f[];
    //TL[] = fL[] > F_ERR ? TL[]/fL[] : 0.;
    //TG[] = ((1. - fL[]) > F_ERR) ? TG[]/(1. - fL[]) : 0.;

#ifdef CONSISTENTPHASE1
    TL[] = fuext[] > F_ERR ? TL[]/fuext[] : 0.;
#else
    TL[] = fu[] > F_ERR ? TL[]/fu[] : 0.;
#endif
#ifdef CONSISTENTPHASE2
    TG[] = ((1. - fuext[]) > F_ERR) ? TG[]/(1. - fuext[]) : 0.;
#else
    TG[] = ((1. - fu[]) > F_ERR) ? TG[]/(1. - fu[]) : 0.;
#endif
    fL[] = f[]; fG[] = 1. - f[];
  }
  boundary({fL,fG,TL,TG});

  /**
  We compute the value of volume fraction *f* on the
  cell-faces using a geometric approach (necessary
  for interface gradients and diffusion equations). */

  face_fraction (fL, fsL);
  face_fraction (fG, fsG);

#ifdef MODIFIED_DIFFUSION
  vector pcs1[], pcs2[];

  foreach() {
    if (fL[] > F_ERR && fL[] < 1.-F_ERR) {
      coord m = mycs (point, fL); 
      double alpha = plane_alpha (fL[], m); 
      coord prel;
      plane_area_center (m, alpha, &prel);
      coord pc; 
      plane_center (m, alpha, fL[], &pc);

      pcs1.x[] = x + pc.x*Delta;
      pcs1.y[] = y + pc.y*Delta;
      pcs1.z[] = z + pc.z*Delta;
    }
    else {
      pcs1.x[] = x;
      pcs1.y[] = y;
      pcs1.z[] = z;
    }

    if (fG[] > F_ERR && fG[] < 1.-F_ERR) {
      coord m = mycs (point, fG); 
      double alpha = plane_alpha (fG[], m); 
      coord prel;
      plane_area_center (m, alpha, &prel);
      coord pc; 
      plane_center (m, alpha, fG[], &pc);

      pcs2.x[] = x + pc.x*Delta;
      pcs2.y[] = y + pc.y*Delta;
      pcs2.z[] = z + pc.z*Delta;
    }
    else {
      pcs2.x[] = x;
      pcs2.y[] = y;
      pcs2.z[] = z;
    }
  }

  face vector corrdist1[], corrdist2[];
  foreach_face() {
    corrdist1.x[] = Delta/fabs (pcs1.x[] - pcs1.x[-1]);
    corrdist2.x[] = Delta/fabs (pcs2.x[] - pcs2.x[-1]);

    // check: Palmore version, I don't think it's correct

    //double magdist1 = 0., magdist2 = 0.;
    //foreach_dimension() {
    //  magdist1 += sq (pcs1.x[] - pcs1.x[-1]);
    //  magdist2 += sq (pcs2.x[] - pcs2.x[-1]);
    //}
    //magdist1 = sqrt (magdist1);
    //magdist2 = sqrt (magdist2);

    //corrdist1.x[] = Delta*(pcs1.x[] - pcs1.x[-1])/magdist1;
    //corrdist2.x[] = Delta*(pcs2.x[] - pcs2.x[-1])/magdist2;
  }
  boundary({corrdist2});
#endif

  /**
  We solve the diffusion equations, confined by means of
  the face fraction fields *fsL* and *fsG*. */

  foreach_face() {
    lambda1f.x[] = lambda1/rho1/cp1*fsL.x[]*fm.x[];
    lambda2f.x[] = lambda2/rho2/cp2*fsG.x[]*fm.x[];
#ifdef MODIFIED_DIFFUSION
    lambda1f.x[] *= corrdist1.x[];  
    lambda2f.x[] *= corrdist2.x[];
#endif
  }
  boundary((scalar *){lambda1f,lambda2f});

  foreach() {
    thetacorr1[] = cm[]*max(fL[], F_ERR);
    thetacorr2[] = cm[]*max(fG[], F_ERR);
  }
  boundary({thetacorr1,thetacorr2});

  foreach() {
    //slTimp[] = -(f[] - f0[])/dt;
    slTimp[] = (f[] - f0[])/dt;
    sgTimp[] = -(f0[] - f[])/dt;
#ifdef AXI
    slTimp[] *= y;
    sgTimp[] *= y;
#endif
  }
  boundary({slTimp,sgTimp});

#ifndef SOLVE_LIQONLY
  //diffusion (TG, dt, D=lambda2f, r=sgT, beta=sgTimp, theta=thetacorr2);
  diffusion (TG, dt, D=lambda2f, r=sgT, theta=thetacorr2);
#endif
#ifndef SOLVE_GASONLY
  //diffusion (TL, dt, D=lambda1f, r=slT, beta=slTimp, theta=thetacorr1);
  diffusion (TL, dt, D=lambda1f, r=slT, theta=thetacorr1);
#endif

  foreach() {
    TL[] *= fL[];
    TG[] *= (1. - fL[]);
  }
  boundary({TL,TG});

  /**
  We reconstruct the one-field temperature field summing
  the two fields $T_L$ and $T_G$ in tracer form. */

  foreach() {
    TL[] = (fL[] > F_ERR) ? TL[] : 0.;
    TG[] = (fG[] > F_ERR) ? TG[] : 0.;
    T[] = TL[] + TG[];
  }
  boundary({TL,TG,T});
}

