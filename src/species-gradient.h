/**
# Species Gradient Phase Change Model
This phase change model is suitable for vaporization driven
by a gradient of chemical species, between the interface
and the sorrounding gas phase. The vaporization rate (per
unit of interface surface) is computed from a mass balance
at the interface, assuming that the mass fraction on the gas
phase side is constant (i.e. pure liquid phase and constant
vapor pressure). Therefore, the vaporization rate is computed
from the following balance:

$$
\dot{m}_{Evap} = -\dfrac{\rho_g\mathcal{D}_g}{1-\hat{Y}_G}
\left.\dfrac{\partial T}{\partial \mathbf{n}_\Gamma}\right\vert_g
$$

where $\mathcal{D}_G$ is the diffusivity coefficient in the gas phase,
while $\hat{Y}_G$ is the saturation mass fraction on the
gas phase side. If *DIFFUSIVE* conditions are defined, the
stefan convection is neglected from the mass balance, and
the vaporization rate is computed with the same formula
but neglecting the denominator:

$$
\dot{m}_{Evap} = -\rho_g\mathcal{D}_g
\left.\dfrac{\partial T}{\partial \mathbf{n}_\Gamma}\right\vert_g
$$

which results in a vaporization rate proportional to the diffusive
flux in the gas phase.
After the calculation of the vaporization rate, this model
solves the transport equation for the mass fraction field $Y$
using a two-field formulation and including the presence of
the phase change, both in the diffusion equation and in the
transport of the mass fraction field which keeps into account
the Stefan convection.
*/

#include "intgrad.h"
#include "fracface.h"
#include "diffusion.h"

/**
## Memory Allocations
This phase change model defines a list of vaporization
rates *mEvapList* which is populated by a single scalar
*mEvap* because a single chemical species in liquid phase is
considered. */

#ifndef PHASECHANGE
scalar mEvap[];
scalar * mEvapList = {mEvap};

/**
The following list allows different methods for the advection
of the tracers in liquid and gas phase to be selected. The
user should not call these lists. */

scalar * gas_tracers = NULL;
scalar * liq_tracers = NULL;
scalar * ftracersorig = NULL;

/**
*fL* and *fG* store the value of volume fractions for the
calculation of the interface gradients, while the vectors
*fsL* and *fsG* contain the face fraction fields computed
using the [fracface.h](/sandbox/ecipriano/src/fracface.h)
module. */

scalar fL[], fG[], f0[];
face vector fsL[], fsG[];

# define PHASECHANGE
#endif

/**
* *Y* one-field mass fraction of the chemical species to be solved
* *YL* liquid-phase mass fraction
* *YG* gas-phase mass fraction
* *YInt* value of mass fraction at the interface (could be avoided)
*/

scalar Y[], YL[], YG[], YInt[];

/**
## User Data
Using this phase change model, the user should define the
following variables:

* *Dmix1* Diffusivity coefficient in the liquid phase
* *Dmix2* Diffusivity coefficient in the gas phase
* *YIntVal* Value of mass fraction at the interface on the
gas phase side. */

extern double Dmix1, Dmix2;
extern double YIntVal;

/**
Allocating diffusivity fields *Dmixf*, volume correction
*theta* and explicit source terms *sgS* and *slS* for
the diffusion equation. */

face vector Dmix1f[], Dmix2f[];
scalar thetacorr1[], thetacorr2[];
scalar sgS[], slS[];
scalar slimp[], sgimp[];

event defaults (i = 0,last)
{
  YL.inverse = false;
  YG.inverse = true;

  gas_tracers = NULL;
  liq_tracers = NULL;

  gas_tracers = {YG};
  liq_tracers = {YL};

  ftracersorig = f.tracers;

#ifdef CONSISTENTPHASE1
  fuext.tracers = list_concat (fuext.tracers, liq_tracers);
  //fu.tracers = list_concat (fu.tracers, liq_tracers);
#else
  fu.tracers = list_concat (fu.tracers, liq_tracers);
  //fuext.tracers = list_concat (fuext.tracers, liq_tracers);
#endif
#ifdef CONSISTENTPHASE2
  //f.tracers = list_concat (f.tracers, gas_tracers);
  fuext.tracers = list_concat (fuext.tracers, gas_tracers);
#else
  fu.tracers = list_concat (fu.tracers, gas_tracers);
#endif
}

event init (i = 0,last)
{
  sgS.nodump = true;
  slS.nodump = true;
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
## Vaporization Rate and Diffusion Equations
In the *phasechange* event, the vaporization rate is computed
and the diffusion step for the mass fraction field (in liquid
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
    YL[] = f[] > F_ERR ? YL[]/f[] : 0.;
    YG[] = ((1. - f[]) > F_ERR) ? YG[]/(1. - f[]) : 0.;
  }
  boundary({fL,fG,YL,YG,f0});

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
    mEvap[] = 0.; YInt[] = 0.;
    if (f[] > F_ERR && f[] < 1.-F_ERR) {
      YInt[] = YIntVal;

      double gtrgrad = ebmgrad (point, YG, fL, fG, fsL, fsG, true, YIntVal, false);

#ifndef DIFFUSIVE
      mEvap[] = -rho2*Dmix2*gtrgrad/min(1. - YIntVal, 0.99);
#else
      mEvap[] = -rho2*Dmix2*gtrgrad;
#endif
    }
  }
  boundary({mEvap});

  /**
  The source term for the diffusion equation of the
  species mass fractions in gas an liquid phase
  are computed here. */

  foreach() {
    sgS[] = 0., slS[] = 0.;
    sgimp[] = 0., slimp[] = 0.;
    if (f[] > F_ERR && f[] < 1.-F_ERR) {
      coord n = facet_normal (point, fL, fsL), p;
      double alpha = plane_alpha (fL[], n);
      double area = plane_area_center (n, alpha, &p);
      normalize (&n);

#ifdef AXI
      sgS[] = -mEvap[]/rho2*area*(y + p.y*Delta)/(Delta*y)*cm[];
      slS[] =  mEvap[]/rho1*area*(y + p.y*Delta)/(Delta*y)*cm[];
      sgimp[] =  mEvap[]/rho2*area*(y + p.y*Delta)/(Delta*y)*cm[];
      slimp[] = -mEvap[]/rho1*area*(y + p.y*Delta)/(Delta*y)*cm[];
#else
      sgS[] = -mEvap[]/rho2*area/Delta*cm[];
      slS[] =  mEvap[]/rho1*area/Delta*cm[];
      sgimp[] =  mEvap[]/rho2*area/Delta*cm[];
      slimp[] = -mEvap[]/rho1*area/Delta*cm[];
#endif
    }
  }

  /**
  We restore the tracer form of the liquid and gas-phase
  mass fraction fields. */

  foreach() {
    YL[] *= f[]*(f[] > F_ERR);
    YG[] *= (1. - f[])*((1. - f[]) > F_ERR);
    Y[]  = YL[] + YG[];
  }
  boundary({YL,YG,Y});
}

event tracer_advection (i++)
{
  /**
  We let the volume fractions *fu* and *fuext* to
  advect the fields YL and YG, as implemented in
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
    YL[] = fuext[] > F_ERR ? YL[]/fuext[] : 0.;
    YG[] = ((1. - fu[]) > F_ERR) ? YG[]/(1. - fu[]) : 0.;
    fL[] = f[]; fG[] = 1. - f[];
  }
  boundary({fL,fG,YL,YG});

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
    Dmix1f.x[] = Dmix1*fsL.x[]*fm.x[];
    Dmix2f.x[] = Dmix2*fsG.x[]*fm.x[];
  }
  boundary((scalar *){Dmix1f,Dmix2f});

  foreach() {
    thetacorr1[] = cm[]*max(fL[], F_ERR);
    thetacorr2[] = cm[]*max(fG[], F_ERR);
  }
  boundary({thetacorr1,thetacorr2});

#ifndef SOLVE_LIQONLY
  diffusion (YG, dt, D=Dmix2f, r=sgS, beta=sgimp, theta=thetacorr2);
#endif
#ifndef SOLVE_GASONLY
  diffusion (YL, dt, D=Dmix1f, r=slS, beta=slimp, theta=thetacorr1);
#endif

  foreach() {

    /**
    Reconstruct mass fractions computing the mass of each
    chemical species in every cell. */

    double totmassliq = rho1*YL[]*fL[]*dv();
    YL[] = (totmassliq > 0.) ? rho1*YL[]*fL[]*dv() / totmassliq : 0.;

    double totmassgas = rho2*YG[]*fG[]*dv() + rho2*(1. - YG[])*fG[]*dv();
    YG[] = (totmassgas > 0.) ? rho2*YG[]*fG[]*dv() / totmassgas : 0.;

    YL[] *= fL[];
    YG[] *= (1. - fL[]);
  }
  boundary({YL,YG});

  /**
  We reconstruct the one-field mass fraction field summing
  the two fields $Y_L$ and $Y_G$ in tracer form. */

  foreach() {
    YL[] = clamp (YL[], 0., 1.);
    YG[] = clamp (YG[], 0., 1.);
    YL[] = (YL[] > F_ERR) ? YL[] : 0.;
    YG[] = (YG[] > F_ERR) ? YG[] : 0.;
    Y[] = YL[] + YG[];
  }
  boundary({YL,YG,Y});
}

