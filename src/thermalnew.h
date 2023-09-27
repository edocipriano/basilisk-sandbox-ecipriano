/**
# Solver for Energy Equation with Compressibility Effects

This model implements a two-field energy equation in order to resolve
the evolution of the temperature field in a domain with variable
material properties. The flow is not considered incompressible,
instead we introduce a source term for the divergence of the velocity
field, which is a function of the temperature field (pressure
gradients are neglected for the calculation of the divergence).
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
scalar frho1[], frho2[];
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

bool is_diffusion = true;
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
  //double x[] = {1.};
  //ts1.T = 300., ts1.P = Pref, ts1.x = x;
  //ts2.T = 300., ts2.P = Pref, ts2.x = x;

  foreach() {
    rho1v[] = rho1;
    rho2v[] = rho2;
    mu1v[] = mu1;
    mu2v[] = mu2;
    cp1v[] = cp1;
    cp2v[] = cp2;
    lambda1v[] = lambda1;
    lambda2v[] = lambda2;

    frho1[] = f[]*rho1;
    frho2[] = (1. - f[])*rho2;
    frhocp1[] = f[]*rho1*cp1;
    frhocp2[] = (1. - f[])*rho2*cp2;
  }
}

#if TREE
void conservative_refine (Point point, scalar s) {
  double cc = f[];
  double scc = s[];
  if (cc <= 0. || cc >= 1.)
    refine_bilinear(point,s);
  else {

    /**
    Otherwise, we reconstruct the interface in the parent cell. */

    coord n = mycs (point, f);
    double alpha = plane_alpha (cc, n);

    /**
    And compute the volume fraction in the quadrant of the coarse cell
    matching the fine cells. We use symmetries to simplify the
    combinations. */

    coord a, b;
    foreach_dimension() {
      a.x = 0.; b.x = 0.5;
    }

    foreach_child() {
      coord nc; 
      foreach_dimension()
        nc.x = child.x*n.x;
      double crefine = rectangle_fraction (nc, alpha, a, b); 
      if (s.inverse)
        s[] = scc/(1. - cc)*(1. - crefine);
      else
        s[] = scc/cc*crefine;
    }
  }
}
#endif

/**
## Defaults

In the defaults event we setup the tracer lists for the
advection of the temperature fields. */

event defaults (i = 0)
{
  for (scalar s in {TG,frho2})
    s.inverse = true;

  for (scalar s in {frho1,frho2})
    s.conservative = true;

  scalar * ftracersold = list_copy (f.tracers);
  free (f.tracers);
  f.tracers = list_concat (ftracersold, {frho1,frho2,TL,TG});

  for (scalar s in f.tracers)
    s.gradient = f.gradient;

#if TREE
  for (scalar s in {frho1, frho2})
    s.refine = s.prolongation = conservative_refine;
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

event cleanup (t = end) {
  free (f.tracers); f.tracers = NULL;
}

/**
## Phase Change

In the *phasechange* event, the vaporization rate is computed
and the diffusion step for the mass fraction field (in liquid
and gas phase) is solved. */

event phasechange (i++) {
  scalar frho[], betarhocp1[], betarhocp2[];
  foreach() {
    clamp (f[], 0., 1.);
    T[] = TL[] + TG[];
    frho[] = frho1[] + frho2[];
    fL[] = f[]; fG[] = 1. - f[];
    TL[] = (fL[] > F_ERR) ? TL[]/fL[] : 0.;
    TG[] = (fG[] > F_ERR) ? TG[]/fG[] : 0.;
  }

  foreach()
    frho[] = frho1[] + frho2[];

  face_fraction (fL, fsL);
  face_fraction (fG, fsG);

  foreach() {
    double rhocpT1 = rho1v[]*cp1*TL[];
    double rhocpT2 = rho2v[]*cp2*TG[];

    betarhocp1[] = (rhocpT1 != 0.) ? 1./rhocpT1 : 0.;
    betarhocp2[] = (rhocpT2 != 0.) ? 1./rhocpT2 : 0.;
  }

  scalar gradT1[], gradT2[];
  foreach() {
    if (f[] > F_ERR & f[] < 1.-F_ERR) {
      vofrecon vr = vof_reconstruction (point, f);
      double ltrgrad = ebmgrad (point, TL, fL, fG, fsL, fsG, false, TIntVal, false);
      double gtrgrad = ebmgrad (point, TG, fL, fG, fsL, fsG, true,  TIntVal, false);

      gradT1[] = lambda1*ltrgrad*vr.dirac;
      gradT2[] = lambda2*gtrgrad*vr.dirac;
    }
    else {
      gradT1[] = 0.;
      gradT2[] = 0.;
    }
  }

  face vector lambdaT1[], lambdaT2[], lambdaT[];
  foreach_face() {
    lambdaT1.x[] = lambda1*face_gradient_x (TL, 0)*fm.x[]*fsL.x[];
    lambdaT2.x[] = lambda2*face_gradient_x (TG, 0)*fm.x[]*fsG.x[];
    double ff = 0.5*(f[] + f[-1]);
    double lambdaf  = ff*lambda1 + (1. - ff)*lambda2;
    lambdaT.x[] = fm.x[]*lambdaf*face_gradient_x (T, 0);
  }

  foreach() {
    double drhodt1 = 0., drhodt2 = 0.;
    foreach_dimension() {
      drhodt1 += (lambdaT1.x[1] - lambdaT1.x[])/Delta;
      drhodt2 += (lambdaT2.x[1] - lambdaT2.x[])/Delta;
    }
    drhodt1 += gradT1[];
    drhodt2 += gradT2[];

    double rhocpT1 = rho1v[]*cp1*TL[];
    double rhocpT2 = rho2v[]*cp2*TG[];

    drhodt1 *= (rhocpT1 != 0.) ? 1./rhocpT1 : 0.;
    drhodt2 *= (rhocpT2 != 0.) ? 1./rhocpT2 : 0.;

    drhodt[] = aavg (f[], drhodt1, drhodt2);
  }

  foreach() {
    TL[] *= f[];
    TG[] *= (1. - f[]);
  }
}

/**
## Tracer Advection

We let the volume fractions *fu* and *fuext* to
advect the fields YL and YG, as implemented in
the tracer_advection event of [evaporation.h](evaporation.h)
*/

event tracer_advection (i++);

/**
## Tracer Diffusion

We solve the diffusion equations for the temperature fields
accounting for the phase change contributions. */

event tracer_diffusion (i++)
{
  foreach() {
    f[] = clamp (f[], 0., 1.);
    fL[] = f[]; fG[] = 1. - f[];
    // Reconstruct density
    rho1v[] = (fL[] > F_ERR) ? frho1[]/fL[] : 1.e-10;
    rho2v[] = (fG[] > F_ERR) ? frho2[]/fG[] : 1.e-10;
    // Reconstruct temperature
    TL[] = (fL[] > F_ERR) ? TL[]/fL[] : 0.;
    TG[] = (fG[] > F_ERR) ? TG[]/fG[] : 0.;
  }

  if (is_diffusion) {

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
        vofrecon vr = vof_reconstruction (point , f);

        double ltrgrad = ebmgrad (point, TL, fL, fG, fsL, fsG, false, TIntVal, false);
        double gtrgrad = ebmgrad (point, TG, fL, fG, fsL, fsG, true, TIntVal, false);

        double lheatflux = lambda1*ltrgrad;
        double gheatflux = lambda2*gtrgrad;

        slT[] = lheatflux/rho1v[]/cp1*vr.dirac;
        sgT[] = gheatflux/rho2v[]/cp2*vr.dirac;
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
      double rho1f = 0.5*(rho1v[] + rho1v[-1]) + 1.e-10;
      double rho2f = 0.5*(rho2v[] + rho2v[-1]) + 1.e-10;
      lambda1f.x[] = lambda1/rho1f/cp1*fsL.x[]*fm.x[];
      lambda2f.x[] = lambda2/rho2f/cp2*fsG.x[]*fm.x[];
    }

    foreach() {
      thetacorr1[] = cm[]*max (fL[], F_ERR);
      thetacorr2[] = cm[]*max (fG[], F_ERR);
    }

    diffusion (TL, dt, D=lambda1f, r=slT, theta=thetacorr1);
    diffusion (TG, dt, D=lambda2f, r=sgT, theta=thetacorr2);
  }

  foreach() {
    TL[] *= f[];
    TG[] *= (1. - f[]);
    T[] = TL[] + TG[];
  }
}

