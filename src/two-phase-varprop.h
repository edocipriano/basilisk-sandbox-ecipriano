/**
# Two-phase interfacial flows

This file helps setup simulations for flows of two fluids separated by
an interface (i.e. immiscible fluids). It is typically used in
combination with a [Navier--Stokes solver](navier-stokes/centered.h). 

The interface between the fluids is tracked with a Volume-Of-Fluid
method. The volume fraction in fluid 1 is $f=1$ and $f=0$ in fluid
2. The densities and dynamic viscosities for fluid 1 and 2 are *rho1*,
*mu1*, *rho2*, *mu2*, respectively. */

#include "vof-varprop.h"

scalar f[], * interfaces = {f};
double rho1 = 1., mu1 = 0., rho2 = 1., mu2 = 0.;
#ifdef VARPROP
extern scalar rho1v, rho2v, mu1v, mu2v;
#endif

/**
Auxilliary fields are necessary to define the (variable) specific
volume $\alpha=1/\rho$ as well as the cell-centered density. */

face vector alphav[];
scalar rhov[];

event defaults (i = 0) {
  alpha = alphav;
  rho = rhov;

  /**
  If the viscosity is non-zero, we need to allocate the face-centered
  viscosity field. */
  
  if (mu1 || mu2)
    mu = new face vector;

  /**
  We add the interface to the default display. */

  display ("draw_vof (c = 'f');");
}

/**
The density and viscosity are defined using arithmetic averages by
default. The user can overload these definitions to use other types of
averages (i.e. harmonic). */

#define aavg(f,v1,v2)(clamp(f,0.,1.)*(v1 - v2) + v2)
#define havg(f,v1,v2)(1./(clamp(f,0,1)*(1./(v1) - 1./(v2)) + 1./(v2)))

/**
We have the option of using some "smearing" of the density/viscosity
jump. */

#ifdef FILTERED
scalar sf[];
#else
# define sf f
#endif

event tracer_advection (i++)
{
  
  /**
  When using smearing of the density jump, we initialise *sf* with the
  vertex-average of *f*. */

#ifndef sf
#if dimension <= 2
  foreach()
    sf[] = (4.*f[] + 
	    2.*(f[0,1] + f[0,-1] + f[1,0] + f[-1,0]) +
	    f[-1,-1] + f[1,-1] + f[1,1] + f[-1,1])/16.;
#else // dimension == 3
  foreach()
    sf[] = (8.*f[] +
	    4.*(f[-1] + f[1] + f[0,1] + f[0,-1] + f[0,0,1] + f[0,0,-1]) +
	    2.*(f[-1,1] + f[-1,0,1] + f[-1,0,-1] + f[-1,-1] + 
		f[0,1,1] + f[0,1,-1] + f[0,-1,1] + f[0,-1,-1] +
		f[1,1] + f[1,0,1] + f[1,-1] + f[1,0,-1]) +
	    f[1,-1,1] + f[-1,1,1] + f[-1,1,-1] + f[1,1,1] +
	    f[1,1,-1] + f[-1,-1,-1] + f[1,-1,-1] + f[-1,-1,1])/64.;
#endif
#endif // !sf

#if TREE
  sf.prolongation = refine_bilinear;
  sf.dirty = true; // boundary conditions need to be updated
#endif
}

event properties (i++)
{
  foreach_face() {
    double ff = (sf[] + sf[-1])/2.;
#ifdef VARPROP
    double rho1vh = 0.5*(rho1v[] + rho1v[-1]);
    double rho2vh = 0.5*(rho2v[] + rho2v[-1]);
#else
    double rho1vh = rho1;
    double rho2vh = rho2;
#endif
    alphav.x[] = fm.x[]/aavg (ff, rho1vh, rho2vh);
    if (mu1 || mu2) {
      face vector muv = mu;
#ifdef VARPROP
      double mu1vh = 0.5*(mu1v[] + mu1v[-1]);
      double mu2vh = 0.5*(mu2v[] + mu2v[-1]);
#else
      double mu1vh = mu1;
      double mu2vh = mu2;
#endif
      muv.x[] = fm.x[]*aavg (ff, mu1vh, mu2vh);
    }
  }
  
  foreach() {
#ifdef VARPROP
    double rho1vh = rho1v[];
    double rho2vh = rho2v[];
#else
    double rho1vh = rho1;
    double rho2vh = rho2;
#endif
    rhov[] = cm[]*aavg (sf[], rho1vh, rho2vh);
  }

#if TREE  
  sf.prolongation = fraction_refine;
  sf.dirty = true; // boundary conditions need to be updated
#endif
}
