/**
# Main Evaporation Module

Include this file every time you need to perform simulations
with phase change. It contains the declaration of the fields
that are useful to every phase change model, it computes the
total evaporation rate and the phase change velocity. This
file must be used in conjuction with a phase change model (i.e.
multicomponent.h temperature-gradient.h, species-gradient.h). */

#include "common-evaporation.h"
#include "fracface.h"

/**
If DIFFUSIVE conditions are defined, the scalar fields are
advected using the VOF-consistent method as f.tracers. In
this case $stefanflow$ is null. If DIFFUSIVE conditions are
not defined, then the default configuration is used, where
$stefanflow$ is computed, the liquid phase is advected using
the vof-consistent method, while the gas phase is advected
by the [tracer-noncons.h] module. Define NOCONSISTENTPHASE
to suppress the use of VOF-consistent advection for the
scalar fields, and use CONSISTENTPHASE2 to apply the VOF-
consistent advection to the gas phase instead of the liquid. */

#ifndef DIFFUSIVE
# ifndef CONSISTENTPHASE2
#   ifndef NOCONSISTENTPHASE
#     define CONSISTENTPHASE1
#   endif
# endif
#else
# define CONSISTENTPHASE1
# define CONSISTENTPHASE2
#endif

/**
By default, the *stefanflow* is shifted toward the liquid-phase
cells close to the interface. If **SHIFT_TO_GAS** is defined, then
the *stefanflow* field will be shifted toward the gas-phase cells. */

#define SHIFTING
#ifndef SHIFT_TO_GAS
# define SHIFT_TO_LIQ
#endif

/**
It **NOSHIFTING** is defined, no *stefanflow* shifting procedure
is applied. */

#ifdef NOSHIFTING
# undef SHIFTING
#endif

scalar mEvapTot[];      // [kg/m2/s]
scalar stefanflow[];
face vector vpc[];
scalar fu[], fuext[];

// Helper fields
face vector uf_save[];

// Externally defined variables
extern scalar * mEvapList;
extern double rho1, rho2;
extern face vector ufext;
extern scalar rho;

/**
Dynamically change shrinking options. */

bool is_shrinking;

event init (i = 0)
{
  is_shrinking = true;
  
  foreach() {
    fu[] = f[];
    fuext[] = f[];
  }
  boundary({fu,fuext});
  
}

/**
Main phase-change event. Computes the total vaporization
flowrate [kg/m2/s] and interpolates the cell-centered normals
and mEvap in order to compute a face change velocity $vpc$,
defined as a face vector. */

event phasechange (i++)
{
  event ("source_mEvap");

  /**
  We store the value of *f* in the *fu* and *fuext* fields
  at the beginning of the time step. */
  
  foreach() {
    fu[] = f[];
    fuext[] = f[];
  }
  boundary({fu,fuext});
  
  /**
  Update vaporization mass-flowrate. */

  foreach() {
    mEvapTot[] = 0.;
    for (scalar mEvap in mEvapList)
      mEvapTot[] += mEvap[];
  }

  /**
  Update phase-change regression velocity. */

  foreach_face() {
    vpc.x[] = 0.;
    coord nf;
    foreach_dimension()
      nf.x = 0.;
    double mEvapf = 0.;
    if (f[] > 0. && f[] < 1.) {
      nf = normal (point, f);
      mEvapf = mEvapTot[];
      if (f[-1] > 0. && f[-1] < 1.) {
        coord np = normal (neighborp(-1), f);
        foreach_dimension()
          nf.x = (nf.x + np.x)/2.;
        mEvapf = (mEvapf + mEvapTot[-1])/2.;
      }
    }
    else if (f[-1] > 0. && f[-1] < 1.) {
      mEvapf = mEvapTot[-1];
      nf = normal (neighborp(-1), f);
    }
#ifdef BYRHOGAS
    vpc.x[] = fm.x[]*(-mEvapf/rho2)*nf.x;
#else
    vpc.x[] = fm.x[]*(-mEvapf/rho1)*nf.x;
#endif
  }

  /**
  Update stefan flow source term for continuity equation. */

  foreach() {
    stefanflow[] = 0.;
#ifndef DIFFUSIVE
    double alpha = 0.;
    double segment = 0.;
    if (f[] > F_ERR && f[] < 1.-F_ERR) {
      coord m = mycs (point, f);
      alpha = plane_alpha (f[], m);
      coord prel;
      segment = plane_area_center (m, alpha, &prel);
#ifdef AXI
      stefanflow[] = cm[]*mEvapTot[]*segment*(y + prel.y*Delta)*(1./rho2 - 1./rho1)/(Delta*y);
#else
      stefanflow[] = mEvapTot[]*segment*(1./rho2 - 1./rho1)/Delta*cm[];
#endif
    }
#endif
  }
}

#ifdef SHIFTING
event shiftsf (i++) {

  scalar avg[];

#if TREE
  avg.refine = avg.prolongation = refinement_avg;
  avg.restriction = no_restriction;
  avg.dirty = true;
#endif

  // Compute avg
  foreach() {
    avg[] = 0.;
    if (f[] > F_ERR && f[] < 1. - F_ERR) {
      int count = 0;
      foreach_neighbor(1) {
#ifndef SHIFT_TO_LIQ
        if (f[] < F_ERR) // Number of pure-gas cells close to the interfacial cell
#else
        if (f[] > 1.-F_ERR) // Number of pure-liquid cells close to the interfacial cell
#endif
          count ++;
      }
      avg[] = count;
    }
  }

  scalar sf0[];
  foreach() {
    sf0[] = stefanflow[];
    stefanflow[] = 0.;
  }
  boundary({sf0,stefanflow});

  // Compute m
  foreach() {
#ifndef SHIFT_TO_LIQ
    if (f[] < F_ERR) { // Move toward pure-gas
#else
    if (f[] > 1.-F_ERR) { // Move toward pure-liquid
#endif
      double val = 0.;
      foreach_neighbor(1) {
        if (f[] > F_ERR && f[] < 1. - F_ERR && avg[] > 0) {
          val += sf0[]/avg[];
        }
      }
      stefanflow[] = val;
    }
    else {
      stefanflow[] = 0.;
    }
  }
}
#endif

/**
Overloading VOF event in order to include the phase-change velocity. */

static scalar * interfaces1 = NULL;

event vof (i++)
{
  foreach_face() {
    uf_save.x[] = uf.x[];
#ifndef INT_USE_UF
    uf.x[] = ufext.x[];
#endif
#ifndef CONSTANT_VOLUME
    if (is_shrinking)
      uf.x[] -= vpc.x[];
#endif
  }
  boundary((scalar*){uf});
  //event ("stability");

  vof_advection ({f}, i);
  boundary ({f});

  /**
  We restore the value of the $\mathbf{uf}$ velocity field. */

  foreach_face()
    uf.x[] = uf_save.x[];
  boundary((scalar *){uf});

  /**
  We set the list of interfaces to NULL so that the default *vof()*
  event does nothing (otherwise we would transport $f$ twice). */

  interfaces1 = interfaces, interfaces = NULL;
}

event tracer_advection (i++)
{
  /**
  If the phase tracers are not advected with the same velocity
  of the volume fraction field, the advection is performed
  here exploiting the geometric advection procedure but
  using a dummy volume fraction field *fu* and *fuext*. */

  foreach_face() {
    uf_save.x[] = uf.x[];
    uf.x[] = ufext.x[];
  }
  boundary((scalar *){uf});

  vof_advection ({fuext}, i);
  boundary({fuext});

  foreach_face()
    uf.x[] = uf_save.x[];
  boundary({uf});

  /**
  First we store the velocity field $\mathbf{uf}$ which contains
  the Stefan convection and we add the interface regression
  velocity contribution $\mathbf{vpc}$. */

  /**
  We call the vof_advection function that transports the
  volume fraction fu and the associated tracers. */

  vof_advection ({fu}, i);
  boundary({fu});

  /**
  We restore the original list of interfaces, which is
  required by [tension.h]. */

  interfaces = interfaces1;
}

event logfile (i++,last);
event movie (i++,last);

