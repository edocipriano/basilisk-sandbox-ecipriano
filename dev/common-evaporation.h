/**
# Common Evaporation

We define tolerances and functions useful
for evaporation models. */

#include "curvature.h"
#include "adapt_wavelet_leave_interface.h"

#ifndef F_ERR
# define F_ERR 1.e-10
#endif
#ifndef P_ERR
# define P_ERR 1.e-10
#endif
#ifndef T_ERR
# define T_ERR 0.
#endif

#ifndef I_TOL
# define I_TOL 1.e-8
#endif

macro foreach_interfacial (scalar f, double tol = 1e-10,
    Reduce reductions = None)
{
  foreach(0, reductions) {
    POINT_VARIABLES();
    if (f[] > tol && f[] < 1. - tol)
    {...}
  }
}

macro foreach_interfacial_plic (scalar f, double tol = 1e-10,
    Reduce reductions = None)
{
  foreach_interfacial (f, tol, reductions) {
    coord m = mycs (point, f);
    double alpha = plane_alpha (f[], m);
    coord prel;
    double area = plane_area_center (m, alpha, &prel);
#if AXI
    double dirac = area*(y + prel.y*Delta)/(Delta*y)*cm[];
#else
    double dirac = area/Delta*cm[];
#endif
    NOT_UNUSED (dirac);
    {...}
  }
}

/**
We define a macro describing a radial
profile. It can be useful for field initialization. */

#define radialprofile(r,r1,r2,T1,T2)(-r1*r2/(r*(r1-r2))*(T1 - T2) + (r1*T1 - r2*T2)/(r1-r2))

/**
## Compute the normal in an interfacial cell
*/

coord normal (Point point, scalar c)
{
  coord n = mycs (point, c);
  double nn = 0.;
  foreach_dimension()
    nn += sq(n.x);
  nn = sqrt(nn);
  foreach_dimension()
    n.x /= nn;
  return n;
}

/**
## *avg_neighbor()*: Compute the average value of a scalar field Y in a 3x3 stencil around the current cell
* *point*: current cell location
* *Y*: field to average
* *f*: vof volume fraction field
*/

double avg_neighbor (Point point, scalar Y, scalar f) {
  double fYnei = 0., fnei = 0.;
  foreach_neighbor(1) {
    double ff = Y.inverse ? 1. - f[] : f[];
    fYnei += ff*Y[];
    fnei  += ff;
  }
  return fYnei / fnei;
}

/**
## *avg_interface()*: Compute the average value of a scalar field Y in a region around the interface
* *Y*: fields to average
* *f*: vof volume fraction field
*/

double avg_interface (scalar Y, scalar f, double tol=F_ERR) {
  double Yavg = 0.;
  int counter = 0;
  foreach(reduction(+:Yavg) reduction(+:counter)) {
    if (f[] > tol && f[] < 1.-tol) {
      counter++;
      Yavg += Y[];
    }
  }
  return Yavg / counter;
}

/**
## *vof_reconstruction()*: VOF reconstruction step

The VOF reconstruction step must be frequently performed to compute
the Dirac delta which allows surface integrals to be transformed
into volume intergrals. Used for the evaporation source terms.

* *point*: current cell in a *foreach()* loop
* *f*: vof field
*/

typedef struct {
  coord m, prel;
  double alpha, area, dirac;
} vofrecon;

vofrecon vof_reconstruction (Point point, scalar f) {
  coord m = mycs (point, f);
  double alpha = plane_alpha (f[], m);
  coord prel;
  double area = plane_area_center (m, alpha, &prel);
#if AXI
  double dirac = area*(y + prel.y*Delta)/(Delta*y)*cm[];
#else
  double dirac = area/Delta*cm[];
#endif

  return (vofrecon){m, prel, alpha, area, dirac};
}

/**
## *shift_field()*: Shift a field localized at the interface toward the closest pure gas or liquid cells

* *fts*: field to shift
* *f*: vof volume fraction field
* *dir*: shifting direction: 1 liquid, gas otherwise
*/

void shift_field (scalar fts, scalar f, int dir) {

  scalar avg[];
//#if TREE
//  avg.refine = avg.prolongation = refinement_avg;
//  avg.restriction = no_restriction;
//  avg.dirty = true;
//#endif

  // Compute avg
  foreach() {
    avg[] = 0.;
    if (f[] > F_ERR && f[] < 1. - F_ERR) {
      if (dir == 1) {
        int count = 0;
        foreach_neighbor (1) {
          if (f[] > 1.-F_ERR) // Number of pure-liquid cells close to the interfacial cell
            count ++;
        }
        avg[] = count;
      }
      else {
        int count = 0;
        foreach_neighbor (1) {
          if (f[] < F_ERR) // Number of pure-gas cells close to the interfacial cell
            count ++;
        }
        avg[] = count;
      }
    }
  }

  scalar sf0[];
  foreach() {
    sf0[] = fts[];
    if (f[] > F_ERR && f[] < 1. - F_ERR)
      fts[] = 0.;
  }

  // Compute m
  foreach() {
    if (dir == 1) {
      if (f[] > 1.-F_ERR) { // Move toward pure-liquid
        double val = 0.;
        foreach_neighbor (1) {
          if (f[] > F_ERR && f[] < 1. - F_ERR && avg[] > 0) {
            val += sf0[]/avg[];
          }
        }
        fts[] += val;
      }
    }
    else {
      if (f[] < F_ERR) { // Move toward pure-gas
        double val = 0.;
        foreach_neighbor (1) {
          if (f[] > F_ERR && f[] < 1. - F_ERR && avg[] > 0) {
            val += sf0[]/avg[];
          }
        }
        fts[] += val;
      }
    }
  }
}

/**
## *distribute_field()*: Distribute a field localized at the interface in the closest pure gas or liquid cells (5x5 stencil)
*/

void distribute_field (scalar ftd, scalar f, int dir) {
  scalar mdisttot[], df[];
  vector mf[];
  foreach() {
    mdisttot[] = 0.;
    if (f[] > F_ERR && f[] < 1.-F_ERR) {
      coord m = mycs (point, f);
      coord n = normal (point, f);
      double alpha = plane_alpha (f[], m);
      coord prel;
      plane_area_center (m, alpha, &prel);
      foreach_dimension()
        mf.x[] = m.x;

      double mdisttothere = 0.;
      foreach_neighbor (2) {
        if ( (f[] > 1.-F_ERR && dir == 1) || (f[] < F_ERR && dir != 1 ) ) {
          coord dist;
          foreach_dimension()
            dist.x = prel.x - x;
          double magdist = 0.;
          foreach_dimension()
            magdist += sq (dist.x);
          magdist = sqrt (magdist);
          double epsi = 0.;
          foreach_dimension()
            epsi += fabs (dist.x*n.x);
          mdisttothere += epsi/magdist;
        }
      }
      mdisttot[] = mdisttothere;
    }
  }

  foreach() {
    df[] = 0.;
    if ( (f[] > 1.-F_ERR && dir == 1) || (f[] < F_ERR && dir != 1 ) ) {
      double dfhere = 0.;
      foreach_neighbor (2) {
        if (f[] > F_ERR && f[] < 1.-F_ERR) {
          coord m, n;
          foreach_dimension() {
            m.x = mf.x[];
            n.x = mf.x[];
          }
          double nnorm = 0.;
          foreach_dimension()
            nnorm += sq (n.x);
          nnorm = sqrt (nnorm);
          foreach_dimension()
            n.x /= nnorm;
          double alpha = plane_alpha (f[], m);
          coord prel;
          plane_area_center (m, alpha, &prel);

          coord dist;
          foreach_dimension()
            dist.x = prel.x - x;
          double magdist = 0.;
          foreach_dimension()
            magdist += sq (dist.x);
          magdist = sqrt (magdist);
          double epsi = 0.;
          foreach_dimension()
            epsi += fabs (dist.x*n.x);
          double mdist = epsi/magdist;

          dfhere += ftd[]*mdist/mdisttot[];
        }
      }
      df[] += dfhere;
    }
  }

  foreach()
    ftd[] = df[];
}

