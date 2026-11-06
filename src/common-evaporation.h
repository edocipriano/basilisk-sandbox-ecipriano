/**
# Common Evaporation

We define tolerances and functions useful
for evaporation models. */

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

//#include "curvature.h"
#include "adapt_wavelet_leave_interface.h"
#include "mass_refine_prolongation.h"

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
    coord m = interface_normal (point, f);
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

macro number radialprofile (double r,
    double r1, double r2,
    double T1, double T2)
{
  return (r <= r1) ? T1
       : (r >= r2) ? T2
       : (-r1*r2/(r*(r1-r2))*(T1 - T2) + (r1*T1 - r2*T2)/(r1-r2));
}

/**
## Compute the normal in an interfacial cell
*/

coord normal (Point point, scalar c)
{
  coord n = interface_normal (point, c);
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
  coord m = interface_normal (point, f);
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

scalar avg[];

trace
void shift_field (scalar fts, scalar f, int dir) {

  // scalar avg[]; // fixme: if avg is local and not global, there is a problem.
  // My impression is that `no_restriction` is not working - but it should be
  // investigated more deeply
  avg.c = f;
#if TREE
  avg.refine = refinement_avg;
  set_prolongation (avg, refinement_avg);
  set_restriction (avg, no_restriction);
#endif

  // Compute avg
  foreach() {
    avg[] = 0.;
    if (f[] > F_ERR && f[] < 1. - F_ERR) {
      if (dir == 1) {
        int count = 0;
        foreach_neighbor (1) {
          if (f[] > 1.-F_ERR && cm[]) // Number of pure-liquid cells close to the interfacial cell
            count ++;
        }
        avg[] = count;
      }
      else {
        int count = 0;
        foreach_neighbor (1) {
          if (f[] < F_ERR && cm[]) // Number of pure-gas cells close to the interfacial cell
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
      if (f[] > 1.-F_ERR && cm[]) { // Move toward pure-liquid
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
      if (f[] < F_ERR && cm[]) { // Move toward pure-gas
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

#include "diffusion.h"

trace
void shift_diffusion (scalar fts, scalar f,
    double cells = 0.1, double tol = 1e-3)
{
  double smax = 0.;
  foreach (reduction(max:smax))
    if (fabs (fts[]) > smax)
      smax = fabs (fts[]);
  if (smax == 0.)
    return;

  foreach()
    fts[] = (cm[] > 0) ? fts[]/cm[] : 0.;
  double totold = statsf (fts).sum;

  double delta = L0/(1 << grid->maxdepth);
  double diff = cells*sq(delta);

  scalar theta[];
  foreach()
    theta[] = cm[];

  face vector D[];
  foreach_face() {
    D.x[] = diff*fm.x[];
    if (f[-1] > F_ERR && f[-1] < 1.-F_ERR && fts[] != 0.)
      D.x[] *= 0.;
    else if (f[] > F_ERR && f[] < 1.-F_ERR && fts[-1] != 0.)
      D.x[] *= 0.;
  }

  double TOLERANCE_BACKUP = TOLERANCE;
  TOLERANCE = tol*smax;
  diffusion (fts, 1., D = D, theta = theta);
  TOLERANCE = TOLERANCE_BACKUP;

  foreach()
    if (fabs (fts[]*cm[]) < tol*smax)
      fts[] = 0.;
  double totnew = statsf (fts).sum;
  foreach()
    fts[] *= totnew ? cm[]*totold/totnew : 0.;
}

/**
## *copy_bcs()*: Copy the boundary conditions from a target field to a list of fields

* *fts*: field to shift
* *f*: vof volume fraction field
* *dir*: shifting direction: 1 liquid, gas otherwise
*/

void copy_bcs (scalar * dest, scalar orig) {
  for (int bid = 0; bid < nboundary; bid++) {
    for (scalar s in dest) {
      s.boundary[bid] = orig.boundary[bid];
      s.boundary_homogeneous[bid] = orig.boundary_homogeneous[bid];
    }
  }
}

#if TREE
void density_refine (Point point, scalar rhov) {
  refine_bilinear (point, rhov);
  double rhou = 0.;
  foreach_child()
    rhou += cm[]*rhov[];
  double drho = rhov[] - rhou/((1 << dimension)*(cm[] + SEPS));
  foreach_child()
    rhov[] += drho;
}

void density_restriction (Point point, scalar rhov) {
  double rhou = 0.;
  foreach_child()
    rhou += cm[]*rhov[];
  rhov[] = rhou/((1 << dimension)*(cm[] + SEPS));
  //restriction_volume_average (point, rhov);
}

void restriction_scalar_stepwise (Point point, scalar s) {
  scalar c = s.c;
  foreach_blockf (s) {
    double num = 0., den = 0.;
    foreach_child() {
      double w = cm[]*(s.inverse ? 1. - c[] : c[]);
      num += s[]*w;
      den += w;
    }
    s[] = (den > 0.) ? num/den : 0.;
  }
}

foreach_dimension()
static double stepwise_gradient_x (Point point, scalar c, scalar t)
{
  double cl = c[-1], cr = c[1];
  if (t.inverse)
    cl = 1. - cl, cr = 1. - cr;
  if (cl > 0. && cr > 0.)
    return (t[1] - t[-1])/(2.*Delta);
  else if (cr > 0.)
    return (t[1] - t[])/Delta;
  else if (cl > 0.)
    return (t[] - t[-1])/Delta;
  return 0.;
}

void prolongation_scalar_stepwise (Point point, scalar s) {
  scalar f = s.c;
  coord g;
  foreach_dimension()
    g.x = Delta*stepwise_gradient_x (point, f, s);
  double sc = s[], cmc = 4.*cm[];
  foreach_child() {
    s[] = sc;
    if (cmc > 0.)
      foreach_dimension()
        s[] += child.x*g.x*cm[-child.x]/cmc;
  }
}
#endif
