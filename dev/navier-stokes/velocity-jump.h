#define VELOCITY_JUMP 1

#include "common.h"
#include "utils.h"
#include "aslam.h"

scalar jump[];
face vector jumpf[];
vector n[], ur[];
face vector nf[];
scalar pg[];
bool _boiling = false;

extern scalar f;
extern double mEvapVal, rho1, rho2;

vector ughostl[], ughostg[];
face vector ufghostl[], ufghostg[];

static void setjump (
    scalar f,
    vector ul,
    vector ug,
    scalar jump,
    vector n
)
{
  foreach() {
    if (f[] < 0.5)
      foreach_dimension()
        ul.x[] = ughostl.x[];
    else
      foreach_dimension()
        ug.x[] = ughostg.x[];
  }
}

static void setjumpf (
    scalar f,
    face vector ulf,
    face vector ugf,
    face vector jumpf,
    face vector nf
)
{
  foreach_face() {
    double ff = face_value (f, 0);
    if (ff < 0.5)
      ulf.x[] = ufghostl.x[];
    else
      ugf.x[] = ufghostg.x[];
  }
}

/**
## Projection Function

We define the function that performs the projection step with the
volume expansion term due to the phase change or due to density
changes. */

#include "poisson.h"

extern scalar * drhodtlist, * plist;
extern face vector * uflist;

trace
mgstats project_jump (face vector * uflist, scalar * plist,
     (const) face vector alpha = unityf,
     double dt = 1.,
     int nrelax = 4)
{
  scalar p1 = plist[1], p2 = plist[0];
  face vector uf1 = uflist[1], uf2 = uflist[0];

  /**
  We allocate a local scalar field and compute the divergence of
  $\mathbf{u}_f$. The divergence is scaled by *dt* so that the
  pressure has the correct dimension. */

  scalar div[];
  scalar drhodt1 = drhodtlist[1], drhodt2 = drhodtlist[0];
  foreach() {
    double div1 = 0., div2 = 0.;
    foreach_dimension() {
      div1 += uf1.x[1] - uf1.x[];
      div2 += uf2.x[1] - uf2.x[];
    }
    div1 /= dt*Delta;
    div2 /= dt*Delta;
    div1 += drhodt1[]/dt;
    div2 += drhodt2[]/dt;
    div[] = (f[] > 0.5) ? div1 : div2;
  }

  /**
  We solve the Poisson problem. The tolerance (set with *TOLERANCE*) is
  the maximum relative change in volume of a cell (due to the divergence
  of the flow) during one timestep i.e. the non-dimensional quantity 
  $$
  |\nabla\cdot\mathbf{u}_f|\Delta t 
  $$ 
  Given the scaling of the divergence above, this gives */

  mgstats mgp = poisson (p1, div, alpha,
       tolerance = TOLERANCE/sq(dt), nrelax = nrelax);

  foreach()
    p2[] = p1[];

  /**
  And compute $\mathbf{u}_f^{n+1}$ using $\mathbf{u}_f$ and $p$. */

  foreach_face() {
    uf1.x[] -= dt*alpha.x[]*face_gradient_x (p1, 0);
    uf2.x[] -= dt*alpha.x[]*face_gradient_x (p2, 0);
  }

  return mgp;
}


trace
mgstats project_liquid (face vector uf, scalar p,
     (const) face vector alpha = unityf,
     double dt = 1.,
     int nrelax = 4)
{
  
  /**
  We allocate a local scalar field and compute the divergence of
  $\mathbf{u}_f$. The divergence is scaled by *dt* so that the
  pressure has the correct dimension. */

  scalar div[], drhodt = drhodtlist[1];
  foreach() {
    div[] = 0.;
    foreach_dimension()
      div[] += uf.x[1] - uf.x[];
    div[] /= dt*Delta;
    div[] += drhodt[]/dt;
  }

  /**
  We solve the Poisson problem. The tolerance (set with *TOLERANCE*) is
  the maximum relative change in volume of a cell (due to the divergence
  of the flow) during one timestep i.e. the non-dimensional quantity 
  $$
  |\nabla\cdot\mathbf{u}_f|\Delta t 
  $$ 
  Given the scaling of the divergence above, this gives */

  mgstats mgp = poisson (p, div, alpha,
       tolerance = TOLERANCE/sq(dt), nrelax = nrelax);

  /**
  And compute $\mathbf{u}_f^{n+1}$ using $\mathbf{u}_f$ and $p$. */

  foreach_face()
    uf.x[] -= dt*alpha.x[]*face_gradient_x (p, 0);

  return mgp;
}

#define project(...) project_jump(__VA_ARGS__)
#include "navier-stokes/low-mach.h"
#undef project

event defaults (i = 0) {
  for (int d = 0; d < nboundary; d++) {
    pg.boundary[d] = p.boundary[d];
    pg.boundary_homogeneous[d] = p.boundary_homogeneous[d];
  }
}

event defaults (i = 0) {
  nv = 2;
}

scalar ls[];
extern scalar d;

event advection_term (i++,last) {
  foreach()
    ls[] = d[];

  gradients ({ls}, {n});
  foreach() {
    double mag = 0.;
    foreach_dimension()
      mag += sq (n.x[]);
    mag = sqrt (mag);
    foreach_dimension()
      n.x[] /= -(mag + 1e-10);
  }

  foreach_face()
    nf.x[] = face_value (n.x, 0);

  setjump (f, ulist[1], ulist[0], jump, n);
  setjumpf (f, uflist[1], uflist[0], jumpf, nf);
}

event viscous_term (i++,last) {
  setjump (f, ulist[1], ulist[0], jump, n);
  setjumpf (f, uflist[1], uflist[0], jumpf, nf);
}

event acceleration (i++,last) {
  setjump (f, ulist[1], ulist[0], jump, n);
  setjumpf (f, uflist[1], uflist[0], jumpf, nf);
}

event projection (i++,last) {
  if (_boiling) {
    face vector ufg = uflist[0];
    foreach_face()
      ufghostg.x[] = ufg.x[];
    project_liquid (ufghostg, pg, alpha, dt);
  } else {
    face vector ufl = uflist[1];
    foreach_face()
      ufghostl.x[] = ufl.x[];
    project_liquid (ufghostl, pg, alpha, dt);
  }
}

event end_timestep (i++,last) {
  if (_boiling) {
    face vector ufg = uflist[0];
    foreach_face()
      ufghostl.x[] = ufg.x[] + jumpf.x[]*nf.x[]*fm.x[];

    vector ug = ulist[0];
    foreach()
      foreach_dimension() {
        ughostl.x[] = ug.x[] + jump[]*n.x[];
        ughostg.x[] = 0.5*(ufghostg.x[1] + ufghostg.x[]);
      }
  }
  else {
    face vector ufw[], ufl = uflist[1];
    foreach_face() {
      double ff = face_value (f, 0);
      ufw.x[] = (ff > 0.5) ? ufl.x[] : ufghostl.x[];
    }
    project_liquid (ufw, pg, alpha, dt);

    foreach_face()
      ufl.x[] = ufw.x[];

    foreach_face()
      ufghostg.x[] = ufl.x[] - jumpf.x[]*nf.x[]*fm.x[];

    vector ul = ulist[1];
    foreach()
      foreach_dimension() {
        ughostl.x[] = 0.5*(ufghostl.x[1] + ufghostl.x[]);
        ughostg.x[] = ul.x[] - jump[]*n.x[];
      }
  }

  vector u1 = ulist[1], u2 = ulist[0];
  foreach()
    foreach_dimension()
      ur.x[] = u1.x[]*f[] + u2.x[]*(1. - f[]);
}

