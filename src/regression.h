/**
# Interface Regression Velocity

This module implements the interface velocity with phase change. In this case,
the advection of the indicator function is obtained from the solution of the
following transport equation:

$$
  \dfrac{DH}{Dt} = \dfrac{\partial H}{\partial t}
  + \mathbf{u}_\Gamma\cdot\nabla H = 0
$$

where the interface velocity is obtained from the Rankine-Hugoniot conidition:

$$
  \mathbf{u}_\Gamma
  = \mathbf{u}_l - \dfrac{\dot{m}}{\rho_l}\mathbf{n}_\Gamma
  = \mathbf{u}_g - \dfrac{\dot{m}}{\rho_g}\mathbf{n}_\Gamma
  = \mathbf{u}_k - \dfrac{\dot{m}}{\rho_k}\mathbf{n}_\Gamma
$$

where $\dot{m}$ is the vaporization rate per unit of surface. Splitting the
liquid and gas phase velocities from the phase change contributions, we obtain:

$$
  \dfrac{\partial H}{\partial t} + \mathbf{u}_k\cdot\nabla H
  = \textcolor{blue}{\dfrac{\dot{m}}{\rho_k}\delta_\Gamma}
$$

The terms on the LHS are aproximated by the [vof](/src/vof.h) module, while in
this module we approximate the last term on the RHS (in blue) either by
computing the interface velocity, or by discretizing the source term.
*/

#include "aslam.h"

/**
## Interface velocity

Apply the phase change velocity $\mathbf{u}_{pc} =
\dot{m}/\rho_k\mathbf{n}_\Gamma$ as an additional velocity which transports the
interface. The discontinuous field $\dot{m}$ is interpolated from cells to edges
in two possible ways: i) using the (fast) VOF-based approach, ii) adopting
PDE-based extrapolation techniques.
*/

void vof_advection_phasechange (
    scalar f,                   // VOF volume fraction
    scalar mEvapTot,            // Total evaporation rate per unit of surface
    scalar rho,                 // Density field
    int i,                      // Iteration index
    bool extrapolation = false  // Use PDE-based extrapolation
)
{
  face vector upc[];

  /**
  The default antisymmetry boundary conditions would compute the wrong flux at
  the boundary. */

  vector n[];
  foreach_dimension()
    for (int b = 0; b < nboundary; b++)
      n.x.boundary[b] = n.x.boundary_homogeneous[b] = symmetry;

  foreach() {
    coord m = interface_normal (point, f);
    normalize (&m);
    foreach_dimension()
      n.x[] = m.x;
    foreach_dimension()
      n.x[] *= (rho[] > 0.) ? mEvapTot[]/rho[] : 0.;
  }

  if (extrapolation) {
#if VELOCITY_JUMP
    scalar fext[];
    foreach()
      fext[] = f[];
    vof_to_ls (fext, ls, imax = 5);

    foreach_dimension() {
      constant_extrapolation (n.x, ls, 0.5, 10, c=fext, nl=0,
          nointerface=true, inverse=false);
      constant_extrapolation (n.x, ls, 0.5, 10, c=fext, nl=0,
          nointerface=true, inverse=true);
    }

    foreach_face()
      upc.x[] = face_value (n.x, 0)*fm.x[];
#endif
  }
  else {
    scalar interf[];
    foreach()
      interf[] = (f[] > 0. && f[] < 1.) ? 1. : 0.;

    foreach_face() {
      upc.x[] = (interf[] && interf[-1]) ? 0.5*(n.x[] + n.x[-1])
        : (interf[]) ? n.x[] : n.x[-1];
      upc.x[] *= fm.x[];
    }
  }

  face vector ufsave[];
  foreach_face() {
    ufsave.x[] = uf.x[];
    uf.x[] = upc.x[];
  }
  boundary ((scalar *){uf});  // do not remove

  /**
  The tracers are transported with first-order upwinding: the liquid removed
  by the phase change leaves the cell with the mean value of the cell, as with
  the source term. */

  int nt = list_len (f.tracers), k = 0;
  double (* gsave[nt + 1]) (double, double, double);
  for (scalar t in f.tracers)
    gsave[k++] = t.gradient, t.gradient = zero;

  vof_advection ({f}, i);

  k = 0;
  for (scalar t in f.tracers)
    t.gradient = gsave[k++];

  foreach_face()
    uf.x[] = ufsave.x[];
  boundary ((scalar *){uf});  // do not remove
}

/**
## Explicit source term

We directly introduce the phase change regression velocity as an explicit source
term in the VOF transport equation.
*/

void vof_expl_sources (
    scalar f,               // VOF volume fraction
    scalar mEvapTot,        // Total evaporation rate per unit of surface
    scalar rho,             // Density field
    double dt               // Time step
)
{
  scalar df[];
  foreach_interfacial_plic (f, F_ERR)
    df[] = (rho[] > 0.) ? dt*mEvapTot[]/rho[]*dirac/cm[] : 0.;

  foreach_interfacial (f, F_ERR) {
    double f0 = f[];
    f[] = clamp (f[] + df[], 0., 1.);

    /**
    The tracers of $f$ are volume-averaged (i.e. $f\phi$ and $(1 - f)\phi$),
    therefore we rescale them with the new volume of their phase: the interface
    motion does not modify the intensive values $\phi$, which are updated with
    the boundary conditions of the diffusion step. */

    for (scalar t in f.tracers)
      t[] *= t.inverse ? (1. - f[])/(1. - f0) : f[]/f0;
  }
}

/**
## Plane shifting

We apply the interface regression by shifting the PLIC approximation of the
interface along its normal direction.
*/

void vof_plane_shifting (
    scalar f,               // VOF volume fraction
    scalar mEvapTot,        // Total evaporation rate per unit of surface
    scalar rho,             // Density field
    double dt               // Time step
)
{
  /**
  The new volume fraction is computed before updating $f$, otherwise the
  PLIC reconstruction would use already updated neighbors. */

  scalar fnew[];
  foreach_interfacial_plic (f, F_ERR) {
    double delta_alpha = 0.;
    double magn = sqrt (sq (m.x) + sq (m.y) + sq (m.z));
    delta_alpha = (rho[] > 0.) ? dt*mEvapTot[]*magn/rho[]/Delta : 0.;
    double ff = plane_volume (m, alpha + delta_alpha);
    fnew[] = (ff > F_ERR && ff < 1.-F_ERR) ? ff : (ff <= F_ERR) ? 0. : 1.;
  }

  foreach_interfacial (f, F_ERR) {
    double f0 = f[];
    f[] = fnew[];

    /**
    The tracers are rescaled with the new volume of their phase, as in
    *vof_expl_sources()*. */

    for (scalar t in f.tracers)
      t[] *= t.inverse ? (1. - f[])/(1. - f0) : f[]/f0;
  }
}

