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
    int i,                      // Iteration index
    bool byrhogas = false,      // Use the gas phase density instead of liquid
    bool extrapolation = false  // Use PDE-based extrapolation
)
{
  vector n[];
  face vector upc[];

  foreach() {
    coord m = mycs (point, f);
    double nn = 0.;
    foreach_dimension()
      nn += sq (m.x);
    foreach_dimension()
      n.x[] = m.x/sqrt (nn);
    foreach_dimension()
      if (byrhogas)
        n.x[] *= (rho2 > 0.) ? mEvapTot[]/rho2 : 0.;
      else
        n.x[] *= (rho1 > 0.) ? mEvapTot[]/rho1 : 0.;
  }

  if (extrapolation) {
#if VELOCITY_JUMP
    scalar fext[], dext[];
    foreach() {
      fext[] = f[];
      dext[] = -d[];
    }

    foreach_dimension() {
      constant_extrapolation (n.x, dext, 0.5, 10, c=fext, nl=0,
          nointerface=true, inverse=false);
      constant_extrapolation (n.x, dext, 0.5, 10, c=fext, nl=0,
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

  vof_advection ({f}, i);

  foreach_face()
    uf.x[] = ufsave.x[];
}

/**
## Explicit source term

We directly introduce the phase change regression velocity as an explicit source
term in the VOF transport equation.
*/

void vof_expl_sources (
    scalar f,               // VOF volume fraction
    scalar mEvapTot,        // Total evaporation rate per unit of surface
    double dt,              // Time step
    bool byrhogas = false   // Use the gas phase density instead of liquid
)
{
  foreach_interfacial_plic (f, F_ERR) {
    if (byrhogas)
      f[] += (rho2 > 0.) ? dt*mEvapTot[]/rho2*area/Delta : 0.;
    else
      f[] += (rho1 > 0.) ? dt*mEvapTot[]/rho1*area/Delta : 0.;
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
    double dt,              // Time step
    bool byrhogas = false   // Use the gas phase density instead of liquid
)
{
  foreach_interfacial_plic (f, F_ERR) {
    double delta_alpha = 0.;
    double magn = sqrt (sq (m.x) + sq (m.y) + sq (m.z));
    if (byrhogas)
      delta_alpha = (rho2 > 0.) ? dt*mEvapTot[]*magn/rho2/Delta : 0.;
    else
      delta_alpha = (rho1 > 0.) ? dt*mEvapTot[]*magn/rho1/Delta : 0.;

    double ff = plane_volume (m, alpha + delta_alpha);
    f[] = (ff > F_ERR && ff < 1.-F_ERR) ? ff : (ff <= F_ERR) ? 0. : 1.;
  }
}

