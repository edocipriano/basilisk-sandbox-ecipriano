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
## Interfacial values

A phase tracer (e.g. the temperature or a mass fraction of a phase) can be
associated with the field of its interfacial value $\phi_\Gamma$ (e.g. the
interface temperature) using the `interfacial` attribute. As for the other
scalar attributes, the default index 0 means that the field is not set. */

attribute {
  scalar interfacial;
}

/**
### *phasechange_tracers()*: tracers update after phase change volume change

This function applies the phase change volume change to the list of tracers.
It is performed separately from the VOF phase change step, where the phase
tracers are not transported with the volume fraction: the VOF fluxes would
carry the value of the tracer in the cell, while the volume converted by phase
change must carry the interfacial value.

To understand the mathematics involved, consider a two-field scalar ($\phi$)
advection in non-conservative form: $H\,D\phi/Dt = 0$, where the material
derivative $D/Dt = \partial/\partial t + \mathbf{u}_k\cdot\nabla$ uses the
velocity of the phase. This is not the form that we resolve using the
VOF-consistent tracers transport, which instead resolves for the extensive
quantity: $D(H\phi)/Dt = 0$.

Manipulating the initial equation to obtain the form that matches the transport
of VOF-tracers, we obtain:

$$
  \dfrac{DH\phi}{Dt} = \phi\dfrac{DH}{Dt} = \phi_\Gamma\dfrac{DH}{Dt}
$$

where $\phi$ is written as $\phi_\Gamma$ because the volume fraction changes are
localized at the interface. With phase change, the interface does not move
with the velocity of the phase, and the material derivative of $H$ is
$DH/Dt = \dot{m}\,\delta_\Gamma/\rho$, which allows us to write:

$$
  \dfrac{DH\phi}{Dt} = \phi_\Gamma\dfrac{\dot{m}}{\rho}\delta_\Gamma
$$

This equation is resolved in two steps, where $H$ is replaced by its cell
average, the volume fraction $c$ of the phase ($c = f$ or $c = 1 - f$ for
`inverse` tracers). First, the phase change step, without transport:

$$
  \dfrac{\partial c\phi}{\partial t} = \phi_\Gamma\dfrac{\partial c}{\partial t}
$$

where the change of $c$ is computed by one of the functions below. Integrating
over the time step, with $\phi_\Gamma$ constant in the cell, gives:

$$
  \Delta (c\phi) = \phi_\Gamma \Delta c
$$

which is essentially Eq. (43) from [Malan et al., 2021](#malan2021). Then,
$c\phi$ is advected with the velocity of the phase, $D(c\phi)/Dt = 0$, by
[vof.h](/src/vof.h). */

void phasechange_tracers (scalar * tracers, scalar f0, scalar f,
    double tol = F_ERR)
{
  foreach()
    for (scalar t in tracers) {
      double ff0 = t.inverse ? 1. - f0[] : f0[];
      double ff = t.inverse ? 1. - f[] : f[];
      if (t.interfacial.i)
        t[] += t.interfacial[]*(ff - ff0);
      // Adjust the tracer without changing its intensive value
      else
        t[] = (ff0 > tol) ? t[]*ff/ff0 : 0.;
    }

  // phase appearing without interfacial value (ff0 <= tol and ff > tol)
  foreach()
    for (scalar t in tracers)
      if (!t.interfacial.i) {
        double ff0 = t.inverse ? 1. - f0[] : f0[];
        double ff = t.inverse ? 1. - f[] : f[];
        if (ff0 <= tol && ff > tol) {
          double tn = 0., ffn = 0.;
          foreach_neighbor (1) {
            double ff0n = t.inverse ? 1. - f0[] : f0[];
            double ff1n = t.inverse ? 1. - f[] : f[];
            if (ff0n > tol && ff1n > tol) {
              tn += t[];
              ffn += ff1n;
            }
          }
          t[] = (ffn > 0.) ? ff*tn/ffn : 0.;
        }
      }

  /**
  Where the phase disappears (ff <= tol) the tracer keeps a residual, which is
  $ff_0(\phi_0 - \phi_\Gamma)$ with the interfacial value, and negligible (at
  most $tol\,\phi_0$) otherwise. It is moved to the neighbors of the same
  phase, proportionally to their volume of the phase. */

  scalar * wl = NULL;
  for (scalar t in tracers) {
    scalar w = new scalar;
    wl = list_append (wl, w);
  }

  foreach() {
    scalar t, w;
    for (t, w in tracers, wl) {
      double ff = t.inverse ? 1. - f[] : f[];
      double wsum = 0.;
      if (ff <= tol && t[])
        foreach_neighbor (1) {
          double ffn = t.inverse ? 1. - f[] : f[];
          wsum += (ffn > tol) ? ffn*dv() : 0.;
        }
      w[] = (wsum > 0.) ? t[]*dv()/wsum : 0.;
    }
  }

  foreach() {
    scalar t, w;
    for (t, w in tracers, wl) {
      double ff = t.inverse ? 1. - f[] : f[];
      if (ff > tol) {
        double add = 0.;
        foreach_neighbor (1)
          add += w[];
        t[] += ff*add;
      }
      else
        t[] = 0.;
    }
  }

  delete (wl), free (wl);
}

/**
## Interface velocity

Apply the phase change velocity $\mathbf{u}_{pc} =
\dot{m}/\rho_k\mathbf{n}_\Gamma$ as an additional velocity which transports the
interface. The discontinuous field $\dot{m}$ is interpolated from cells to edges
in two possible ways: i) using the (fast) VOF-based approach ([Cipriano et al.,
2024](#cipriano2023), Sec. 3.1, Eq. 20 and Fig. 4), ii) adopting PDE-based
extrapolation techniques.
*/

void vof_advection_phasechange (
    scalar f,                   // VOF volume fraction
    scalar mEvapTot,            // Total evaporation rate per unit of surface
    scalar rho,                 // Density field
    int i,                      // Iteration index
    bool extrapolation = false  // Use PDE-based extrapolation
)
{
  vector n[];
  face vector upc[];

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

  vof_advection ({f}, i);

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
  foreach_interfacial_plic (f, F_ERR)
    f[] += (rho[] > 0.) ? dt*mEvapTot[]/rho[]*area/Delta : 0.;
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
  foreach_interfacial_plic (f, F_ERR) {
    double delta_alpha = 0.;
    double magn = sqrt (sq (m.x) + sq (m.y) + sq (m.z));
    delta_alpha = (rho[] > 0.) ? dt*mEvapTot[]*magn/rho[]/Delta : 0.;
    double ff = plane_volume (m, alpha + delta_alpha);
    f[] = (ff > F_ERR && ff < 1.-F_ERR) ? ff : (ff <= F_ERR) ? 0. : 1.;
  }
}

/**
## References

~~~bib
@hal{cipriano2023, hal-04252476}

@article{malan2021,
  title = {A geometric VOF method for interface resolved phase change and
           conservative thermal energy advection},
  author = {Malan, L. C. and Malan, A. G. and Zaleski, S. and Rousseau, P. G.},
  journal = {Journal of Computational Physics},
  volume = {426},
  pages = {109920},
  year = {2021},
  doi = {10.1016/j.jcp.2020.109920}
}
~~~
*/
