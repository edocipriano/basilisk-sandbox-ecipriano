/**
# Incompressible Navier--Stokes solver with divergence source term (centered formulation)

This extension of the [centered.h](/src/navier-stokes/centered.h)
Navier--Stokes equations solver considers situations where the
velocity divergence is not null. Examples are phase change simulations,
or Low Mach systems with variable physical properties.

The modified set of equations reads:

$$
\partial_t\mathbf{u}+\nabla\cdot(\mathbf{u}\otimes\mathbf{u}) = 
\frac{1}{\rho}\left[-\nabla p + \nabla\cdot(2\mu\mathbf{D})\right] + 
\mathbf{a}
$$
$$
\nabla\cdot\mathbf{u} =
  {\color{blue} \dot{m} \left(\dfrac{1}{\rho_g}
  - \dfrac{1}{\rho_l}\right)\delta_\Gamma
  -\color{blue} \dfrac{1}{\rho}\dfrac{D\rho}{Dt}}
$$
with the deformation tensor 
$\mathbf{D}=[\nabla\mathbf{u} + (\nabla\mathbf{u})^T]/2$.
*/

/**
## Field Allocations

We define scalar fields with two possible sources of divergence:
`intexp` is the phase change source term, localized at the
gas-liquid interface, while `drhodt` considers density changes. */

#define LOW_MACH 1

scalar drhodt[], * drhodtlist = NULL;
scalar intexp[], * intexplist = NULL;

bool no_advection_div = false, closed = false;

/**
## Closed systems

If the domain is closed, the total volume cannot change, and the net expansion
must be balanced by a variation of the thermodynamic pressure $P_0$, which is
uniform in space for a low Mach number system. Splitting the density variations
into a thermodynamic pressure contribution and into the remaining effects
(temperature, composition, and phase change), the divergence reads:
$$
  \nabla\cdot\mathbf{u} = S - \chi_T\dfrac{dP_0}{dt}
$$
where $\chi_T = 1/\rho\left(\partial\rho/\partial P\right)_{T,\omega}$ is the
isothermal compressibility. Integrating over the closed domain, where
$\int_\Omega\nabla\cdot\mathbf{u}\,dV = 0$, the pressurization rate is obtained:
$$
  \dfrac{dP_0}{dt} = \dfrac{\int_\Omega S\,dV}{\int_\Omega \chi_T\,dV}
$$
The compensation is weighted on the *local* compressibility: since the gas phase
is about five orders of magnitude more compressible than the liquid phase, the
ullage absorbs almost the entire volume variation, while the liquid remains
still. The field `chiT` must be filled by the module that computes the material
properties; if it is left to zero everywhere, the net expansion is redistributed
uniformly over the domain, which is the behaviour of an incompressible closed
system.

The thermodynamic pressure `P0` is integrated in time at the end of every time
step, and it should be initialized by the user (or by the phase change model) to
the initial pressure of the system. The pressurization rate is computed by
`project_lowmach()`, therefore it is not available when the velocity jump
formulation is used. */

scalar chiT[];
double P0 = 0., dP0dt = 0.;

/**
## Projection Function

We define the function that performs the projection step with the volume
expansion term due to the phase change or due to density changes. */

#include "poisson.h"

trace
mgstats project_lowmach (face vector uf, scalar p,
     (const) face vector alpha = unityf,
     double dt = 1.,
     int nrelax = 4)
{
  
  /**
  We allocate a local scalar field and compute the divergence of
  $\mathbf{u}_f$. The divergence is scaled by *dt* so that the
  pressure has the correct dimension. */

  scalar div[];
  foreach() {
    div[] = 0.;
    foreach_dimension()
      div[] += uf.x[1] - uf.x[];
    div[] /= dt*Delta;
  }

  /**
  We add the volume expansion contribution. */

  extern int nv;
  static int inv = 0;
  scalar intexp = intexplist[inv];
  scalar drhodt = drhodtlist[inv];

  double volume = 0., srcsum = 0., chisum = 0.;
  if (closed) {
    foreach (reduction(+:volume) reduction(+:srcsum) reduction(+:chisum)) {
      volume += dv();
      srcsum += (intexp[] + drhodt[])*dv();
      chisum += chiT[]*dv();
    }
    /**
    The pressurization rate is computed just once, using the fields of the
    whole domain, because the thermodynamic pressure is a property of the
    system and not of the single velocity field. */

    if (inv == 0)
      dP0dt = (chisum > 0.) ? -srcsum/chisum : 0.;
  }

  /**
  If the field `chiT` is zero (by default), the pressurization is evenly
  distributed in the domain. */

  foreach() {
    div[] += (intexp[] + drhodt[])/dt;
    if (closed)
      div[] += ((chisum > 0.) ? chiT[]*dP0dt : -srcsum/volume)/dt;
  }
  inv++;
  inv = (inv == nv) ? 0 : inv;

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

/**
## Advection Function

We overwrite the advection function for non-incompressible flows, by removing
the divergence of the velocity from the convective term of the momentum
equation, in order to be consistent with the non-conservative formulation.

TODO: embed compatibility */

#include "utils.h"
#include "bcg.h"

void advection_div (scalar * tracers, face vector u, double dt,
		scalar * src = NULL)
{
  /**
  If *src* is not provided we set all the source terms to zero. */
  
  scalar * psrc = src;
  if (!src)
    for (scalar s in tracers) {
      const scalar zero[] = 0.;
      src = list_append (src, zero);
    }
  assert (list_len (tracers) == list_len (src));

  scalar f, source;
  for (f,source in tracers,src) {
    face vector flux[];
    tracer_fluxes (f, u, flux, dt, source);
#if !EMBED
    foreach() {
      double fold = f[];
      NOT_UNUSED (fold);
      foreach_dimension() {
        if (no_advection_div)
          f[] += dt*(flux.x[] - flux.x[1] + fold*(u.x[1] - u.x[]))/(Delta*cm[]);
        else
          f[] += dt*(flux.x[] - flux.x[1])/(Delta*cm[]);
      }
    }
#else // EMBED
    update_tracer (f, u, flux, dt);
#endif // EMBED
  }

  if (!psrc)
    free (src);
}

/**
We overwrite the functions `project()` and `advection()` in
[centered.h](/src/navier-stokes/centered.h) in order to call `project_lowmach()`
and `advection_div()` instead, accounting for the divergcence source terms. */

#if VELOCITY_JUMP
#define advection(...) advection_div(__VA_ARGS__)
#include "navier-stokes/centered-new.h"
#undef advection
#else
#define project(...) project_lowmach(__VA_ARGS__)
#define advection(...) advection_div(__VA_ARGS__)
#include "navier-stokes/centered-new.h"
#undef advection
#undef project
#endif

/**
## Thermodynamic Pressure

The pressurization rate is computed by `project_lowmach()` during the
`projection` event, therefore the thermodynamic pressure is integrated in time
afterwards, at the end of the time step. It stays constant unless the system is
`closed`. */

event end_timestep (i++) {
  if (closed)
    P0 += dP0dt*dt;
}

/**
We set the default divergence source term to zero (for the liquid phase) */

event defaults (i = 0) {

  /**
  The pressurization rate is used explicitly by the source terms of the next
  time step, therefore it must be reset at the beginning of every simulation.
  Otherwise, consecutive calls to `run()` (e.g. a convergence study) would start
  from the rate of the previous simulation. The thermodynamic pressure `P0` is
  not reset here, because it is initialized by the user or by the phase change
  model. */

  dP0dt = 0.;

  drhodtlist = list_add (drhodtlist, drhodt);
  intexplist = list_add (intexplist, intexp);

  for (int i = 1; i < nv; i++) {
    // Density changes
    {
      scalar drhodt = new scalar;
      char name[80];
      sprintf (name, "drhodt%d", i);
      free (drhodt.name);
      drhodt.name = strdup (name);
      drhodtlist = list_add (drhodtlist, drhodt);
    }

    // Stefan flow
    {
      scalar intexp = new scalar;
      char name[80];
      sprintf (name, "intexp%d", i);
      free (intexp.name);
      intexp.name = strdup (name);
      intexplist = list_add (intexplist, intexp);
    }
  }
}

event cleanup (t = end) {
  for (int i = 1; i < nv; i++) {
    scalar drhodt = drhodtlist[i];
    scalar intexp = intexplist[i];

    delete ({drhodt,intexp});
  }
  free (drhodtlist), drhodtlist = NULL;
  free (intexplist), intexplist = NULL;
}

