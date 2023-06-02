/**
# Incompressible Navier--Stokes solver (MAC formulation)

We wish to approximate numerically the incompressible Navier--Stokes
equations
$$
\partial_t\mathbf{u}+\nabla\cdot(\mathbf{u}\otimes\mathbf{u}) = 
-\nabla p + \nabla\cdot(\nu\nabla\mathbf{u})
$$
$$
\nabla\cdot\mathbf{u} = 0
$$

We will use the generic time loop, a CFL-limited timestep and we will
need to solve a Poisson problem. */

#include "run.h"
#include "timestep.h"
#include "poisson.h"

/**
The Markers-And-Cells (MAC) formulation was first described in the
pioneering paper of [Harlow and Welch,
1965](/src/references.bib#harlow1965). It relies on a *face*
discretisation of the velocity components `u.x` and `u.y`, relative to
the (centered) pressure `p`. This guarantees the consistency of the
discrete gradient, divergence and Laplacian operators and leads to a
stable (mode-free) integration. */

scalar p[];
face vector uf[];

/**
In the case of variable density, the user will need to define both the
face and centered specific volume fields ($\alpha$ and $\alpha_c$
respectively) i.e. $1/\rho$. If not specified by the user, these
fields are set to one i.e. the density is unity.

Viscosity is set by defining the face dynamic viscosity $\mu$; default
is zero.

The face field $\mathbf{a}$ defines the acceleration term; default is
zero.

The statistics for the (multigrid) solution of the pressure Poisson
problems and implicit viscosity are stored in *mgp*, *mgpf*, *mgu*
respectively. 

If *stokes* is set to *true*, the velocity advection term
$\nabla\cdot(\mathbf{u}\otimes\mathbf{u})$ is omitted. This is a
reference to [Stokes flows](http://en.wikipedia.org/wiki/Stokes_flow)
for which inertia is negligible compared to viscosity. */

(const) face vector mu = zerof, a = zerof, alpha = unityf;
(const) scalar rho = unity;

/**
The volume expansion term is declared in
[evaporation.h](/sandbox/ecipriano/src/evaporation.h). */

extern scalar stefanflow;
#ifdef VARPROP
extern scalar drhodt;
#endif

/**
## Helper functions

We define the function that performs the projection
step with the volume expansion term due to the phase
change. */

trace
mgstats project_sf (struct Project q)
{
  face vector uf = q.uf;
  scalar p = q.p;
  (const) face vector alpha = q.alpha.x.i ? q.alpha : unityf;
  double dt = q.dt ? q.dt : 1.;
  int nrelax = q.nrelax ? q.nrelax : 4;
  
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

  foreach() {
    div[] += stefanflow[]/dt;
#ifdef VARPROP
    div[] += drhodt[]/dt;
#endif
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

/**
The only parameter is the viscosity coefficient $\nu.y[-1,0]*nu$.

The statistics for the (multigrid) solution of the Poisson problem are
stored in `mgp`. */

face vector nu[];
mgstats mgp;

#if EMBED
# define neumann_pressure(i) (alpha.n[i] ? a.n[i]*fm.n[i]/alpha.n[i] :	\
			      a.n[i]*rho[]/(cm[] + SEPS))
#else
# define neumann_pressure(i) (a.n[i]*fm.n[i]/alpha.n[i])
#endif

p[right] = neumann (neumann_pressure(ghost));
p[left]  = neumann (- neumann_pressure(0));

#if AXI
uf.n[bottom] = 0.;
uf.t[bottom] = dirichlet(0); // since uf is multiplied by the metric which
                             // is zero on the axis of symmetry
p[top]    = neumann (neumann_pressure(ghost));
#else // !AXI
#  if dimension > 1
p[top]    = neumann (neumann_pressure(ghost));
p[bottom] = neumann (- neumann_pressure(0));
#  endif
#  if dimension > 2
p[front]  = neumann (neumann_pressure(ghost));
p[back]   = neumann (- neumann_pressure(0));
#  endif
#endif // !AXI

/**
## Time integration

### Advection--Diffusion 

In a first step, we compute $\mathbf{u}_*$
$$
\frac{\mathbf{u}_* - \mathbf{u}_n}{dt} = \nabla\cdot\mathbf{S}
$$
with $\mathbf{S}$ the symmetric tensor
$$
\mathbf{S} = - \mathbf{u}\otimes\mathbf{u} + \nu\nabla\mathbf{u} =
\left(\begin{array}{cc}
- u_x^2 + 2\nu\partial_xu_x & - u_xu_y + \nu(\partial_yu_x + \partial_xu_y)\\
\ldots & - u_y^2 + 2\nu\partial_yu_y
\end{array}\right)
$$ 

The timestep for this iteration is controlled by the CFL condition
(and the timing of upcoming events). */

double dtmax;

event set_dtmax (i++,last) dtmax = DT;

event stability (i++,last) {
  dt = dtnext (timestep (uf, dtmax));
}

/**
If we are using VOF or diffuse tracers, we need to advance them (to
time $t+\Delta t/2$) here. Note that this assumes that tracer fields
are defined at time $t-\Delta t/2$ i.e. are lagging the
velocity/pressure fields by half a timestep. */

event vof (i++,last);
event tracer_advection (i++,last);
event tracer_diffusion (i++,last);

/**
The fluid properties such as specific volume (fields $\alpha$ and
$\alpha_c$) or dynamic viscosity (face field $\mu_f$) -- at time
$t+\Delta t/2$ -- can be defined by overloading this event. */

event properties (i++,last);

event advance (i++,last)
{
  /**
  We update the dynamic viscosity. */

  foreach_face()
    nu.x[] = mu.x[]*alpha.x[];

  /**
  We reset the acceleration field (if it is not a constant). */

  if (!is_constant(a.x)) {
    face vector af = a;
    trash ({af});
    foreach_face()
      af.x[] = 0.;
  }

  /**
  We allocate a local symmetric tensor field. To be able to compute the
  divergence of the tensor at the face locations, we need to
  compute the diagonal components at the center of cells and the
  off-diagonal component at the vertices. 
  
  ![Staggering of $\mathbf{u}$ and $\mathbf{S}$](/src/figures/Sxx.svg) */
  
  symmetric tensor S[]; // fixme: this does not work on trees

  /**
  We average the velocity components at the center to compute the
  diagonal components. */

  foreach()
    foreach_dimension()
      S.x.x[] = - sq(uf.x[] + uf.x[1,0])/4. + 2.*(nu.x[1,0]*uf.x[1,0] - nu.x[]*uf.x[])/Delta;

  /**
  We average horizontally and vertically to compute the off-diagonal
  component at the vertices. */

  foreach_vertex()
    S.x.y[] = 
      - (uf.x[] + uf.x[0,-1])*(uf.y[] + uf.y[-1,0])/4. +
      (nu.x[]*uf.x[] - nu.x[0,-1]*uf.x[0,-1] + nu.y[]*uf.y[] - nu.y[-1,0]*uf.y[-1,0])/Delta;

  /**
  Finally we compute
  $$
  \mathbf{u}_* = \mathbf{u}_n + dt\nabla\cdot\mathbf{S}
  $$ */

  foreach_face()
    uf.x[] += dt*(S.x.x[] - S.x.x[-1,0] + S.x.y[0,1] - S.x.y[])/Delta;


  //vertex scalar uv[];
  //foreach_vertex()
  //  uv[] = 0.25*(uf.x[] + uf.x[0,-1])*(uf.y[] + uf.y[-1,0]);

  //vector uc[];
  //foreach()
  //  foreach_dimension()
  //    uc.x[] = 0.25*sq(uf.x[1] + uf.x[]);

  //foreach_face(x) {
  //  double ue = uc.x[];
  //  double uw = uc.x[-1];
  //  double un = uv[1];
  //  double us = uv[];

  //  double Aij = (ue - uw + un - us)/Delta;

  //  double de = 0.5*(nu.x[1,0] + nu.x[])/Delta*((uf.x[1,0] - uf.x[])/Delta);
  //  double dw = 0.5*(nu.x[] + nu.x[-1,0])/Delta*((uf.x[] - uf.x[-1,0])/Delta);
  //  double dn = 0.5*(nu.x[0,1] + nu.x[])/Delta*((uf.x[0,1] - uf.x[])/Delta);
  //  double ds = 0.5*(nu.x[] + nu.x[0,-1])/Delta*((uf.x[] - uf.x[0,-1])/Delta);

  //  double Dij = (de - dw + dn - ds);

  //  uf.x[] += dt*(-Aij+Dij);
  //}

  //foreach_face(y) {

  //  double un = uc.y[];
  //  double us = uc.y[-1];
  //  double ue = uv[1];
  //  double uw = uv[];

  //  double Aij = (ue - uw + un - us)/Delta;

  //  double de = 0.5*(nu.y[1,0] + nu.y[])/Delta*((uf.y[1,0] - uf.y[])/Delta);
  //  double dw = 0.5*(nu.y[] + nu.y[-1,0])/Delta*((uf.y[] - uf.y[-1,0])/Delta);
  //  double dn = 0.5*(nu.y[0,1] + nu.y[])/Delta*((uf.y[0,1] - uf.y[])/Delta);
  //  double ds = 0.5*(nu.y[] + nu.y[0,-1])/Delta*((uf.y[] - uf.y[0,-1])/Delta);

  //  double Dij = (de - dw + dn - ds);

  //  uf.y[] += dt*(-Aij+Dij);
  //}
}

/**
### Projection 

In a second step we compute
$$
\mathbf{u}_{n+1} = \mathbf{u}_* - \Delta t\nabla p
$$
with the condition
$$
\nabla\cdot\mathbf{u}_{n+1} = 0
$$
This gives the Poisson equation for the pressure
$$
\nabla\cdot(\nabla p) = \frac{\nabla\cdot\mathbf{u}_*}{\Delta t}
$$ */

event projection (i++,last)
{
  mgp = project_sf (uf, p, alpha, dt, mgp.nrelax);
}

/**
### Acceleration term

We add the acceleration terms. */

event acceleration (i++,last)
{
  foreach_face()
    uf.x[] += dt*alpha.x[]*a.x[];
}

/**
Some derived solvers need to hook themselves at the end of the
timestep. */

event end_timestep (i++, last);

