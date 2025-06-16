/**
# Incompressible Navier--Stokes solver (centered formulation)

We wish to approximate numerically the incompressible,
variable-density Navier--Stokes equations
$$
\partial_t\mathbf{u}+\nabla\cdot(\mathbf{u}\otimes\mathbf{u}) = 
\frac{1}{\rho}\left[-\nabla p + \nabla\cdot(2\mu\mathbf{D})\right] + 
\mathbf{a}
$$
$$
\nabla\cdot\mathbf{u} = 0
$$
with the deformation tensor 
$\mathbf{D}=[\nabla\mathbf{u} + (\nabla\mathbf{u})^T]/2$.

The scheme implemented here is close to that used in Gerris ([Popinet,
2003](/src/references.bib#popinet2003), [Popinet,
2009](/src/references.bib#popinet2009), [Lagrée et al,
2011](/src/references.bib#lagree2011)).

We will use the generic time loop, a CFL-limited timestep, the
Bell-Collela-Glaz advection scheme and the implicit viscosity
solver. If embedded boundaries are used, a different scheme is used
for viscosity. */

#include "run.h"
#include "timestep.h"
#include "bcg.h"
#if EMBED
# include "viscosity-embed.h"
#else
# include "viscosity.h"
#endif

/**
The primary variables are the centered pressure field $p$ and the
centered velocity field $\mathbf{u}$. The centered vector field
$\mathbf{g}$ will contain pressure gradients and acceleration terms.

We will also need an auxilliary face velocity field $\mathbf{u}_f$ and
the associated centered pressure field $p_f$. */

int nv = 1;

scalar p[], * plist = NULL;
vector u[], * ulist = NULL;
vector g[], * glist = NULL;
scalar pf[], * pflist = NULL;
face vector uf[], * uflist = NULL;

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
mgstats mgp = {0}, mgpf = {0}, mgu = {0};
mgstats * mgplist = NULL, * mgpflist = NULL, * mgulist = NULL;
bool stokes = false;

/**
## Boundary conditions

For the default symmetric boundary conditions, we need to ensure that
the normal component of the velocity is zero after projection. This
means that, at the boundary, the acceleration $\mathbf{a}$ must be
balanced by the pressure gradient. Taking care of boundary orientation
and staggering of $\mathbf{a}$, this can be written */

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
For [embedded boundaries on trees](/src/embed-tree.h), we need to
define the pressure gradient for prolongation of pressure close to
embedded boundaries. */

#if TREE && EMBED
void pressure_embed_gradient (Point point, scalar p, coord * g)
{
  foreach_dimension()
    g->x = rho[]/(cm[] + SEPS)*(a.x[] + a.x[1])/2.;
}
#endif // TREE && EMBED

/**
## Initial conditions */

event defaults (i = 0)
{
  /**
  We fill the lists of velocities and pressures with the
  default fields. */

  plist = list_add (plist, p);
  ulist = vectors_add (ulist, u);
  glist = vectors_add (glist, g);
  pflist = list_add (pflist, pf);
  uflist = vectors_add (uflist, uf);

  /**
  If multiple velocity fields are considered, they are created
  and added to the respective lists. */

  for (int i = 1; i < nv; i++) {
    struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};

    // Pressures
    {
      scalar p = new scalar;
      char name[80];
      sprintf (name, "p%d", i);
      free (p.name);
      p.name = strdup (name);
      plist = list_add (plist, p);
    }

    // Velocities
    {
      vector u = new vector;
      foreach_dimension() {
        char name[80];
        sprintf (name, "u%d%s", i, ext.x);
        free (u.x.name);
        u.x.name = strdup (name);
      }
      ulist = vectors_add (ulist, u);
    }

    // Acceleration + pressure grads
    {
      vector g = new vector;
      foreach_dimension() {
        char name[80];
        sprintf (name, "g%d%s", i, ext.x);
        free (g.x.name);
        g.x.name = strdup (name);
      }
      glist = vectors_add (glist, g);
    }

    // Fractional pressure
    {
      scalar pf = new scalar;
      char name[80];
      sprintf (name, "pf%d", i);
      free (pf.name);
      pf.name = strdup (name);
      pflist = list_add (pflist, pf);
    }

    // Face velocities
    {
      face vector uf = new face vector;
      foreach_dimension() {
        char name[80];
        sprintf (name, "uf%d%s", i, ext.x);
        free (uf.x.name);
        uf.x.name = strdup (name);
      }
      uflist = vectors_add (uflist, uf);
    }
  }

  // pf inherits BCs from p
  for (int d = 0; d < nboundary; d++) {
    pf.boundary[d] = p.boundary[d];
    pf.boundary_homogeneous[d] = p.boundary_homogeneous[d];
  }

  // Copy the boundary conditions for pressure
  for (int i = 1; i < nv; i++) {
    scalar p = plist[i], p0 = plist[0];
    scalar pf = pflist[i], pf0 = pflist[0];
    vector u = ulist[i], u0 = ulist[0];
    vector g = glist[i], g0 = glist[0];
    face vector uf = uflist[i], uf0 = uflist[0];
#if TREE
    p.refine = p0.refine;
    pf.refine = pf0.refine;
    p.prolongation = p0.prolongation ;
    pf.prolongation = pf0.prolongation ;

    foreach_dimension() {
      u.x.refine = u0.x.refine;
      uf.x.refine = uf0.x.refine;
      u.x.prolongation = u0.x.prolongation ;
      uf.x.prolongation = uf0.x.prolongation ;
    }
#endif
    for (int d = 0; d < nboundary; d++) {
      p.boundary[d] = p0.boundary[d];
      p.boundary_homogeneous[d] = p0.boundary_homogeneous[d];
      pf.boundary[d] = pf0.boundary[d];
      pf.boundary_homogeneous[d] = pf0.boundary_homogeneous[d];

      foreach_dimension() {
        u.x.boundary[d] = u0.x.boundary[d];
        u.x.boundary_homogeneous[d] = u0.x.boundary_homogeneous[d];
        g.x.boundary[d] = g0.x.boundary[d];
        g.x.boundary_homogeneous[d] = g0.x.boundary_homogeneous[d];
        uf.x.boundary[d] = uf0.x.boundary[d];
        uf.x.boundary_homogeneous[d] = uf0.x.boundary_homogeneous[d];
      }
    }
  }

  // Reset fields
  for (scalar p in plist)
    foreach()
      p[] = 0.;

  for (scalar pf in pflist)
    foreach()
      pf[] = 0.;

  for (vector u in ulist)
    foreach()
      foreach_dimension()
        u.x[] = 0.;

  for (vector g in glist)
    foreach()
      foreach_dimension()
        g.x[] = 0.;

  for (face vector uf in uflist)
    foreach_face()
      foreach_dimension()
        uf.x[] = 0.;

#if 0 // Debug
  list_print (plist, stdout);
  list_print ((scalar *)ulist, stdout);
  list_print ((scalar *)glist, stdout);
  list_print (pflist, stdout);
  list_print ((scalar *)uflist, stdout);
#endif

  /**
  We create the lists with convergence statistics. */

  mgplist = (mgstats *)malloc (nv*sizeof (mgstats));
  mgpflist = (mgstats *)malloc (nv*sizeof (mgstats));
  mgulist = (mgstats *)malloc (nv*sizeof (mgstats));

  /**
  We reset the multigrid parameters to their default values. */
  
  mgp = (mgstats){0};
  mgpf = (mgstats){0};
  mgu = (mgstats){0};  

  mgplist[0] = mgp;
  mgpflist[0] = mgpf;
  mgulist[0] = mgu;

  for (int i = 1; i < nv; i++) {
    mgstats mgp = mgplist[i];
    mgstats mgpf = mgpflist[i];
    mgstats mgu = mgulist[i];

    mgp = (mgstats){0};
    mgpf = (mgstats){0};
    mgu = (mgstats){0};

    NOT_UNUSED (mgp);
    NOT_UNUSED (mgpf);
    NOT_UNUSED (mgu);
  }

  
  CFL = 0.8;

  /**
  The pressures are never dumped. */

  scalar p, pf;
  for (p,pf in plist,pflist)
    p.nodump = pf.nodump = true;
  
  /**
  The default density field is set to unity (times the metric and the
  solid factors). */

  if (alpha.x.i == unityf.x.i) {
    alpha = fm;
    rho = cm;
  }
  else if (!is_constant(alpha.x)) {
    face vector alphav = alpha;
    foreach_face()
      alphav.x[] = fm.x[];
  }

  /**
  On trees, refinement of the face-centered velocity field needs to
  preserve the divergence-free condition. */

#if TREE
  for (face vector uf in uflist)
    uf.x.refine = refine_face_solenoidal;

  /**
  When using [embedded boundaries](/src/embed.h), the restriction and
  prolongation operators need to take the boundary into account. */

#if EMBED // TODO
  uf.x.refine = refine_face;
  foreach_dimension()
    uf.x.prolongation = refine_embed_face_x;
  for (scalar s in {p, pf, u, g}) {
    s.restriction = restriction_embed_linear;
    s.refine = s.prolongation = refine_embed_linear;
    s.depends = list_add (s.depends, cs);
  }
  for (scalar s in {p, pf})
    s.embed_gradient = pressure_embed_gradient;
#endif // EMBED
#endif // TREE

  /**
  We set the dimensions of the velocity field. */

  foreach()
    for (scalar u in ulist)
      foreach_dimension()
        dimensional (u.x[] == Delta/t);
}


/**
We had some objects to display by default. */

event default_display (i = 0)
  display ("squares (color = 'u.x', spread = -1);");

/**
After user initialisation, we initialise the face velocity and fluid
properties. */

double dtmax;

event init (i = 0)
{
  vector u;
  face vector uf;
  for (u,uf in ulist,uflist) {
    trash ({uf});
    foreach_face()
      uf.x[] = fm.x[]*face_value (u.x, 0);
  }

  /**
  We update fluid properties. */

  event ("properties");

  /**
  We set the initial timestep (this is useful only when restoring from
  a previous run). */

  dtmax = DT;
  event ("stability");
}

event cleanup (t = end) {
  free (mgplist), mgplist = NULL;
  free (mgpflist), mgpflist = NULL;
  free (mgulist), mgulist = NULL;

  for (int i = 1; i < nv; i++) {
    scalar p = plist[i];
    vector u = ulist[i];
    vector g = glist[i];
    scalar pf = pflist[i];
    face vector uf = uflist[i];

    delete ({p});
    delete ((scalar *){u});
    delete ((scalar *){g});
    delete ({pf});
    delete ((scalar *){uf});
  }
  free (plist), plist = NULL;
  free (ulist), ulist = NULL;
  free (glist), glist = NULL;
  free (pflist), pflist = NULL;
  free (uflist), uflist = NULL;
}

/**
## Time integration

The timestep for this iteration is controlled by the CFL condition,
applied to the face centered velocity field $\mathbf{u}_f$; and the
timing of upcoming events. */

event set_dtmax (i++,last) dtmax = DT;

event stability (i++,last) {
  for (face vector uf in uflist)
    dt = dtnext (stokes ? dtmax : timestep (uf, dtmax));
}

/**
If we are using VOF or diffuse tracers, we need to advance them (to
time $t+\Delta t/2$) here. Note that this assumes that tracer fields
are defined at time $t-\Delta t/2$ i.e. are lagging the
velocity/pressure fields by half a timestep. */

event vof (i++,last);
event vof_sources (i++,last);
event tracer_advection (i++,last);
event tracer_diffusion (i++,last);

/**
The fluid properties such as specific volume (fields $\alpha$ and
$\alpha_c$) or dynamic viscosity (face field $\mu_f$) -- at time
$t+\Delta t/2$ -- can be defined by overloading this event. */

event properties (i++,last);

/**
### Predicted face velocity field

For second-order in time integration of the velocity advection term
$\nabla\cdot(\mathbf{u}\otimes\mathbf{u})$, we need to define the face
velocity field $\mathbf{u}_f$ at time $t+\Delta t/2$. We use a version
of the Bell-Collela-Glaz [advection scheme](/src/bcg.h) and the
pressure gradient and acceleration terms at time $t$ (stored in vector
$\mathbf{g}$). */

void prediction (vector u, vector g, face vector uf)
{
  vector du;
  foreach_dimension() {
    scalar s = new scalar;
    du.x = s;
  }

  if (u.x.gradient)
    foreach()
      foreach_dimension() {
#if EMBED
        if (!fs.x[] || !fs.x[1])
	  du.x[] = 0.;
	else
#endif
	  du.x[] = u.x.gradient (u.x[-1], u.x[], u.x[1])/Delta;
      }
  else
    foreach()
      foreach_dimension() {
#if EMBED
        if (!fs.x[] || !fs.x[1])
	  du.x[] = 0.;
	else
#endif
	  du.x[] = (u.x[1] - u.x[-1])/(2.*Delta);
    }

  trash ({uf});
  foreach_face() {
    double un = dt*(u.x[] + u.x[-1])/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    uf.x[] = u.x[i] + (g.x[] + g.x[-1])*dt/4. + s*(1. - s*un)*du.x[i]*Delta/2.;
    #if dimension > 1
    if (fm.y[i,0] && fm.y[i,1]) {
      double fyy = u.y[i] < 0. ? u.x[i,1] - u.x[i] : u.x[i] - u.x[i,-1];
      uf.x[] -= dt*u.y[i]*fyy/(2.*Delta);
    }
    #endif
    #if dimension > 2
    if (fm.z[i,0,0] && fm.z[i,0,1]) {
      double fzz = u.z[i] < 0. ? u.x[i,0,1] - u.x[i] : u.x[i] - u.x[i,0,-1];
      uf.x[] -= dt*u.z[i]*fzz/(2.*Delta);
    }
    #endif
    uf.x[] *= fm.x[];
  }

  delete ((scalar *){du});
}

/**
### Advection term

We predict the face velocity field $\mathbf{u}_f$ at time $t+\Delta
t/2$ then project it to make it divergence-free. We can then use it to
compute the velocity advection term, using the standard
Bell-Collela-Glaz advection scheme for each component of the velocity
field. */

event advection_term (i++,last)
{
  if (!stokes) {
    vector u, g;
    face vector uf;

    for (u,uf,g in ulist,uflist,glist)
      prediction (u, g, uf);
#if VELOCITY_JUMP
    mgpf = project (uflist, pflist, alpha, dt/2., mgpf.nrelax);
#else
    for (int i = 0; i < nv; i++) {
      scalar pf = pflist[i];
      uf = uflist[i];
      mgpflist[i] = project (uf, pf, alpha, dt/2., mgpflist[i].nrelax);
    }
#endif
    for (u,uf,g in ulist,uflist,glist)
      advection ((scalar *){u}, uf, dt, (scalar *){g});
  }
}

/**
### Viscous term

We first define a function which adds the pressure gradient and
acceleration terms. */

static void correction (double dt)
{
  vector u, g;
  for (u,g in ulist,glist)
    foreach()
      foreach_dimension()
        u.x[] += dt*g.x[];
}

/**
The viscous term is computed implicitly. We first add the pressure
gradient and acceleration terms, as computed at time $t$, then call
the implicit viscosity solver. We then remove the acceleration and
pressure gradient terms as they will be replaced by their values at
time $t+\Delta t$. */

event viscous_term (i++,last)
{
  if (constant(mu.x) != 0.) {
    correction (dt);
    for (int i = 0; i < nv; i++) {
      vector u = ulist[i];
      mgu = mgulist[i];
      mgu = viscosity (u, mu, rho, dt, mgu.nrelax);
    }
    correction (-dt);
  }

  /**
  We reset the acceleration field (if it is not a constant). */

  if (!is_constant(a.x)) {
    face vector af = a;
    trash ({af});
    foreach_face()
      af.x[] = 0.;
  }
}

/**
### Acceleration term

The acceleration term $\mathbf{a}$ needs careful treatment as many
equilibrium solutions depend on exact balance between the acceleration
term and the pressure gradient: for example Laplace's balance for
surface tension or hydrostatic pressure in the presence of gravity.

To ensure a consistent discretisation, the acceleration term is
defined on faces as are pressure gradients and the centered combined
acceleration and pressure gradient term $\mathbf{g}$ is obtained by
averaging. 

The (provisionary) face velocity field at time $t+\Delta t$ is
obtained by interpolation from the centered velocity field. The
acceleration term is added. */

event acceleration (i++,last)
{
  vector u, uf;
  for (u,uf in ulist, uflist) {
  trash ({uf});
  foreach_face()
    uf.x[] = fm.x[]*(face_value (u.x, 0) + dt*a.x[]);
  }
}

/**
## Approximate projection

This function constructs the centered pressure gradient and
acceleration field *g* using the face-centered acceleration field *a*
and the cell-centered pressure field *p*. */

void centered_gradient (scalar p, vector g)
{

  /**
  We first compute a face field $\mathbf{g}_f$ combining both
  acceleration and pressure gradient. */

  face vector gf[];
  foreach_face()
    gf.x[] = fm.x[]*a.x[] - alpha.x[]*(p[] - p[-1])/Delta;

  /**
  We average these face values to obtain the centered, combined
  acceleration and pressure gradient field. */

  trash ({g});
  foreach()
    foreach_dimension()
      g.x[] = (gf.x[] + gf.x[1])/(fm.x[] + fm.x[1] + SEPS);
}

/**
To get the pressure field at time $t + \Delta t$ we project the face
velocity field (which will also be used for tracer advection at the
next timestep). Then compute the centered gradient field *g*. */

event projection (i++,last)
{
  scalar p;
  vector g;

#if VELOCITY_JUMP
  mgp = project (uflist, plist, alpha, dt, mgp.nrelax);
#else
  for (int i = 0; i < nv; i++) {
    face vector uf = uflist[i];
    p = plist[i];
    mgplist[i] = project (uf, p, alpha, dt, mgplist[i].nrelax);
  }
#endif
  for (p,g in plist,glist)
    centered_gradient (p, g);

  /**
  We add the gradient field *g* to the centered velocity field. */

  correction (dt);
}

/**
Some derived solvers need to hook themselves at the end of the
timestep. */

event end_timestep (i++, last);

/**
## Adaptivity

After mesh adaptation fluid properties need to be updated. When using
[embedded boundaries](/src/embed.h) the fluid fractions and face
fluxes need to be checked for inconsistencies. */

#if TREE
event adapt (i++,last) {
#if EMBED
  fractions_cleanup (cs, fs);
  foreach_face()
    for (face vector uf in uflist)
      if (uf.x[] && !fs.x[])
        uf.x[] = 0.;
#endif
  event ("properties");
}
#endif

/**
## See also

* [Double projection](double-projection.h)
* [Performance monitoring](perfs.h)
*/
