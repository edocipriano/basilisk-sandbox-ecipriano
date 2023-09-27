/**
# Velocity Potential Model

We want to obtain a divergence-free velocity extrapolation from the
field velocity with phase change. A possible approach, proposed in
different works with small differences ([Scapin et al. 2020](#scapin2020volume),
[Malan et al. 2021](malan2021geometric), [Palmore et al. 2019](#palmore2019volume)),
consists in solving an additional Poisson equation that allows a velocity
potential $\phi$ to be computed:
$$
\nabla \cdot \left( \alpha \nabla \phi \right)
=
\dfrac{\dot{m}}{\Delta t} \left( \frac{1}{\rho_g} - \frac{1}{\rho_l} \right) \delta_\Gamma \\
$$
The stefan velocity can be obtained from the velocity potential as:
$$
\mathbf{u}^S = -\Delta t \alpha \nabla \phi
$$
We then calculate the extended velocity by subtracting the stefan velocity
from the field velocity. The resulting extended velocity field will be
divergence-free by construction:
$$
\mathbf{u}^E = \mathbf{u} - \mathbf{u}^S
$$
*/

#include "common-evaporation.h"

extern scalar f;

/*

We write a multigrid solver for solving the velocity potential poisson projection step on a subdomain. 
To improve the convergence of the solver we use a flux corrected weighted interpolation similar to ____.

For the poisson projection problem we apply a homogenous neumann boundary condition on the pressure in the gas and no flux condition on the velocity in the liquid.
We apply these immersed boundary conditions by modifying both the relax and residual functions as well the coefficients.

*/
static inline double bilinear_weighted(Point point, scalar s, scalar q){
  #if dimension == 1
    return (3.*coarse(s)/coarse(q) + coarse(s,child.x)/coarse(q,child.x))/(3/coarse(q) + 1/coarse(q,child.x));
  #elif dimension == 2
    return (9.*coarse(s)/coarse(q)
       + 3.*(coarse(s,child.x)/coarse(q, child.x) + coarse(s,0,child.y)/coarse(q, 0,child.y))
       + coarse(s,child.x,child.y)/coarse(q, child.x,child.y))/
       (9/coarse(q) + 3*(1/coarse(q, child.x) + 1/coarse(q, 0,child.y)) + 1/coarse(q, child.x,child.y));
  #else // dimension == 3
    return (27.*coarse(s)/coarse(q)
      + 9.*(coarse(s,child.x)/coarse(q,child.x) + coarse(s,0,child.y)/coarse(q,0,child.y) + coarse(s,0,0,child.z)/coarse(q,0,0,child.z))
      + 3.*(coarse(s,child.x,child.y)/coarse(q,child.x,child.y)
      + coarse(s,child.x,0,child.z)/coarse(q,child.x,0,child.z) + coarse(s,0,child.y,child.z)/coarse(q,0, child.y,child.z))
      + coarse(s,child.x,child.y,child.z)/coarse(q, child.x,child.y,child.z))/
      (27./coarse(q)
        + 9.*(1/coarse(q,child.x) + 1/coarse(q,0,child.y) + 1/coarse(q,0,0,child.z))
        + 3.*(1/coarse(q,child.x,child.y) + 1/coarse(q,child.x,0,child.z) + 1/coarse(q,0, child.y,child.z))
        + 1/coarse(q, child.x,child.y,child.z));
  #endif
}

void mg_cycle_weighted (scalar * a, scalar * res, scalar * da,
         void (* relax) (scalar * da, scalar * res, 
             int depth, void * data),
         void * data,
         int nrelax, int minlevel, int maxlevel, scalar rho){

  /**
  We first define the residual on all levels. */

  restriction (res);

  /**
  We then proceed from the coarsest grid (*minlevel*) down to the
  finest grid. */

  minlevel = min (minlevel, maxlevel);
  for (int l = minlevel; l <= maxlevel; l++) {

    /**
    On the coarsest grid, we take zero as initial guess. */

    if (l == minlevel)
      foreach_level_or_leaf (l)
  for (scalar s in da)
    foreach_blockf (s)
      s[] = 0.;

    /**
    On all other grids, we take as initial guess the approximate solution
    on the coarser grid bilinearly interpolated onto the current grid. */

    else
      foreach_level (l)
  for (scalar s in da)
    foreach_blockf (s)
      s[] = bilinear_weighted(point, s, rho);
    
    /**
    We then apply homogeneous boundary conditions and do several
    iterations of the relaxation function to refine the initial guess. */

    boundary_level (da, l);
    for (int i = 0; i < nrelax; i++) {
      relax (da, res, l, data);
      boundary_level (da, l);
    }
  }

  /**
  And finally we apply the resulting correction to *a*. */

  foreach() {
    scalar s, ds;
    for (s, ds in a, da)
      foreach_blockf (s)
  s[] += ds[];
  }
}

/**
## Multigrid solver

The multigrid solver itself uses successive calls to the multigrid
cycle to refine an initial guess until a specified tolerance is
reached. 

The maximum number of iterations is controlled by *NITERMAX* and the
tolerance by *TOLERANCE* with the default values below. */

/**
The user needs to provide a function which computes the residual field
(and returns its maximum) as well as the relaxation function. The
user-defined pointer *data* can be used to pass arguments to these
functions. The optional number of relaxations is *nrelax* (default is
one) and *res* is an optional list of fields used to store the
residuals. The minimum level of the hierarchy can be set (default is
zero i.e. the root cell). */

struct MGSolve_weighted {
  scalar * a, * b;
  double (* residual) (scalar * a, scalar * b, scalar * res,
           void * data);
  void (* relax) (scalar * da, scalar * res, int depth, 
      void * data);
  void * data;
  
  int nrelax;
  scalar * res;
  int minlevel;
  double tolerance;
  scalar rho;
};

mgstats mg_solve_weighted (struct MGSolve_weighted p)
{

  /**
  We allocate a new correction and residual field for each of the scalars
  in *a*. */

  scalar * da = list_clone (p.a), * res = p.res;
  if (!res)
    res = list_clone (p.b);

  /**
  The boundary conditions for the correction fields are the
  *homogeneous* equivalent of the boundary conditions applied to
  *a*. */

  for (int b = 0; b < nboundary; b++)
    for (scalar s in da)
      s.boundary[b] = s.boundary_homogeneous[b];
  
  /**
  We initialise the structure storing convergence statistics. */

  mgstats s = {0};
  double sum = 0.;
  foreach (reduction(+:sum))
    for (scalar s in p.b)
      sum += s[];
  s.sum = sum;
  s.nrelax = p.nrelax > 0 ? p.nrelax : 4;
  
  /**
  Here we compute the initial residual field and its maximum. */

  double resb;
  resb = s.resb = s.resa = p.residual (p.a, p.b, res, p.data);

  /**
  We then iterate until convergence or until *NITERMAX* is reached. Note
  also that we force the solver to apply at least one cycle, even if the
  initial residual is lower than *TOLERANCE*. */

  if (p.tolerance == 0.)
    p.tolerance = TOLERANCE;
  for (s.i = 0;
       s.i < NITERMAX && (s.i < NITERMIN || s.resa > p.tolerance);
       s.i++) {
    mg_cycle_weighted(p.a, res, da, p.relax, p.data,
        s.nrelax,
        p.minlevel,
        grid->maxdepth, p.rho);
    s.resa = p.residual (p.a, p.b, res, p.data);

    /**
    We tune the number of relaxations so that the residual is reduced
    by between 2 and 20 for each cycle. This is particularly useful
    for stiff systems which may require a larger number of relaxations
    on the finest grid. */

#if 1
    if (s.resa > p.tolerance) {
      if (resb/s.resa < 1.2 && s.nrelax < 100)
  s.nrelax++;
      else if (resb/s.resa > 10 && s.nrelax > 2)
  s.nrelax--;
    }
#else
    if (s.resa == resb) /* convergence has stopped!! */
      break;
    if (s.resa > resb/1.1 && p.minlevel < grid->maxdepth)
      p.minlevel++;
#endif

    resb = s.resa;
  }
  s.minlevel = p.minlevel;
  
  /**
  If we have not satisfied the tolerance, we warn the user. */

  if (s.resa > p.tolerance){
    scalar v = p.a[0];
    fprintf (ferr, 
       "WARNING: convergence for %s not reached after %d iterations\n"
       "  res: %g sum: %g nrelax: %d\n", v.name,
       s.i, s.resa, s.sum, s.nrelax), fflush (ferr);
  }
    
  /**
  We deallocate the residual and correction fields and free the lists. */

  if (!p.res)
    delete (res), free (res);
  delete (da), free (da);

  return s;
}

/**
## Application to the Poisson--Helmholtz equation

We now apply the generic multigrid solver to the Poisson--Helmholtz equation
$$
\nabla\cdot (\alpha\nabla a) + \lambda a = b
$$
We first setup the data structure required to pass the extra
parameters $\alpha$ and $\lambda$. We define $\alpha$ as a face
vector field because we need values at the face locations
corresponding to the face gradients of field $a$. 

*alpha* and *lambda* are declared as *(const)* to indicate that the
function works also when *alpha* and *lambda* are constant vector
(resp. scalar) fields. If *tolerance* is set, it supersedes the
default *TOLERANCE* of the multigrid solver, *nrelax* controls the
initial number of relaxations (default is one), *minlevel* controls
the minimum level of the hierarchy (default is one) and *res* is an
optional list of fields used to store the final residual (which can be
useful to monitor convergence). 


To modify the solver to solve the problem on a subdomain*/

struct Poisson_weighted {
  scalar a, b;
  (const) face vector alpha;
  (const) scalar lambda;
  double tolerance;
  int nrelax, minlevel;
  scalar * res;
#if EMBED
  double (* embed_flux) (Point, scalar, vector, double *);
#endif
  scalar rho;
  scalar subdomain;
};

static void relax2 (scalar * al, scalar * bl, int l, void * data)
{
  scalar a = al[0], b = bl[0];
  struct Poisson_weighted * p = (struct Poisson_weighted *) data;
  (const) face vector alpha = p->alpha;
  (const) scalar lambda = p->lambda;
  (const) scalar subdomain = p->subdomain;

  /**
  We use either Jacobi (under)relaxation or we directly reuse values
  as soon as they are updated. For Jacobi, we need to allocate space
  for the new field *c*. Jacobi is useful mostly as it gives results
  which are independent of the order in which the cells are
  traversed. This is not the case for the simple traversal, which
  means for example that results will depend on whether a tree or
  a multigrid is used (because cells will be traversed in a different
  order). The same comment applies to OpenMP or MPI parallelism. In
  practice however Jacobi convergence tends to be slower than simple
  reuse. */
  
#if JACOBI
  scalar c[];
#else
  scalar c = a;
#endif
  
  /**
  We use the face values of $\alpha$ to weight the gradients of the
  5-points Laplacian operator. We get the relaxation function. */

  foreach_level_or_leaf (l) {
    double n = - sq(Delta)*b[], d = - lambda[]*sq(Delta);
    foreach_dimension() {
      n += alpha.x[1]*a[1] + alpha.x[]*a[-1];
      d += alpha.x[1] + alpha.x[];
    }
#if EMBED
    if (p->embed_flux) {
      double c, e = p->embed_flux (point, a, alpha, &c);
      n -= c*sq(Delta);
      d += e*sq(Delta);
    }
    if (!d)
      c[] = b[] = 0.;
    else
#endif // EMBED
      c[] = (d > 0.0 && subdomain[] > 0.0) ? n/d : 0.0;
  }

  /**
  For weighted Jacobi we under-relax with a weight of 2/3. */
  
#if JACOBI
  foreach_level_or_leaf (l)
    a[] = (a[] + 2.*c[])/3.;
#endif
  
#if TRASH
  scalar a1[];
  foreach_level_or_leaf (l)
    a1[] = a[];
  trash ({a});
  foreach_level_or_leaf (l)
    a[] = a1[];
#endif
}


static double residual2 (scalar * al, scalar * bl, scalar * resl, void * data)
{
  scalar a = al[0], b = bl[0], res = resl[0];
  struct Poisson_weighted * p = (struct Poisson_weighted *) data;
  (const) face vector alpha = p->alpha;
  (const) scalar lambda = p->lambda;
  (const) scalar subdomain = p->subdomain;

  double maxres = 0.;
#if TREE
  /* conservative coarse/fine discretisation (2nd order) */
  face vector g[];
  foreach_face()
    g.x[] = alpha.x[]*face_gradient_x (a, 0);
  foreach (reduction(max:maxres), nowarning) {
    res[] = b[] - lambda[]*a[];
    foreach_dimension()
      res[] -= (g.x[1] - g.x[])/Delta;
    res[] *= subdomain[];
#if EMBED
    if (p->embed_flux) {
      double c, e = p->embed_flux (point, a, alpha, &c);
      res[] += c - e*a[];
    }
#endif // EMBED    
    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
  }
#else // !TREE
  /* "naive" discretisation (only 1st order on trees) */
  foreach (reduction(max:maxres), nowarning) {
    res[] = b[] - lambda[]*a[];
    foreach_dimension()
      res[] += (alpha.x[0]*face_gradient_x (a, 0) - alpha.x[1]*face_gradient_x (a, 1))/Delta;
    res[] *= subdomain[];  
#if EMBED
    if (p->embed_flux) {
      double c, e = p->embed_flux (point, a, alpha, &c);
      res[] += c - e*a[];
    }
#endif // EMBED
    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
  }
#endif // !TREE
  return maxres;
}

/**
## User interface

Finally we provide a generic user interface for a Poisson--Helmholtz
equation of the form
$$
\nabla\cdot (\alpha\nabla a) + \lambda a = b
$$ */

mgstats poisson_weighted (struct Poisson_weighted p)
{

  /**
  If $\alpha$ or $\lambda$ are not set, we replace them with constant
  unity vector (resp. zero scalar) fields. Note that the user is free to
  provide $\alpha$ and $\beta$ as constant fields. */

  if (!p.alpha.x.i)
    p.alpha = unityf;
  if (!p.lambda.i)
    p.lambda = zeroc;

  /**
  We need $\alpha$ and $\lambda$ on all levels of the grid. */

  face vector alpha = p.alpha;
  scalar lambda = p.lambda;
  scalar rho = p.rho;
  scalar subdomain = p.subdomain;
  restriction ({alpha,lambda, rho, subdomain});

  /**
  If *tolerance* is set it supersedes the default of the multigrid
  solver. */

  double defaultol = TOLERANCE;
  if (p.tolerance)
    TOLERANCE = p.tolerance;

  scalar a = p.a, b = p.b;
#if EMBED
  if (!p.embed_flux && a.boundary[embed] != symmetry)
    p.embed_flux = embed_flux;
#endif // EMBED
  mgstats s = mg_solve_weighted ({a}, {b}, residual2, relax2,
      &p, p.nrelax, p.res, minlevel = max(1, p.minlevel), rho = rho);

  /**
  We restore the default. */

  if (p.tolerance)
    TOLERANCE = defaultol;

  return s;
}

struct Project_weighted {
  face vector uf;
  scalar p;
  face vector alpha;
  scalar rho;
  scalar f;
  bool liq;
  double dt;         // optional: default one
  int nrelax;        // optional: default four
};


/**
## Fields

We initialize the velocity potential field *ps*, the stefan velocity *ufs*,
and the divergence-free extended velocity *ufext*. */

scalar ps[], subdomain[];
face vector ufs[], ufext[];
mgstats mgpsf;

/**
The volume expansion term is declared in
[evaporation.h](/sandbox/ecipriano/src/evaporation.h). */

extern scalar stefanflow;

/**
## Helper functions

We define the function that performs the projection
of the stefan velocity onto a field with divergence
equal to the volume expansion term. */

face vector alpha_subdomain[];
scalar rho_subdomain[];

trace
mgstats project_sv (face vector uf, scalar p,
    (const) face vector alpha = unityf,
    double dt = 1.,
    int nrelax = 4)
{
  face vector ufs = q.uf;
  scalar ps = q.p;
  
  (const) face vector alpha = q.alpha.x.i ? q.alpha : unityf;
  (const) scalar rho = q.rho.i ? q.rho : unity;
  (const) scalar f = q.f;
  bool liq = q.liq;
  double dt = q.dt ? q.dt : 1.;
  int nrelax = q.nrelax ? q.nrelax : 4;


  foreach_dimension(){
    alpha_subdomain.x.restriction = alpha.x.restriction;
    #if TREE
    alpha_subdomain.x.refine = alpha.x.refine;
    #endif
    alpha_subdomain.x.prolongation = alpha.x.prolongation;
  }

  rho_subdomain.restriction = rho.restriction;
  rho_subdomain.prolongation = rho.prolongation;
  #if TREE
  rho_subdomain.refine = rho.refine;
  subdomain.refine = refine_injection;
  #endif
  subdomain.prolongation = refine_injection;
  

  foreach() {
    int count = 0;
    foreach_neighbor(2){
      if(f[] > F_ERR && f[] < 1-F_ERR)
        count++;
    }
    if (liq)
      subdomain[] = (count > 0 && f[] < 1-F_ERR) ? 1 : 0;
    else
      subdomain[] = (count > 0 && f[] > F_ERR) ? 1 : 0;
  }
  boundary({subdomain});


  foreach(){
    rho_subdomain[] = rho[];
    if (liq){
      if (subdomain[] == 0.0 && f[] > 1-F_ERR){
        rho_subdomain[] = 1e32;
        ps[] = 0.0;
      }
      if (subdomain[] == 0.0 && f[] < F_ERR){
        rho_subdomain[] = rho[]/2;
        ps[] = 0.0; 
      }
    }
    else{

       if (subdomain[] == 0.0 && f[] < 1-F_ERR){
        rho_subdomain[] = 1e32;
        ps[] = 0.0;
      }
      if (subdomain[] == 0.0 && f[] > F_ERR){
        rho_subdomain[] = rho[]/2;
        ps[] = 0.0; 
      }

    }

  }

  boundary((scalar *) {rho_subdomain, ps});

  foreach_face(){
    if (liq){
      alpha_subdomain.x[] = alpha.x[];
      if(subdomain[-1]*subdomain[0] == 0.0 && clamp(f[],0,1)*clamp(f[-1],0,1) > 0.0)
        alpha_subdomain.x[] = 0.0;
      if(subdomain[-1]*subdomain[0] == 0.0 && clamp(f[],0,1)*clamp(f[-1],0,1) <= 0.0)
        alpha_subdomain.x[] = 2*alpha.x[];
      ufs.x[] = 0.0;
    }
    else{

      alpha_subdomain.x[] = alpha.x[];
      if(subdomain[-1]*subdomain[0] == 0.0 && clamp(f[],0,1)*clamp(f[-1],0,1) <= 0.0)
        alpha_subdomain.x[] = 0.0;
      if(subdomain[-1]*subdomain[0] == 0.0 && clamp(f[],0,1)*clamp(f[-1],0,1) > 0.0)
        alpha_subdomain.x[] = 2*alpha.x[];
      ufs.x[] = 0.0;

    }
  }

  boundary((scalar *) {alpha_subdomain});

  scalar div[];
  foreach()
    div[] = stefanflow[]/dt;


  mgstats mgp = poisson_weighted (ps, div, alpha_subdomain,
       tolerance = TOLERANCE/sq(dt), nrelax = nrelax, rho = rho_subdomain, subdomain=subdomain);

  foreach_face()
    ufs.x[] = -dt*alpha_subdomain.x[]*face_gradient_x (ps, 0);
  boundary((scalar*){ufs});

  return mgp;
}

/**
## Boundary conditions

It is important to impose for *ps* the same boundary conditions and the
same tolerance used in the Poisson equation for the pressure *p*. */

ps[right] = neumann (neumann_pressure(ghost));
ps[left]  = neumann (- neumann_pressure(0));

#if AXI
ufs.n[bottom] = 0.;
ufs.t[bottom] = dirichlet(0);
ps[top]    = neumann (neumann_pressure(ghost));
#else // !AXI
#  if dimension > 1
ps[top]    = neumann (neumann_pressure(ghost));
ps[bottom] = neumann (- neumann_pressure(0));
#  endif
#  if dimension > 2
ps[front]  = neumann (neumann_pressure(ghost));
ps[back]   = neumann (- neumann_pressure(0));
#  endif
#endif // !AXI

/**
## Extended velocity

We perform the projection of the stefan velocity
by overloading the event end_timestep, defined in
[navier-stokes/centered.h](/src/navier-stokes/centered.h). */

event end_timestep (i++, last)
{

  /**
  We solve the Poisson equation using the multigrid solver. */

  mgpsf = project_sv (ufs, ps, alpha, rho, f, true, dt, mgpsf.nrelax);
  boundary ((scalar *){ufs});
  fprintf(ferr, "%d %d\n", mgpsf.i, mgpsf.nrelax); 

  /**
  We compute a divergence-free extended velocity by subtracting
  the stefan velocity from the field velocity. */

  foreach_face()
    ufext.x[] = uf.x[] - ufs.x[];
  boundary ((scalar *){ufext});

  double div_test1, div_test2, div_test3, fsum;

  div_test1 = 0.;
  div_test2 = 0.;
  div_test3 = 0.;
  fsum = 0.0;


  foreach(reduction(+:div_test1), reduction(+:div_test2), reduction(+:div_test3), reduction(+:fsum)) {

    foreach_dimension(){
      div_test1 += (uf.x[1] - uf.x[])/dt*dv();
      div_test2 += (ufs.x[1] - ufs.x[])/dt*dv();
      div_test3 += (ufext.x[1] - ufext.x[])/dt*dv();
    }
    fsum += (1-f[])*dv();
  }

  fprintf(fout, "%12e %12e %12e %12e %12e\n", t, div_test1, div_test2, div_test3, fsum); 

}

/**
## Notes
This approach works fine when the field velocity is larger than
the velocity due to the phase change. An example is a falling
droplet, or a droplet in forced convective conditions. Static
droplets, evaporating in reduced gravity conditions, can be
simulated using this method just for small vaporization rates
or for small density ratio values. If these conditions are not
met, the double pressure velocity coupling approach should be
preferred.

## References

~~~bib
@article{scapin2020volume,
  title={A volume-of-fluid method for interface-resolved simulations of phase-changing two-fluid flows},
  author={Scapin, Nicol{\`o} and Costa, Pedro and Brandt, Luca},
  journal={Journal of Computational Physics},
  volume={407},
  pages={109251},
  year={2020},
  publisher={Elsevier}
}

@article{malan2021geometric,
  title={A geometric VOF method for interface resolved phase change and conservative thermal energy advection},
  author={Malan, LC and Malan, Arnaud G and Zaleski, St{\'e}phane and Rousseau, PG},
  journal={Journal of Computational Physics},
  volume={426},
  pages={109920},
  year={2021},
  publisher={Elsevier}
}

@article{palmore2019volume,
  title={A volume of fluid framework for interface-resolved simulations of vaporizing liquid-gas flows},
  author={Palmore Jr, John and Desjardins, Olivier},
  journal={Journal of Computational Physics},
  volume={399},
  pages={108954},
  year={2019},
  publisher={Elsevier}
}
~~~

*/


