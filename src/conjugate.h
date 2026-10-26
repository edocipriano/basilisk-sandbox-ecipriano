/**
# Conjugate Heat Transfer

This module implements Conjugate Heat Transfer (CHT) in a naive manner, which
does not require embedded boundaries. One of the boundaries of the (square)
domain is considered as the solid-fluid interface; the solid is resolved on the
same grid as the fluid, and the boundary conditions on the solid-fluid interface
are computed and set in order to achieve the heat exchange between the different
phases.

The main advantage of this method is its simplicity, and the fact that it can
exploit the default ghost cells method to impose the boundary conditions.
Therefore, it does not require dealing with interfacial gradients as tipically
done with the embedded boundaries. However, the cost is the use of simple
geometries, and the introduction of additional cells into the solution of the
fluid phase. Conversely from methods using a master-slave coupling approach,
this approach automatically guarantees MPI parallelization, and it does not
introduce the overhead due to the interpolation of values from another domain,
potentially in another processor. In many systems (e.g., droplet combustion or
nucleate boiling problems) the solid heater is simple and can be effectively
described by this approach.

## Usage

Include this module after the Navier--Stokes and two-phase solver, there is not
strict rule about it. Set the global variables which store the solid material
properties.  Set additional heat sources such as the joule effect. Other source
terms can be introduced by directly modifying the field `QSexp` introduced at
the next section.  Finally, decide the shape of the solid by setting the volume
and surface fractions `cw` and `fw`, this can be done using the `solid()`
function but consider that the solid-fluid boundary must coicinde with one of the
domain boundaries. The solid-fluid boundary is the left one by default, or the
bottom one with `CONJUGATE_BOTTOM`. The boundary conditions of the phase
temperatures on the solid-fluid boundary are set by this module. No boundary
condition is needed on `TS` (the default symmetry condition is used, the heat
exchanged with the fluid is a source term).

The properties of the solid are stored in the fields `rhos`, `cps`, and
`lambdas`, initialized with the uniform values `rho3`, `cp3`, and `lambda3`.
They can be modified in the simulation file to describe non-uniform solids
(e.g. a thin heater film on the substrate).

~~~literatec
#include "conjugate.h"

int main() {
  ...
  rho3 = 8000., cp3 = 1080., lambda3 = 0.2;
  qjoule = 1e4;
  ...
}

event init (i = 0) {
  solid (cw, fw, x - LW);
}
~~~

## Improvements

1. Consider using a better adaptation method for `TS`.
1. Improve handling of Dirichlet boundary condition on the solid (embedded)
boundary which is not in contact with a domain boundary. */

/**
## Variables

We introduce the *wall* metric factors `cw` and `fw` which are used to correct
the integration of the solid heat equation. */

scalar cw[];
face vector fw[];

/**
The solid temperature `TS` is resolved in this file. The heat source/sink terms
can be introduced using the field `QSexp` [W/m$^3$]. The field `qw` stores the
heat flux [W/m$^2$] from the solid to the fluid through the wall face of the
cells adjacent to the solid-fluid boundary (zero elsewhere). */

scalar TS[], QSexp[], qw[];

/**
### Physical properties

We define the density, heat capacity, and thermal conductivity of the solid.
The fields are initialized with the uniform values `rho3`, `cp3`, and
`lambda3`, which are also used as reference values to scale the solid heat
equation. The thermal power of the Joule effect can be set using the respective
variable, initialized as null. The initial solid temperature is set through
`TS0`. */

scalar rhos[], cps[], lambdas[];

double rho3 = 1., cp3 = 1., lambda3 = 1., qjoule = 0., TS0 = 300.;

/**
### Heat transfer resistance

We introduce the values of heat transfer resistance (from Hertz-Knudsen theory) and
the accommodation coefficient. The former is a function of the material
properties, while the latter is a case-dependent factor. */

double RInt = 0., acc = 1.;

/**
### Solid-fluid boundary

The normal direction `cwn` of the solid-fluid boundary, and the condition
identifying the cells adjacent to it. */

#ifndef CONJUGATE_BOTTOM
# define CONJUGATE_BOTTOM 0
#endif

#if CONJUGATE_BOTTOM
# define cwn y
# define conjugate_wall_cell() (y - Y0 < Delta)
#else
# define cwn x
# define conjugate_wall_cell() (x - X0 < Delta)
#endif

/**
## Initialization

We specify an initial solid temperature and the function used to refine or
coarsen the solid temperature and the solid volume fraction. */

event defaults (i = 0) {
  foreach() {
    TS[] = TS0;
    qw[] = 0.;
    rhos[] = rho3, cps[] = cp3, lambdas[] = lambda3;
  }
#if TREE
  cw.refine = fraction_refine;
  set_prolongation (cw, fraction_refine);
  for (scalar t in {TS}) {
    // fixme: solid tracers deserve a better refinement method
    t.refine = refine_injection;
    set_prolongation (t, refine_injection);
    //t.c = cw;
    //t.depends = list_add (t.depends, cw);
  }
#endif
}

event reset_sources (i++) {
  foreach()
    QSexp[] = 0.;
}

/**
## Solid-Fluid coupling

The phase temperatures and the solid temperature are solved sequentially: first
the fluid, using the solid temperature of the previous time step, then the
solid, which receives the heat flux absorbed by the fluid. The heat flux
leaving the solid is equal to the one entering the fluid, therefore the energy
is conserved at the solid-fluid boundary.

### Fluid side

The boundary temperature of each phase is computed from the balance of the
heat fluxes, considering the half-cell resistance $R_f = \Delta/(2\lambda)$ of
the fluid cell, and the external resistance $R_{ext}$ (the contact resistance
plus the half-cell resistance of the solid):
$$
T_B = \dfrac{R_f T_S + R_{ext} T}{R_f + R_{ext}}
\qquad
R_{ext} = acc\,R_{Int} + \dfrac{\Delta}{2\lambda_s}
$$
This is a Robin condition, $T_B + \lambda R_{ext}\,\partial_n T = T_S$, imposed
using [navier()](/src/grid/cartesian-common.h#navier). It is implicit with
respect to the phase temperature $T$, while the solid temperature $T_S$ is the
one of the previous time step. Where the wall face is not in contact with the
solid ($f_w = 0$) the boundary is adiabatic. */

static double conjugate_rext (Point point) {
  return acc*RInt + Delta/(2.*lambdas[]);
}

double conjugate_length (Phase * phase, Point point = point) {
  scalar lambda = phase->lambda;
  return lambda[]*conjugate_rext (point);
}

event init (i = 0) {
  scalar TL = liq->T, TG = gas->T;
#if CONJUGATE_BOTTOM
  TL[bottom] = fw.y[] > 0. ? navier (TS[], conjugate_length (liq)) : neumann (0.);
  TG[bottom] = fw.y[] > 0. ? navier (TS[], conjugate_length (gas)) : neumann (0.);
#else
  TL[left] = fw.x[] > 0. ? navier (TS[], conjugate_length (liq)) : neumann (0.);
  TG[left] = fw.x[] > 0. ? navier (TS[], conjugate_length (gas)) : neumann (0.);
#endif
}

/**
The wall metric `fw` is updated before the solution of the phase temperatures.
The fields used by the boundary conditions (including the conductivities of the
phases) are restricted, since the conditions are applied also on the coarse
levels of the multigrid solver. Otherwise the coarse levels would see a null
Robin length (i.e. a Dirichlet condition) and the solver could diverge. */

event tracer_diffusion (i++) {
  face_fraction (cw, fw);
  scalar lambdal = liq->lambda, lambdag = gas->lambda;
  restriction ((scalar *){TS,lambdas,fw,lambdal,lambdag});
}

/**
### Wall heat flux

After the solution of the phase temperatures, we compute the heat flux absorbed
by each phase through the wall face, using the same face fractions of
[phase_diffusion()](phase.h#phase_diffusion) and the same boundary condition:
$$
q_w = \sum_{k = l,g} f_{s,k}\dfrac{T_S - T_k}{R_{f,k} + R_{ext}}
$$
This is done in the `properties` event, before the phase temperatures are
modified by other modules. */

static void conjugate_face_fraction (Phase * phase, face vector fs) {
  scalar ff[];
  scalar_clone (ff, f);
  ff.inverse = phase->inverse;
  foreach()
    ff[] = phase->inverse ? 1. - f[] : f[];
  face_fraction (ff, fs);
}

static double conjugate_phase_flux (Phase * phase, face vector fs,
    Point point)
{
  scalar T = phase->T, lambda = phase->lambda;
  double ff = phase->inverse ? 1. - f[] : f[];
  if (ff <= 0. || phase->isothermal)
    return 0.;
  double Rf = Delta/(2.*lambda[]);
  return fs.cwn[]*(TS[] - T[]/ff)/(Rf + conjugate_rext (point));
}

event properties (i++) {
  face vector fsl[], fsg[];
  conjugate_face_fraction (liq, fsl);
  conjugate_face_fraction (gas, fsg);
  foreach() {
    qw[] = 0.;
    if (conjugate_wall_cell() && fw.cwn[] > 0.)
      qw[] = conjugate_phase_flux (liq, fsl, point) +
        conjugate_phase_flux (gas, fsg, point);
  }
}

/**
The temperatures at the solid and at the fluid side of the wall face, which
differ because of the contact resistance, can be computed for post-processing
from the wall heat flux. */

double conjugate_wall_solid (Point point) {
  return TS[] - qw[]*Delta/(2.*lambdas[]);
}

double conjugate_wall_fluid (Point point) {
  return conjugate_wall_solid (point) - acc*RInt*qw[];
}

/**
### Dirichlet condition on the solid boundary

We use the interface gradients calculation in order to impose Dirichlet BCs on
the solid heater. */

double conjugate_dirichlet = 0.;

static void conjugate (void) {
  if (conjugate_dirichlet) {
    scalar fl[], fg[];
    foreach() {
      fl[] = cw[];
      fg[] = 1. - cw[];
    }
    fg.inverse = true;

    face vector fsl[], fsg[];
    face_fraction (fl, fsl);
    foreach_face()
      fsg.x[] = fsl.x[];

    foreach_interfacial_plic (cw, F_ERR) {
      double bc = conjugate_dirichlet;
      double intgrad = plic_gradient (point, TS, fl, fsl, bc, true, NULL);
      QSexp[] += lambda3*intgrad*dirac;
    }
  }
}

/**
## Heat equation

We solve the heat equation for the solid phase:

$$
\rho_s cp_s \dfrac{\partial T_s}{\partial t} =
\nabla\cdot\left(\lambda_s\nabla T_s\right) + \dot{Q}_s
$$

considering metric factors and heat source/sink terms. Explicit (`e`) and
implicit (`i`) source terms to the solid diffusion equation can be added
exploiting the function pointer `energy_sources` defined below. */

void no_sources (scalar i, scalar e) {
  return;
}

void (* energy_sources) (scalar i, scalar e) = no_sources;

/**
The solid temperature is solved at the end of the time step, after the phase
temperatures. The heat flux $q_w$ leaving the solid through the wall face is a
sink term in the cells adjacent to the solid-fluid boundary, with the face
metric used by the fluid. The conductivity on the faces is the harmonic mean of
the cell values (layered solids). The equation is divided by the reference
value `rho3*cp3`. Otherwise, the residual would be of order
$\rho_s c_{p,s} T/\Delta t$, and the multigrid solver might not reach
`TOLERANCE` because of round-off errors. */

event end_timestep (i++) {
  conjugate();

  scalar cmm[];
  foreach()
#if AXI
    cmm[] = max (fabs (Y0 - (y - Y0)), 1./HUGE)*cw[];
#else
    cmm[] = cm[]*cw[];
#endif

  face vector fmm[];
  foreach_face()
#if AXI
    fmm.x[] = max (fabs (Y0 - (y - Y0)), 1./HUGE)*fw.x[];
#else
    fmm.x[] = fm.x[]*fw.x[];
#endif

  double rhocp0 = rho3*cp3;

  scalar theta[];
  foreach()
    theta[] = max (cmm[]*rhos[]*cps[]/rhocp0, 1e-10);

  face vector D[];
  foreach_face() {
    double lambdaf = 2.*lambdas[]*lambdas[-1]/(lambdas[] + lambdas[-1] + 1e-30);
    D.x[] = fmm.x[]*lambdaf/rhocp0;
  }

  energy_sources (zeroc, QSexp);

  scalar r[];
  foreach() {
    r[] = QSexp[];
    if (conjugate_wall_cell())
      r[] -= qw[]*fm.cwn[]/Delta;
    r[] /= rhocp0;
  }

  diffusion (TS, dt, D = D, theta = theta, r = r);
}

/**
## Solid properties

The materials of the solid heaters/suspenders are pretty much the same in all
experimental works. Therefore, for common materials we provide a functions that
helps initializing solid properties. */

enum solid_props {
  QUARTZ, SiC, SAPPHIRE
};

void solid_properties (enum solid_props s) {
  switch (s) {
    case QUARTZ:
      rho3 = 2220, cp3 = 760., lambda3 = 1.5;
      break;
    case SiC:
      rho3 = 2740., cp3 = 670., lambda3 = 5.2;
      break;
    case SAPPHIRE:
      rho3 = 4890., cp3 = 410., lambda3 = 10.9;
      break;
  }
}
