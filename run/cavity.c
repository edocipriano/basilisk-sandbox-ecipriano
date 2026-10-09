/**
# Natural Convection in a Differentially Heated Cavity

A square cavity filled with air has its vertical walls at the temperatures
$T_h = T_0(1 + \epsilon)$ and $T_c = T_0(1 - \epsilon)$, with $\epsilon = 0.6$,
and adiabatic horizontal walls. The temperature difference is too large for the
Boussinesq approximation: the density follows the ideal gas law, the transport
properties follow Sutherland's law, and, since the cavity is closed, the
thermodynamic pressure is set by the mass of gas. This is the benchmark of
[Le Quéré et al., 2005](#lequere2005modelling), which provides the steady
Nusselt number and thermodynamic pressure to six digits.

Compared to [pressurization.c](pressurization.c), this case has no analytic
solution, but it couples the pressurization model of
[low-mach.h](../src/navier-stokes/low-mach.h) with a buoyancy-driven flow.

![Steady temperature field on the finest grid](cavity/temperature.png)
*/

/**
## Simulation Setup

The cavity contains only the gas phase: the phase change model is used with a
null vaporization rate and a null volume fraction. The boundary layers cover
the whole height of the walls, therefore we use uniform grids. */

#include "grid/multigrid.h"
#include "navier-stokes/low-mach.h"
#define P_ERR 1.e-6
#include "two-phase-varprop.h"
#include "basilisk-properties.h"
#include "two-phase.h"
#include "gravity.h"
#include "phasechange.h"
#include "fixedflux.h"
#include "view.h"

/**
### Case Parameters

The benchmark proposes three cases: `T1` with constant transport properties,
`T2` (default) and `T3` with Sutherland's laws. For a matter of computational
time, the test runs only the grids from $32^2$ to $128^2$: the convergence
towards the reference solution can be observed by increasing the level of
refinement (e.g. `-DLEVELMAX=9`). */

#ifndef CASE
# define CASE 2
#endif

#if CASE == 1
# define RAYLEIGH 1.e6
# define SUTHERLAND 0
#elif CASE == 2
# define RAYLEIGH 1.e6
# define SUTHERLAND 1
#else
# define RAYLEIGH 1.e7
# define SUTHERLAND 1
#endif

#ifndef LEVELMIN
# define LEVELMIN 5
#endif
#ifndef LEVELMAX
# define LEVELMAX 7
#endif

#define P00     101325.   // Initial thermodynamic pressure [Pa]
#define RGAS    287.      // Specific gas constant [J/kg/K]
#define PRANDTL 0.71      // Prandtl number [-]
#define GAMMA   1.4       // Heat capacity ratio [-]
#define EPS     0.6       // Non-dimensional temperature difference [-]

int maxlevel;
double Thot, Tcold, T0ref, mu0, lambda0, alpha0, tend;

/**
### Material Properties

The viscosity follows Sutherland's law, and the thermal conductivity is given
by the constant Prandtl number, $\lambda = \mu c_p/Pr$. */

double sutherland (double T) {
  return 1.68e-5*pow (T/273., 1.5)*(273. + 110.5)/(T + 110.5);
}

double mu_law (double T) {
  return SUTHERLAND ? sutherland (T) : mu0;
}

double lambda_law (double T) {
  return mu_law (T)*cp2/PRANDTL;
}

double gasprop_density_cavity (void * p) {
  ThermoState * ts = p;
  return ts->P/(RGAS*ts->T);
}

double gasprop_viscosity_cavity (void * p) {
  ThermoState * ts = p;
  return mu_law (ts->T);
}

double gasprop_conductivity_cavity (void * p) {
  ThermoState * ts = p;
  return lambda_law (ts->T);
}

/**
### Boundary Conditions

The walls are no-slip. The vertical walls are isothermal, the horizontal ones
are adiabatic (default). */

u.t[left]   = dirichlet (0.);
u.t[right]  = dirichlet (0.);
u.t[top]    = dirichlet (0.);
u.t[bottom] = dirichlet (0.);

T[left]  = dirichlet (Thot);
T[right] = dirichlet (Tcold);

int main (void) {
  NGS = 1, NLS = 1;

  T0ref = 600.;
  Thot = T0ref*(1. + EPS), Tcold = T0ref*(1. - EPS);
  P0 = P00, TG0 = T0ref, TL0 = T0ref;

  cp1 = cp2 = GAMMA*RGAS/(GAMMA - 1.);
  mu0 = sutherland (T0ref);
  lambda0 = lambda_law (T0ref);
  double rho0 = P00/(RGAS*T0ref);
  alpha0 = lambda0/(rho0*cp2);

  /**
  The liquid phase is never used, but its properties must be finite. */

  rho1 = 1000., mu1 = 1.e-3, lambda1 = 0., Dmix1 = 0., dhev = 0.;
  rho2 = rho0, mu2 = mu0, lambda2 = lambda0, Dmix2 = 0.;

  closed = true;
  pcm.isomassfrac = true;
  pcm.divergence = true;
  mEvapVal = 0.;

  /**
  The size of the cavity follows from the Rayleigh number
  $Ra = Pr\,g\rho_0^2(T_h - T_c)L^3/(T_0\mu_0^2)$. The maximum time step only
  limits the first steps, then the CFL condition takes over. The Nusselt number
  and the thermodynamic pressure are steady to five digits after
  $0.25L^2/\alpha_0$, and we stop the simulation at $0.3L^2/\alpha_0$. */

  G.y = -9.81;
  L0 = cbrt (RAYLEIGH*T0ref*sq(mu0) /
      (PRANDTL*fabs (G.y)*sq(rho0)*(Thot - Tcold)));
  DT = 0.1*L0/sqrt (fabs (G.y)*L0*(Thot - Tcold)/T0ref);
  tend = 0.3*sq(L0)/alpha0;

  for (maxlevel = LEVELMIN; maxlevel <= LEVELMAX; maxlevel++) {
    init_grid (1 << maxlevel);
    run();
  }
}

/**
### Thermodynamic Pressure

The thermodynamic pressure is computed from the conservation of the mass of
gas, i.e. of the gas initially at rest at $T_0$ and $P_{00}$:
$$
  P_0 = P_{00}\dfrac{\int_\Omega dV/T_0}{\int_\Omega dV/T}
$$
Integrating $dP_0/dt$ instead accumulates the residual imbalance between the
discrete heat fluxes of the two walls, and the pressure drifts away during the
thousands of time steps needed to reach the steady state. */

void pressure_from_mass (void) {
  scalar TG = gas->T;
  double invT = 0.;
  foreach (reduction(+:invT))
    invT += dv()/TG[];
  P0 = P00*sq(L0)/(T0ref*invT);
}

event pressurization (i++) {
  pressure_from_mass();
}

/**
### Initial Conditions

Each level starts from the gas at rest at $T_0$. The conductivity at the walls
is set to its value at the wall temperature: otherwise, with Sutherland's law,
the wall flux uses the conductivity of the first cell, which unbalances the
heat fluxes of the two walls (35% instead of 4% at level 5). */

event init (i = 0) {
  ThermoState tsl, tsg;
  tsl.T = T0ref, tsl.P = P0, tsl.x = (double[]){1.};
  tsg.T = T0ref, tsg.P = P0, tsg.x = (double[]){1.};
  phase_set_thermo_state (liq, &tsl);
  phase_set_thermo_state (gas, &tsg);
  phase_set_properties (liq, MWs = (double[]){R_GAS*1.e3/RGAS});
  phase_set_properties (gas, MWs = (double[]){R_GAS*1.e3/RGAS});

  tp2.rhov    = gasprop_density_cavity;
  tp2.muv     = gasprop_viscosity_cavity;
  tp2.lambdav = gasprop_conductivity_cavity;
  tp2.betaT   = gasprop_thermal_expansion;
  tp2.chiT    = gasprop_isothermal_compressibility;

  scalar TG = gas->T, lambdaG = gas->lambda;
  TG[left]  = dirichlet (Thot);
  TG[right] = dirichlet (Tcold);
  lambdaG[left]  = dirichlet (lambda_law (Thot));
  lambdaG[right] = dirichlet (lambda_law (Tcold));
}

/**
## Post-Processing

The average Nusselt number of a vertical wall is
$$
  \overline{Nu} = \dfrac{1}{\lambda_0\left(T_h-T_c\right)}
    \int_0^L \lambda\dfrac{\partial T}{\partial x}\bigg|_w dy
$$
computed with the face gradient between the first cell and its ghost value,
i.e. with the wall flux of the discrete equations. */

double nusselt (bool hot) {
  scalar TG = gas->T;
  boundary ({TG});
  double flux = 0.;
  if (hot)
    foreach_boundary (left, reduction(+:flux))
      flux += lambda_law (Thot)*(TG[] - TG[-1]);
  else
    foreach_boundary (right, reduction(+:flux))
      flux += lambda_law (Tcold)*(TG[1] - TG[]);
  return fabs (flux)/(lambda0*(Thot - Tcold));
}

/**
We write the Nusselt number of the two walls and the pressure ratio
$P_0/P_{00}$ at fixed fractions of the simulation time. */

event output_data (i++) {
  static FILE * fp = NULL;
  static double tout;
  if (i == 0) {
    char name[80];
    sprintf (name, "OutputData-%d", maxlevel);
    if (fp)
      fclose (fp);
    fp = fopen (name, "w");
    tout = 0.;
  }
  if (t < tout)
    return 0;
  tout += tend/100.;

  fprintf (fp, "%g %g %g %g\n", t*alpha0/sq(L0), nusselt (true),
      nusselt (false), P0/P00);
  fflush (fp);
}

event stop (t = tend);

/**
We write the steady values (for testing) and the temperature field of the
finest grid. The pressure is computed again because
[low-mach.h](../src/navier-stokes/low-mach.h) resets it at the end of the
run. */

event logger (t = end) {
  pressure_from_mass();
  double nuh = nusselt (true), nuc = nusselt (false);
  fprintf (stderr, "%d %.3f %.3f %.5f\n", maxlevel, nuh, nuc, P0/P00);

  if (maxlevel == LEVELMAX) {
    clear();
    box();
    view (tx = -0.5, ty = -0.5);
    squares ("T", min = Tcold, max = Thot, linear = true,
        map = blue_white_red);
    isoline ("T", n = 25, lw = 2.);
    save ("temperature.png");
  }
}

/**
## Results

The steady values of the default case `T2` are compared with the reference
solution $\overline{Nu} = 8.6866$, $P_0/P_{00} = 0.924487$:

| Grid | $\overline{Nu}_h$ | $\overline{Nu}_c$ | $P_0/P_{00}$ | error $\overline{Nu}_h$ | error $P_0$ |
|------|-------|-------|---------|--------|--------|
| $32^2$  | 8.401 | 8.819 | 0.97966 | -3.29% | +5.97% |
| $64^2$  | 8.681 | 8.926 | 0.93522 | -0.06% | +1.16% |
| $128^2$ | 8.670 | 8.764 | 0.92658 | -0.19% | +0.23% |

The thermodynamic pressure converges monotonically towards the reference value.
The convergence of the Nusselt number is instead not monotonic: the error of the
hot wall almost vanishes on the $64^2$ grid and it increases again on the
$128^2$ grid, while the error of the cold wall increases from the $32^2$ to the
$64^2$ grid. The heat fluxes of the two walls, which must be equal at steady
state, differ by 5%, 2.8% and 1.1% on the three grids. This energy conservation
error is a known feature of segregated algorithms that advance the temperature
in non-conservative form, with the density computed from the equation of state
(algorithm A1 of [Knikker, 2011](#knikker2011comparative)). It is added to the
discretization error of the Nusselt number with a different sign, resulting in
a non-monotonic convergence on coarse grids. As in Knikker (2011), the cold
wall, whose boundary layer is thinner, is the most affected. Energy conserving
algorithms remove this error.

The steady state also depends on the time step: the BCG scheme of
[vof.h](/src/vof.h) adds a numerical diffusion $u^2\Delta t/2$ along the
streamlines, which does not vanish at steady state.

~~~gnuplot Convergence of the Nusselt number towards the steady state
reset
set grid
set key bottom right
set xlabel "t α_0/L^2 [-]"
set ylabel "Nu [-]"
set yrange [7:10]
set size square
set label "hot wall" at 0.12,7.55 center
set arrow from 0.12,7.65 to 0.12,8.33 head filled
set label "cold wall" at 0.12,9.65 center
set arrow from 0.12,9.55 to 0.12,9.0 head filled
levels = system("awk '{print $1}' log")
plot for [l in levels] "OutputData-".l u 1:2 w l lw 2 lc (l - 4) \
       t sprintf("%d^2", 2**(l + 0)), \
     for [l in levels] "OutputData-".l u 1:3 w l lw 2 lc (l - 4) \
       dt 2 notitle, \
     8.6866 w l lc rgb "black" dt 2 t "Le Quéré et al. 2005"
~~~

~~~gnuplot Convergence of the thermodynamic pressure towards the steady state
reset
set grid
set key top right
set xlabel "t α_0/L^2 [-]"
set ylabel "P_0/P_{00} [-]"
set size square
levels = system("awk '{print $1}' log")
plot for [l in levels] "OutputData-".l u 1:4 w l lw 2 lc (l - 4) \
       t sprintf("%d^2", 2**(l + 0)), \
     0.924487 w l lc rgb "black" dt 2 t "Le Quéré et al. 2005"
~~~

## References

~~~bib
@article{lequere2005modelling,
  title={Modelling of natural convection flows with large temperature
         differences: a benchmark problem for low Mach number solvers.
         Part 1. Reference solutions},
  author={Le Qu{\'e}r{\'e}, Patrick and Weisman, Catherine and
          Paill{\`e}re, Henri and Vierendeels, Jan and Dick, Erik and
          Becker, Roland and Braack, Malte and Locke, James},
  journal={ESAIM: Mathematical Modelling and Numerical Analysis},
  volume={39},
  number={3},
  pages={609--616},
  year={2005},
  publisher={EDP Sciences}
}

@article{knikker2011comparative,
  title={A comparative study of high-order variable-property segregated
         algorithms for unsteady low Mach number flows},
  author={Knikker, Ronnie},
  journal={International Journal for Numerical Methods in Fluids},
  volume={66},
  number={4},
  pages={403--427},
  year={2011},
  doi={10.1002/fld.2261}
}
~~~
*/
