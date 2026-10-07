/**
# Passive Pressurization of a Closed Tank

A closed tank, partially filled with liquid, receives a uniform heat flux
$q_w$ from the environment through its walls. This configuration reproduces
the passive pressurization of cryogenic tanks (benchmark case 2.3 of the HASTA
project), and it extends the test case [pressurization.c](pressurization.c): the
heat is not supplied by a volumetric source, but it enters the domain by
conduction across the walls of the tank.

The thermal conductivity is not zero, therefore the heat flux supplied by the
walls is transported by the Fourier law in the gas and in the liquid. The
temperature of the ullage is not uniform, and the divergence of the conductive
flux makes the velocity divergence non-null. This volumetric expansion is
balanced by the compression due to the increase of the thermodynamic pressure
$P_0$. The phase change is switched off and the interface is adiabatic,
therefore the energy of the ullage changes only because of the heat $Q_g$
supplied by the walls wetted by the gas. The liquid is incompressible, the
volume of the ullage $V_g$ is constant, and integrating the low Mach number
energy equation over the ullage gives:
$$
  \dfrac{dP_0}{dt} = \left(\gamma - 1\right)\dfrac{Q_g}{V_g}
$$
This solution does not depend on the temperature distribution in the ullage.
The total mass of the ullage must be conserved, while the gas density changes
locally with the temperature.

The configuration is computed in two tanks. In the square box of side $L_0$,
the heat flux is imposed on the boundaries of the domain. The interface is at
a height $H_l$ from the bottom, the ullage has height $H_g = L_0 - H_l$, and
(per unit depth) $Q_g = q_w \left(L_0 + 2H_g\right)$ and $V_g = L_0 H_g$:
$$
  P_0(t) = P_0^0 + \left(\gamma - 1\right)
    \dfrac{q_w \left(L_0 + 2H_g\right)}{L_0 H_g}t
$$
The heat flux through the boundaries of the domain is included in the face
fluxes, and the interface meets the walls at a right angle. This case
therefore verifies the pressurization model independently of the treatment of
the embedded boundaries.

The second tank is a circle of radius $R$ obtained with an embedded boundary
(`-DSOLID=1`, [heatflux-embed.c](heatflux-embed.c)). The interface is at a
distance $h$ above the center of the tank, and the ullage is a circular
segment. Defining $\theta = \arccos (h/R)$, the wall wetted by the gas has
length $2R\theta$ and $V_g = R^2\theta - h\sqrt{R^2 - h^2}$, therefore:
$$
  P_0(t) = P_0^0 + \left(\gamma - 1\right)
    \dfrac{2 q_w R\theta}{R^2\theta - h\sqrt{R^2 - h^2}}t
$$
which reduces to $P_0^0 + 2\left(\gamma - 1\right)q_w t/R$ for $h = 0$. The
heat flux is imposed as a Neumann condition on the embedded boundary. The flux
through the embedded boundary is computed by `phase_embed_flux()` (see
[phase.h](/sandbox/ecipriano/src/phase.h)), both in the temperature equation and
in the velocity divergence. In the cut cells where the interface meets the wall,
this function splits the flux between the phases according to the length of the
wall wetted by each phase (see the utilities for embedded boundaries in
[common-evaporation.h](/sandbox/ecipriano/src/common-evaporation.h)), while the
contact angle is not imposed. Comparing the two tanks shows which inaccuracies
are due to the embedded boundaries.

![Evolution of the temperature field (square box)](heatflux/movie.mp4)

![Evolution of the temperature field (embedded boundary)](heatflux-embed/movie.mp4)
*/

/**
## Simulation Setup

We use the low Mach number Navier--Stokes solver, which includes the volumetric
expansion terms in the projection step, together with the variable properties
module. The phase change is switched off using the fixed flux model with a null
vaporization rate, in order to focus exclusively on the pressurization.
*/

#include "grid/multigrid.h"
#if SOLID
# include "embed.h"
#endif
#include "navier-stokes/low-mach.h"
#define P_ERR 1.e-6
#include "two-phase-varprop.h"
#include "basilisk-properties.h"
#include "two-phase.h"
#include "tension.h"
#include "phasechange.h"
#include "fixedflux.h"
#include "view.h"

/**
### Model Data

The molecular weight of the gas phase, the heat flux supplied by the walls, the
radius of the circular tank, and the initial liquid level. */

int maxlevel, minlevel = 2;
double qw = 100., R0, level0 = 0.5;
double gamma0, mg0, P0exact;

#define MWG 29.
#define TEND 1.

/**
The density of the gas phase follows the ideal gas law, which makes the ullage
compressible and sensitive to the variation of the thermodynamic pressure. */

double gasprop_density_idealgas (void * p) {
  ThermoState * ts = p;
  return ts->P*MWG/(R_GAS*1.e3*ts->T);
}

/**
### Boundary Conditions

The heat flux is imposed in the `init` event, where the temperature fields of
the two phases are available: on the boundaries of the domain for the square
box, and on the embedded boundary for the circular tank. In the latter case,
the boundaries of the domain lie outside of the tank. */

int main (void) {

  /**
  We use a single chemical species in each phase. */

  NGS = 1, NLS = 1;

  /**
  We set the material properties of water and air. The species diffusivity is
  null, since the composition of the phases does not change. The gas phase
  density is overwritten by the ideal gas law. */

  rho1 = 1000., rho2 = 1.;
  mu1 = 1.e-3, mu2 = 1.8e-5;
  Dmix1 = 0., Dmix2 = 0.;
  lambda1 = 0.6, lambda2 = 0.026;
  cp1 = 4184., cp2 = 1004.5;
  dhev = 0.;

  P0 = 101325.;
  TG0 = 300., TL0 = 300.;

  /**
  The system is closed: the net expansion is converted into a variation of the
  thermodynamic pressure instead of leaving the domain. */

  closed = true;

  /**
  The composition of the two phases does not change, and the phase change is
  switched off setting a null vaporization rate. */

  pcm.isomassfrac = true;
  pcm.divergence = true;
  mEvapVal = 0.;

  L0 = 0.01;
  R0 = 0.45*L0;
  DT = 1.e-3;

  /**
  The heat capacity ratio of the ullage is used by the analytic solution. */

  gamma0 = cp2/(cp2 - R_GAS*1.e3/MWG);

  /**
  We reduce the tolerance of the Poisson solver: with the default value the
  velocity oscillates in time during the transient. The tolerance of the
  temperature equation is a residual, with the units of
  $\rho c_p T/\Delta t$ (see [phase.h](/sandbox/ecipriano/src/phase.h)): it
  corresponds to a change of the temperature of about $10^{-8}$ K per
  timestep in the ullage. */

  TOLERANCE = 1e-6 [*];
  TOLERANCE_TEMPERATURE = 1e-2;

  for (maxlevel = 4; maxlevel <= 8; maxlevel++) {
    init_grid (1 << maxlevel);
    run();
  }
}

/**
We initialize a flat interface at the height $H_l = 0.5 L_0$. In the embedded
case, the tank is a circle of radius $0.45 L_0$ centered in the box, and the
interface crosses its center ($h = 0$). The interface lies on a grid line on
all the grids, therefore the interfacial cells are initially full or empty.
When the interface cuts the cells (e.g. `level0 = 0.51`), the VOF advection
adds an $O(\Delta)$ error to the volume of the ullage, because the velocity
divergence of the compressible gas in the mixed cells is attributed to the
liquid when $f > 0.5$ (see [vof.h](/src/vof.h)): on the levels where this
happens (6 and 8 for `level0 = 0.51`), the mass of the ullage changes by up to
$10^{-4}$ and the pressurization rate changes accordingly, in both tanks. */

#define circle(x,y,R) (sq(R) - sq(x - 0.5*L0) - sq(y - 0.5*L0))

event init (i = 0) {
#if EMBED
  solid (cs, fs, circle (x, y, R0));
  fractions_cleanup (cs, fs);
#endif
  fraction (f, level0*L0 - y);

  ThermoState tsl, tsg;
  tsl.T = TL0, tsl.P = P0, tsl.x = (double[]){1.};
  tsg.T = TG0, tsg.P = P0, tsg.x = (double[]){1.};

  phase_set_thermo_state (liq, &tsl);
  phase_set_thermo_state (gas, &tsg);

  phase_set_properties (liq, MWs = (double[]){18.});
  phase_set_properties (gas, MWs = (double[]){MWG});

  /**
  The gas phase is an ideal gas: its density, thermal expansion coefficient and
  isothermal compressibility are set to the corresponding analytic functions. */

  tp2.rhov = gasprop_density_idealgas;
  tp2.betaT = gasprop_thermal_expansion;
  tp2.chiT = gasprop_isothermal_compressibility;

  /**
  The heat flux $q_w$ enters both phases through the walls of the tank. The
  boundary condition is the temperature gradient normal to the wall, which is
  positive towards the outside of the tank. The heat supplied to the liquid
  does not contribute to the pressurization, since the liquid is
  incompressible. */

  scalar TL = liq->T, TG = gas->T;
#if EMBED
  TL[embed] = neumann (qw/lambda1);
  TG[embed] = neumann (qw/lambda2);
#else
  for (int bid = 0; bid < nboundary; bid++) {
    TL[bid] = neumann (qw/lambda1);
    TG[bid] = neumann (qw/lambda2);
  }
#endif

  /**
  We store the initial mass of the ullage, which must be conserved. */

  mg0 = 0.;
  foreach (reduction(+:mg0))
    mg0 += tp2.rhov (&tsg)*(1. - f[])*dv();
}

/**
### Pressurization

The thermodynamic pressure is integrated in time using the pressurization rate
computed by the projection step. */

event pressurization (i++) {
  P0 += dt*dP0dt;
}

/**
## Post-Processing

The following lines of code are for post-processing purposes.
*/

/**
### Output Files

We write the thermodynamic pressure and the pressurization rate, together with
the corresponding analytic values, the relative variation of the mass of the
ullage, the average temperature of the ullage, and the maximum velocity. The
temperature of the gas phase is stored as a tracer, i.e. multiplied by the
volume fraction of the ullage. */

event output_data (i++) {
  char name[80];
  sprintf (name, "OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, "w");

#if EMBED
  double h = (level0 - 0.5)*L0, theta = acos (h/R0);
  double dP0dtexact = (gamma0 - 1.)*qw*2.*R0*theta /
    (sq(R0)*theta - h*sqrt (sq(R0) - sq(h)));
#else
  double Hg = (1. - level0)*L0;
  double dP0dtexact = (gamma0 - 1.)*qw*(L0 + 2.*Hg)/(L0*Hg);
#endif
  P0exact = 101325. + dP0dtexact*t;

  scalar TG = gas->T;
  double mg = 0., Tavg = 0., Vg = 0., umax = 0.;
  foreach (reduction(+:mg) reduction(+:Tavg) reduction(+:Vg)
      reduction(max:umax)) {
    double fg = 1. - f[];
    if (fg > F_ERR) {
      double Tg = TG[]/fg;
      mg += P0*MWG/(R_GAS*1.e3*Tg)*fg*dv();
      Tavg += Tg*fg*dv();
      Vg += fg*dv();
    }
    umax = max (umax, norm (u));
  }
  Tavg /= Vg;

  double relerr = fabs (P0 - P0exact)/P0exact;

  fprintf (fp, "%g %g %g %g %g %g %g %g %g\n", t, P0, P0exact, relerr,
      dP0dt, dP0dtexact, (mg - mg0)/mg0, Tavg, umax),
    fflush (fp);
}

/**
### Logger

We output the thermodynamic pressure of the tank in time (for testing). */

event logger (t += 0.1) {
  fprintf (stderr, "%d %.2f %.6g\n", maxlevel, t, P0);
}

/**
### Movie

We write the animation with the evolution of the temperature field. */

event movie (t += 0.01; t <= TEND) {
  if (maxlevel == 6) {
#if EMBED
    clear();
    view (tx = -0.5, ty = -0.5);
    draw_vof ("f", lw = 2.);
    draw_vof ("cs", "fs", filled = -1, fc = {1.,1.,1.});
    draw_vof ("cs", "fs", lw = 4.);
    squares ("T", min = TG0, max = TG0 + 100.);
    save ("movie.mp4");
#else
    clear();
    view (tx = -0.5, ty = -0.5);
    draw_vof ("f", lw = 2.);
    squares ("T", min = TG0, max = TG0 + 100.);
    box (notics = true, lw = 4.);
    save ("movie.mp4");
#endif
  }
}

/**
## Results

In all the following figures, the left panel refers to the square box, while the
right panel refers to the circular tank obtained with the embedded boundary.

The thermodynamic pressure increases linearly in time, following the analytic
solution in both tanks. The heat supplied by the walls is first stored in a
thin thermal boundary layer, and it is then distributed over the whole ullage
by conduction, but the pressurization rate does not depend on the temperature
distribution: after a short transient of about 1 ms, due to the compression
work which is computed using the pressurization rate of the previous time
step, it remains constant.

~~~gnuplot Evolution of the thermodynamic pressure
reset
set term @SVG size 900,450
set multiplot layout 1,2
set grid
set key top left
set xlabel "t [s]"
set ylabel "P_0 [Pa]"
set size square

set title "Square box"
plot "OutputData-5" u 1:3 every 50 w p ps 0.9 pt 6 lc rgb "black" t "Analytic", \
     "OutputData-4" u 1:2 w l lw 2 lc 1 t "LEVEL 4", \
     "OutputData-5" u 1:2 w l lw 2 lc 2 t "LEVEL 5", \
     "OutputData-6" u 1:2 w l lw 2 lc 3 t "LEVEL 6", \
     "OutputData-7" u 1:2 w l lw 2 lc 4 t "LEVEL 7", \
     "OutputData-8" u 1:2 w l lw 2 dt 2 lc 5 t "LEVEL 8"

set title "Embedded boundary"
plot "../heatflux-embed/OutputData-5" u 1:3 every 50 w p ps 0.9 pt 6 \
       lc rgb "black" t "Analytic", \
     "../heatflux-embed/OutputData-4" u 1:2 w l lw 2 lc 1 t "LEVEL 4", \
     "../heatflux-embed/OutputData-5" u 1:2 w l lw 2 lc 2 t "LEVEL 5", \
     "../heatflux-embed/OutputData-6" u 1:2 w l lw 2 lc 3 t "LEVEL 6", \
     "../heatflux-embed/OutputData-7" u 1:2 w l lw 2 lc 4 t "LEVEL 7", \
     "../heatflux-embed/OutputData-8" u 1:2 w l lw 2 dt 2 lc 5 t "LEVEL 8"
unset multiplot
~~~

In the square box, the relative error of the pressurization rate at
$t = 1$ s is $2 \cdot 10^{-5}$ on level 4, and about $4 \cdot 10^{-5}$ on
levels 5 to 8. It does not decrease with the grid resolution, and it grows
slowly in time: it is not a discretization error. The discrete pressurization
rate is equal to $(\gamma - 1)Q_g/V_g$, evaluated with the heat which enters
the ullage through the boundaries of the domain, up to $10^{-8}$, and the
error comes from the analytic solution, which assumes a flat interface. The gas
expands next to the walls and it is compressed in the core of the ullage, and
the pressure of this flow deforms the interface, which is not kept flat by
gravity or by surface tension. At $t = 1$ s the interface is about
$0.4~\mu$m lower next to the side walls and $0.2~\mu$m higher at the center
of the box, on all grids, while the volume of the ullage does not change. The
length of the walls wetted by the gas increases by $0.8~\mu$m, which is
$4 \cdot 10^{-5}$ of $L_0 + 2H_g$, and so does the heat supplied to the
ullage. The deformation is the same with no-slip walls, while it decreases if
the viscosity of the liquid, which opposes the displacement of the interface,
is increased.

In the circular tank, the relative error is $4.0 \cdot 10^{-3}$,
$1.1 \cdot 10^{-3}$, $3.2 \cdot 10^{-4}$, $1.2 \cdot 10^{-4}$ and
$5.1 \cdot 10^{-5}$ on levels 4 to 8. The pressurization rate is equal to
$(\gamma - 1)Q_g/V_g$ evaluated with the discrete heat input of the gas phase,
which is split exactly between the phases on the discrete geometry by
`phase_embed_flux()`. The error is due to the polygonal approximation of the
tank given by the embedded boundary: the length of the wall wetted by the gas
and the volume of the ullage are smaller than the analytic ones by
$2.1 \cdot 10^{-3}$ and $6.0 \cdot 10^{-3}$ on level 4, and both converge with
second order. On the finest grids the error approaches the one of the square
box, since the interface is deformed in the same way.

With the default flux of the embedded boundaries, `embed_flux()`, the heat flux
of the cut cells crossed by the interface is split using the average face
fraction of each phase, and the error is 1.7%, 1.1%, 0.71%, 0.31% and 0.085% on
levels 4 to 8. This error is of $O(\Delta)$, and its magnitude and sign depend
on the position of the interface within the cells at the contact line. For
example, at the beginning of the simulation the interface lies on a grid line,
and the gas cell next to the contact line receives only part of the heat of its
fragment of wall, because its interfacial face does not contain gas.

~~~gnuplot Pressurization rate
reset
set term @SVG size 900,450
set multiplot layout 1,2
set grid
set key bottom right
set xlabel "t [s]"
set ylabel "dP_0/dt [Pa/s]"
set size square

set title "Square box"
set yrange [15500:16500]
plot "OutputData-5" u 1:6 w l lw 2 lc rgb "black" t "Analytic", \
     "OutputData-4" u 1:5 w l lw 2 lc 1 t "LEVEL 4", \
     "OutputData-5" u 1:5 w l lw 2 lc 2 t "LEVEL 5", \
     "OutputData-6" u 1:5 w l lw 2 lc 3 t "LEVEL 6", \
     "OutputData-7" u 1:5 w l lw 2 lc 4 t "LEVEL 7", \
     "OutputData-8" u 1:5 w l lw 2 dt 2 lc 5 t "LEVEL 8"

set title "Embedded boundary"
set yrange [17000:18500]
plot "../heatflux-embed/OutputData-5" u 1:6 w l lw 2 lc rgb "black" \
       t "Analytic", \
     "../heatflux-embed/OutputData-4" u 1:5 w l lw 2 lc 1 t "LEVEL 4", \
     "../heatflux-embed/OutputData-5" u 1:5 w l lw 2 lc 2 t "LEVEL 5", \
     "../heatflux-embed/OutputData-6" u 1:5 w l lw 2 lc 3 t "LEVEL 6", \
     "../heatflux-embed/OutputData-7" u 1:5 w l lw 2 lc 4 t "LEVEL 7", \
     "../heatflux-embed/OutputData-8" u 1:5 w l lw 2 dt 2 lc 5 t "LEVEL 8"
unset multiplot
~~~

The interpolated order of convergence of the circular tank is shown in the
following figure, while the error of the square box is not interpolated, since
it does not depend on the grid. The convergence of the circular tank is of
second order on the coarse grids, and it slows down on the finest grids, where
the error approaches the one of the square box.

~~~gnuplot Convergence of the pressurization rate
reset
stats "<tail -n 1 OutputData-4" u (abs ($5 - $6)/$6) nooutput name "BOX4"
stats "<tail -n 1 OutputData-5" u (abs ($5 - $6)/$6) nooutput name "BOX5"
stats "<tail -n 1 OutputData-6" u (abs ($5 - $6)/$6) nooutput name "BOX6"
stats "<tail -n 1 OutputData-7" u (abs ($5 - $6)/$6) nooutput name "BOX7"
stats "<tail -n 1 OutputData-8" u (abs ($5 - $6)/$6) nooutput name "BOX8"
stats "<tail -n 1 ../heatflux-embed/OutputData-4" u (abs ($5 - $6)/$6) \
  nooutput name "EMBED4"
stats "<tail -n 1 ../heatflux-embed/OutputData-5" u (abs ($5 - $6)/$6) \
  nooutput name "EMBED5"
stats "<tail -n 1 ../heatflux-embed/OutputData-6" u (abs ($5 - $6)/$6) \
  nooutput name "EMBED6"
stats "<tail -n 1 ../heatflux-embed/OutputData-7" u (abs ($5 - $6)/$6) \
  nooutput name "EMBED7"
stats "<tail -n 1 ../heatflux-embed/OutputData-8" u (abs ($5 - $6)/$6) \
  nooutput name "EMBED8"

set print "errors"

print sprintf ("%d %.12f %.12f", 2**4, BOX4_mean, EMBED4_mean)
print sprintf ("%d %.12f %.12f", 2**5, BOX5_mean, EMBED5_mean)
print sprintf ("%d %.12f %.12f", 2**6, BOX6_mean, EMBED6_mean)
print sprintf ("%d %.12f %.12f", 2**7, BOX7_mean, EMBED7_mean)
print sprintf ("%d %.12f %.12f", 2**8, BOX8_mean, EMBED8_mean)

unset print

reset
set xlabel "Resolution"
set ylabel "Relative Error"

set xr[2**3:2**9]
set yr[1e-6:1e-1]
set key bottom left
set size square
set grid

set logscale x 2
set logscale y

f(x) = a*x**-b
fit f(x) "errors" u 1:3 via a,b

ftitle(a,b) = sprintf("%.3f/x^{%4.2f}", a, b)

plot "errors" u 1:2 pt 6 lc 1 title "Square box", \
     "errors" u 1:3 pt 8 lc 2 title "Embedded boundary", \
     f(x) w l lc 2 title ftitle(a, b)
~~~

The mass of the ullage is conserved: although the gas density changes locally
with the temperature, the variation of the total mass decreases slowly in time,
and it remains smaller than $1.3 \cdot 10^{-6}$ in the square box and
$6 \cdot 10^{-6}$ in the circular tank.

~~~gnuplot Relative variation of the mass of the ullage
reset
set term @SVG size 900,450
set multiplot layout 1,2
set grid
set key top right
set xlabel "t [s]"
set ylabel "(m_g - m_g^0)/m_g^0"
set size square

set title "Square box"
plot "OutputData-4" u 1:7 w l lw 2 lc 1 t "LEVEL 4", \
     "OutputData-5" u 1:7 w l lw 2 lc 2 t "LEVEL 5", \
     "OutputData-6" u 1:7 w l lw 2 lc 3 t "LEVEL 6", \
     "OutputData-7" u 1:7 w l lw 2 lc 4 t "LEVEL 7", \
     "OutputData-8" u 1:7 w l lw 2 dt 2 lc 5 t "LEVEL 8"

set title "Embedded boundary"
plot "../heatflux-embed/OutputData-4" u 1:7 w l lw 2 lc 1 t "LEVEL 4", \
     "../heatflux-embed/OutputData-5" u 1:7 w l lw 2 lc 2 t "LEVEL 5", \
     "../heatflux-embed/OutputData-6" u 1:7 w l lw 2 lc 3 t "LEVEL 6", \
     "../heatflux-embed/OutputData-7" u 1:7 w l lw 2 lc 4 t "LEVEL 7", \
     "../heatflux-embed/OutputData-8" u 1:7 w l lw 2 dt 2 lc 5 t "LEVEL 8"
unset multiplot
~~~

Differently from [pressurization.c](pressurization.c), the velocity is not
zero: the gas expands in the thermal boundary layer next to the walls, where
the divergence of the conductive flux is positive, and it is compressed in the
core of the ullage, where the velocity divergence is
$-a = -\left(dP_0/dt\right)/\left(\gamma P_0\right)$. The gas therefore moves
away from the heated walls, towards the core of the ullage. In the circular
tank, mirroring the ullage with respect to the interface gives an axisymmetric
problem, whose solution is a radial velocity field $u_r(r)$, tangent to the
interface. In the core, as long as the thermal boundary layer is thin,
$u_r = -a r/2$, while $u_r = 0$ on the wall. On level 6, the volume average of
the radial component
$\langle \mathbf{u}\cdot\hat{\mathbf{r}}\rangle/\langle|\mathbf{u}|\rangle$
is between $-0.8$ and $-0.94$ for $0.005 \le t \le 0.1$ s, and at $t = 0.05$ s
the radial velocity at $r = R/2$ is $-1.36 \cdot 10^{-4}$ m/s, while $-aR/4 =
-1.41 \cdot 10^{-4}$ m/s. The liquid, which is incompressible, remains almost at
rest ($|\mathbf{u}| < 3 \cdot 10^{-5}$ m/s).

In both tanks, the velocity is largest at the beginning of the simulation, when
the thermal boundary layer is thin, and it decreases exponentially as the
temperature profile approaches a quasi-steady shape. In the circular tank, the
relaxation time is $R^2/(\alpha j_1'^2) \approx 0.06$ s, where $\alpha$ is the
thermal diffusivity of the gas and $j_1' \approx 3.83$ is the first zero of the
derivative of the Bessel function $J_0$; the ullage of the square box is
larger, and the relaxation is slower. The maximum velocity then has a minimum,
at $t \approx 0.48$ s in the square box and at $t \approx 0.29$ s in the
circular tank, where the direction of the flow is reversed, and it reaches a
plateau of about $4 \cdot 10^{-6}$ and $2.2 \cdot 10^{-6}$ m/s respectively,
which is the same on all grids. For an ideal gas
$\nabla\cdot\mathbf{u} = -D\ln\rho/Dt = DT/Dt/T - d\ln P_0/dt$. In the
quasi-steady state the temperature of all the fluid elements increases at the
same rate, therefore the hot gas next to the walls expands less than the cold
gas in the core, while the compression is uniform: the divergence becomes
slightly negative next to the walls, and the gas moves towards the walls. In
the circular tank, the magnitude of this flow is of the order of
$a (\Delta T/T) R/8 \approx 2 \cdot 10^{-6}$ m/s, where
$\Delta T = q_w R/(2\lambda_g) \approx 9$ K is the quasi-steady temperature
difference between the wall and the center of the tank.

~~~gnuplot Maximum velocity
reset
set term @SVG size 900,450
set multiplot layout 1,2
set grid
set key bottom right
set xlabel "t [s]"
set ylabel "max |u| [m/s]"
set logscale y
set format y "10^{%T}"
set size square

set title "Square box"
plot "OutputData-4" u 1:9 w l lw 2 lc 1 t "LEVEL 4", \
     "OutputData-5" u 1:9 w l lw 2 lc 2 t "LEVEL 5", \
     "OutputData-6" u 1:9 w l lw 2 lc 3 t "LEVEL 6", \
     "OutputData-7" u 1:9 w l lw 2 lc 4 t "LEVEL 7", \
     "OutputData-8" u 1:9 w l lw 2 dt 2 lc 5 t "LEVEL 8"

set title "Embedded boundary"
plot "../heatflux-embed/OutputData-4" u 1:9 w l lw 2 lc 1 t "LEVEL 4", \
     "../heatflux-embed/OutputData-5" u 1:9 w l lw 2 lc 2 t "LEVEL 5", \
     "../heatflux-embed/OutputData-6" u 1:9 w l lw 2 lc 3 t "LEVEL 6", \
     "../heatflux-embed/OutputData-7" u 1:9 w l lw 2 lc 4 t "LEVEL 7", \
     "../heatflux-embed/OutputData-8" u 1:9 w l lw 2 dt 2 lc 5 t "LEVEL 8"
unset multiplot
~~~
*/
