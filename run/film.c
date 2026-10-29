/**
# Unresolved Thin Film

A liquid film, thinner than one grid cell, lies on a superheated wall
($T_w$) and evaporates into a vapor phase at saturation temperature ($T_s$).
The test assesses the interface gradients computed when the film is not
resolved by the grid: the liquid side of the stencil reaches the wall
directly.

For small Jakob numbers the temperature profile in the film is linear. A thermal
resistance $R_c$ (interfacial, with accommodation coefficient $\omega$) is in
series with the film, and the heat flux at the interface is:
$$
  q = \dfrac{T_w - T_s}{R_c + h/\lambda_l}
$$
The liquid is at rest, and the film is consumed at $\dot{h} = -\dot{m}/\rho_l
= -q/(\rho_l \Delta h_{ev})$:
$$
  h(t) = \sqrt{(h_0 + \lambda_l R_c)^2 - 2\dfrac{\lambda_l (T_w -
  T_s)}{\rho_l \Delta h_{ev}} t} - \lambda_l R_c
$$
until the dry-out time $t_d = \rho_l \Delta h_{ev} h_0 (R_c +
h_0/(2\lambda_l))/(T_w - T_s) \approx 21$ s.

Three ways of computing the liquid heat flux $q$ (and $\dot{m} = -q/\Delta
h_{ev}$) are compared, one after the other:

1. *ghost*: default [plic_gradient()](/src/gradients.h), the probe reaches the
ghost cell beyond the wall (one point gradient along the probe);
2. *degenerate*: the wall face is closed to the probe ($f_s = 0$), the gradient
falls back to the degenerate case (cell and interface values);
3. *SGS*: sub-grid model with a linear profile in the film, $q = (T_w -
T_s)/(R_c + \delta/\lambda_l)$, with $\delta = f\Delta$. */

#include "navier-stokes/low-mach.h"
#include "two-phase.h"
#include "tension.h"
#include "boiling.h"

/**
### Boundary conditions

Superheated wall, wetted by the liquid, on the left. The resistance is imposed
with the Robin condition $T_B + \lambda_l R_c \partial_n T = T_w$, i.e.
[navier()](/src/grid/cartesian-common.h#navier). */

double Ts = 373., Tw = 383., Rc;
T[left] = navier (Tw, lambda1*Rc);
f[left] = dirichlet (1.);

/**
### Problem Data

The initial film thickness is half a cell. */

int MAXLEVEL = 2, CASE;
const double H0 = 0.5, TEND = 30.;

int main (void) {
  rho1 = 958., rho2 = 0.59;
  mu1 = 2.82e-4, mu2 = 1.23e-6;
  lambda1 = 0.6, lambda2 = 0.026;
  cp1 = 4216., cp2 = 2034.;
  dhev = 2.257e+6;

  TL0 = TG0 = TIntVal = Ts;

  double fac = 0.01;
  Rc = (2. - fac)/(2.*fac)*
    sqrt (2.*pi*8.314*pow (TIntVal, 3.)/0.018)/rho2/sq (dhev);

  fprintf (stderr, "params %g %g %g\n", Rc, lambda1, Tw - Ts);

  nv = 1;
  pcm.consistent = true;
  pcm.expansion = false;

  DT = 1e-3;
  size (2.6e-3);
  for (CASE = 0; CASE < 3; CASE++) {
    init_grid (1 << MAXLEVEL);
    run();
  }
}

/**
In cases 2 and 3 the evaporation rate is the field *mSGS* computed below,
which replaces the default *mEvap* in the total evaporation rate. */

scalar mSGS[], * mSGSList = {mSGS}, * mPLICList = {mEvap};

event init (i = 0) {
  double delta = L0/(1 << MAXLEVEL);
  fraction (f, H0*delta - x);

  scalar TL = liq->T, TG = gas->T;
  copy_bcs ({TL,TG}, T);

  mEvapList = CASE ? mSGSList : mPLICList;
}

/**
## Post-Processing

### Gradient stencils

The field *deg* stores how the liquid interface gradient is computed in
[dirichlet_gradient()](/src/gradients.h) in the interfacial cells (*nodata*
elsewhere): 0 degenerate (cell value and interface value), 1 one point in
the liquid, 2 two points in the liquid (VOF-biased or third order), 3 SGS
model.

This event runs before the one in
[temperature-gradient.h](/src/temperature-gradient.h), which uses the ghost
gradient $g_1$ for the liquid temperature source. In cases 2 and 3 we add the
difference with the gradient $g$ of the case, to be consistent with
$\dot{m}$. */

scalar deg[];

event phasechange (i++) {
  phase_tracers_to_scalars (liq, f, tol = F_ERR);

  scalar fl[];
  scalar_clone (fl, f);
  foreach()
    fl[] = f[];

  face vector fsl[], fsd[];
  face_fraction (fl, fsl);
  foreach_face()
    fsd.x[] = fsl.x[];

  foreach_boundary (left)
    fsd.x[] = 0.;

  scalar TL = liq->T, slT = liq->STexp;
  foreach() {
    deg[] = nodata;
    if (f[] > F_ERR && f[] < 1. - F_ERR) {
      dstats s;
      double g1 = plic_gradient (point, TL, fl, fsl, TIntVal, true, &s);
      double g = CASE == 0 ? g1 :
        CASE == 1 ? plic_gradient (point, TL, fl, fsd, TIntVal, true, &s) :
        (TIntVal - Tw)/(f[]*Delta + lambda1*Rc);
      deg[] = CASE == 2 ? 3. :
        (s.v[0] == nodata) ? 0. : (s.v[1] == nodata) ? 1. : 2.;
      mSGS[] = lambda1*g/dhev;
      slT[] += lambda1*(g - g1)*vof_reconstruction (point, fl).dirac;
    }
  }

  phase_scalars_to_tracers (liq, f);
}

/**
### Exact Solution */

double hexact (double t) {
  double h0 = H0*L0/(1 << MAXLEVEL), lR = lambda1*Rc;
  return max (0., sqrt (max (0., sq (h0 + lR)
          - 2.*lambda1*(Tw - Ts)/(rho1*dhev)*t)) - lR);
}

double qexact (double t) {
  double h = hexact (t);
  return h > 0. ? (Tw - Ts)/(Rc + h/lambda1) : 0.;
}

/**
### Output

We write on the log (tagged by the case) the film thickness (liquid volume over the wall length) and the
interface heat flux $q = -\Delta h_{ev}\dot{m}$ (the evaporation rate is
negative by convention), averaged over the interfacial cells. A line is
written only when $h$ changes by 1\% of $h_0$, or $q$ by 2\%, or every 0.25 s:
the output is dense only close to dry-out. The start-up transient ($t <
5$ ms) is skipped. */

event output (i++) {
  double h = statsf (f).sum/L0;
  double q = 0.; int n = 0;
  foreach (reduction(+:q) reduction(+:n))
    if (f[] > F_ERR && f[] < 1. - F_ERR)
      q -= dhev*mEvapTot[], n++;
  q /= max (n, 1);

  static double tl = HUGE, hl, ql;
  if (t < 0.005)
    return 0;
  if (t < tl || t - tl >= 0.25 || t >= TEND || fabs (q - ql) > 0.02*fabs (ql) ||
      fabs (h - hl) > 0.01*H0*L0/(1 << MAXLEVEL)) {
    fprintf (stderr, "case %d %.4g %.4g %.4g %.4g %.4g\n",
        CASE, t, h, hexact (t), q, qexact (t));
    tl = t, hl = h, ql = q;
  }
}

event logger (t += 1) {
  fprintf (stderr, "%d %.3f %.4g\n", i, t, statsf (f).sum);
}

/**
We stop at dry-out (or at *TEND*). */

event stop (i++; t <= TEND) {
  return statsf (f).sum <= 0.;
}

/**
## Results

~~~gnuplot Film thickness
cases = "Ghost Degenerate SGS"
set xlabel "t [s]"
set ylabel "h [m]"
set key right top
set size square
set grid
plot '<grep "case 0" log' every 2 u 3:5 w p lt 1 t "Exact", \
     for [k=0:2] sprintf ('<grep "case %d" log', k) u 3:4 w l lt k+2 t word(cases,k+1)
~~~

~~~gnuplot Interface heat flux
set ylabel "q [W/m^2]"
set key left top
plot '<grep "case 0" log' every 2 u 3:7 w p lt 1 t "Exact", \
     for [k=0:2] sprintf ('<grep "case %d" log', k) u 3:6 w l lt k+2 t word(cases,k+1)
~~~

The heat flux as a function of the film thickness shows the two regimes: $q
\approx \lambda_l \Delta T/h$ (film-limited) for $h \gg \lambda_l R_c$, and the
saturation $q \to \Delta T/R_c$ (resistance-limited) for $h \ll \lambda_l
R_c$. They cross at $h = \lambda_l R_c$ (vertical line), where the
resistance and the film have the same length scale.

~~~gnuplot Heat flux vs film thickness
stats '<grep "^params" log' u (Rc = $2, lam = $3, dT = $4) nooutput
unset logscale
set logscale xy
set xrange [1e-6:4e-4]
set xlabel "h [m]"
set ylabel "q [W/m^2]"
set key outside
set samples 25
lR = lam*Rc
set arrow from lR, graph 0 to lR, graph 1 nohead dt 2 lc rgb "black"
set label "h = {/Symbol l}_l R_c" at lR, graph 0.04 right offset -0.5,0
set label "{/Symbol D}T/R_c" at graph 0.97, first dT/Rc right offset 0,0.8
plot dT/(Rc + x/lam) w p lt 1 t "Exact", \
     dT/(x/lam) w p pt 2 lc 1 t "Exact, R_c = 0", \
     for [k=0:2] sprintf ('<grep "case %d" log', k) u 4:6 w l lt k+2 t word(cases,k+1), \
     dT/Rc w l dt 2 lc rgb "black" notitle
~~~
*/
