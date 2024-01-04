/**
# Test cases for the Poisson--Helmholtz solver with Ghost Fluid Method

We reproduce the test cases proposed by [Liu et al. 2000](#liu2000boundary) to
verify the implementation of the GFM approach to impose interface boundary
conditions in the solution of Poisson--Helmholtz equations.
*/

#include "grid/multigrid1D.h"
#include "poisson-gfm.h"
#include "fractions.h"
#include "run.h"

int maxlevel = 11;

scalar a[], f[];

int main (void) {
  init_grid (1 << maxlevel);
  run();
}

event init (i = 0) {
  L0 = 1 [0];
  fraction (f, x - 0.5*L0);

  /**
  We define interfacial faces according to the following condition:

  $$
    \left(f_i - 0.5\right)\left(f_{i-1} - 0.5\right) < 0
  $$
  */

  face vector intface[];
  foreach_face()
    intface.x[] = ((f[] - 0.5)*(f[-1] - 0.5) < 0.) ? 1. : 0.;

  /**
  We define the relative position of the interface on the current face as:

  $$
    \lambda = \dfrac{f_{i-1} - 0.5}{f_{i-1} - f_{i}}
  $$
  */

  face vector reldist[];
  foreach_face()
    reldist.x[] = (intface.x[] == 1.) ? (f[-1] - 0.5)/(f[-1] - f[]) : 0.;

  /**
  ## Case 1: Laplace equation

  We solve the following Laplace equation:
  $$
    \nabla^2 a = 0
  $$

  with:
  $$
    \left[a\right] = 1
  $$
  $$
  \left[\nabla a\cdot\mathbf{n}_\Gamma\right]_\Gamma = 0
  $$
  */

  a[left] = dirichlet (0.);
  a[right] = dirichlet (2.);

  scalar b[];
  foreach()
    b[] = 0.;

  face vector ajump[];
  foreach_face()
      ajump.x[] = 1.;

  poisson_ghost (a, b, f = f, ajump = ajump,
      tolerance = TOLERANCE/sq(dt), nrelax = 4);

  for (double x = 0.; x < L0; x += 0.5*L0/(1 << maxlevel))
    fprintf (stderr, "%g %g %g\n", x, interpolate (f, x), interpolate (a, x));
  fprintf (stderr, "\n\n");

  /**
  ~~~gnuplot
  reset
  set xr[0:1]
  set grid
  set xlabel "x"
  set ylabel "a"
  set key top left

  f(x) = (x <= 0.5) ? x : x + 1

  plot "log" index 0 u 1:3 w l lw 2 t "Numerical", \
       f(x) w p pt 34 t "Exact"
  ~~~
  */

  /**
  ## Case 2: Laplace equation

  We solve the following Laplace equation:
  $$
    \nabla^2 a = 0
  $$

  with:
  $$
    \left[a\right] = 0
  $$
  $$
  \left[\nabla a\cdot\mathbf{n}_\Gamma\right]_\Gamma = 1
  $$
  */

  a[left] = dirichlet (0.);
  a[right] = dirichlet (1.5);

  foreach() {
    a[] = 0.;
    b[] = 0.;
  }

  face vector bjump[];
  foreach_face()
    bjump.x[] = 1.;

  poisson_ghost (a, b, f = f, bjump = bjump,
      tolerance = TOLERANCE/sq(dt), nrelax = 4);

  for (double x = 0.; x < L0; x += 0.5*L0/(1 << maxlevel))
    fprintf (stderr, "%g %g %g\n", x, interpolate (f, x), interpolate (a, x));
  fprintf (stderr, "\n\n");

  /**
  ~~~gnuplot
  reset
  set xr[0:1]
  set grid
  set xlabel "x"
  set ylabel "a"
  set key top left

  f(x) = (x <= 0.5) ? x : 2.*x - 0.5

  plot "log" index 1 u 1:3 w l lw 2 t "Numerical", \
       f(x) w p pt 34 t "Exact"
  ~~~
  */

  /**
  ## Case 3: Poisson equation

  We solve the following Poisson equation:
  $$
    \nabla\cdot\left(\beta\nabla a\right) = (8x^2 - 4)e^{-x^2}
  $$

  with:
  $$
    \left[a\right]_\Gamma = \begin{cases}
      -e^{-0.09} & \text{if } x = 0.4 \\
      -e^{-0.36} & \text{if } x = 0.6
    \end{cases}
  $$
  $$
  \left[\beta\nabla a\cdot\mathbf{n}_\Gamma\right]_\Gamma = \begin{cases}
    -1.2e^{-0.09} & \text{if } x = 0.4 \\
     2.4e^{-0.36} & \text{if } x = 0.6
  \end{cases}
  $$
  */

  a[left] = dirichlet (0.);
  a[right] = dirichlet (0.);

  fraction (f, fabs (x - 0.45) - 0.15);

  intfacefun (f);
  intsidefun (f);
  reldistfun (f);

  foreach() {
    a[] = 0.;
    b[] = (x >= 0.3 && x <= 0.6) ? (8.*sq(x) - 4.)*exp(-sq(x)) : 0.;
  }

  foreach_face()
    ajump.x[] = (x <= 0.5) ? -exp(-0.09) : -exp(-0.36);

  double betaL = 1.;  // outer region
  double betaG = 2.;  // inner region

  face vector betac[];
  foreach_face() {
    if (intface.x[]) {
      if (intside.x[] > 0.)
        betac.x[] = betaG*betaL/(reldist.x[]*betaG + (1. - reldist.x[])*betaL);
      else
        betac.x[] = betaG*betaL/(reldist.x[]*betaG + (1. - reldist.x[])*betaL);
    }
    else {
      if (intside.x[] > 0.)
        betac.x[] = betaG;
      else
        betac.x[] = betaL;
    }
  }

  foreach_face() {
    double bjumpval = (x <= 0.5) ? -1.2*exp(-0.09) : 2.4*exp(-0.36);
    if (intface.x[]) {
      if (intside.x[] > 0.)
        bjump.x[] = bjumpval/betaL;
      else
        bjump.x[] = bjumpval/betaG;
    }
    else
      bjump.x[] = 0.;
  }

  poisson_ghost (a, b, f = f, ajump = ajump, bjump = bjump, alpha = betac,
      tolerance = TOLERANCE/sq(dt), nrelax = 4);

  for (double x = 0.; x < L0; x += 0.5*L0/(1 << maxlevel))
    fprintf (stderr, "%g %g %g\n", x, interpolate (f, x), interpolate (a, x));
  fprintf (stderr, "\n\n");

  /**
  ~~~gnuplot
  reset
  set xr[0:1]
  set grid
  set xlabel "x"
  set ylabel "a"
  set key top left

  f(x) = (x >= 0.3 && x <= 0.6) ? exp(-x**2) : 0.

  plot "log" index 2 u 1:3 w l lw 2 t "Numerical", \
       f(x) w p pt 34 t "Exact"
  ~~~
  */
}

/**
## References

~~~bib
@article{liu2000boundary,
  title={A boundary condition capturing method for Poisson's equation on irregular domains},
  author={Liu, Xu-Dong and Fedkiw, Ronald P and Kang, Myungjoo},
  journal={Journal of computational Physics},
  volume={160},
  number={1},
  pages={151--178},
  year={2000},
  publisher={Elsevier}
}
~~~
*/

