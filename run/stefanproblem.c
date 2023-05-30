/**
# Stefan Problem

A thin vapor layer is initialized close to a superheated
wall. The high wall temperature heats up the vapor layer,
which is in contact with a liquid phase that always remains
at saturation temperature leading to the phase change. The
simulation setup used here was adapted from [Malan et al, 2021](#malan2021geometric).
The analytic solution to the problem describes the evolution
of the vapor layer thickness:
$$
  \delta(t) = 2\lambda\sqrt{\alpha_g t}
$$
where $\alpha_g$ is the thermal diffusivity, defined as
$\alpha_g = \lambda_g/\rho_g/cp_g$, while $\lambda$ is
a value, specific to the physical properties under
investigation, which can be found from the solution
of the transcendental equation:
$$
  \lambda\exp(\lambda^2) \text{erf}(\lambda) =
  \dfrac{cp_g(T_{wall}-T_{sat})}{\Delta h_{ev}\sqrt{\pi}}
$$
The animation shows that, at the beginning of the simulation,
the superheated wall heats up the vapor layer which becomes
hotter than the liquid phase, leading to the evaporation
which establishes the interface velocity jump.
![Evolution of the interface and the temperature field](stefanproblem/movie.mp4)(width="400" height="400")
*/

/**
## Phase Change Setup

We move the interface using the velocity *uf*, with the
expansion term shifted toward the gas-phase. In this way
*uf* is divergence-free at the interface. The double
pressure velocity couping is used to obtain an extended
velocity, used to transport the gas phase tracers. */

#define SHIFT_TO_GAS
#define INT_USE_UF
#define CONSISTENTPHASE2

/**
## Simulation Setup

We use the centered solver with the divergence source term,
and the extended velocity is obtained using the centered-doubled
procedure. The evaporation model is used in combination
with the temperature-gradient mechanism, which manages
the solution of the temperature field.
*/

#include "grid/multigrid.h"
#include "navier-stokes/centered-evaporation.h"
#include "navier-stokes/centered-doubled.h"
#include "two-phase.h"
#include "tension.h"
#include "evaporation.h"
#include "temperature-gradient.h"
#include "view.h"

/**
We declare the variables required by the
temperature-gradient model. */

double lambda1, lambda2, cp1, cp2, dhev;
double TL0, TG0, TIntVal, Tsat, Twall;

/**
### Boundary conditions

Outflow boundary conditions for velocity and pressure are
imposed on the left wall, while the temperature on the
superheated wall is set to the superheated value.
*/

u.n[left] = neumann (0);
u.t[left] = neumann (0);
p[left] = dirichlet (0);
uext.n[left] = neumann (0);
uext.t[left] = neumann (0);
pext[left] = dirichlet (0);

TG[left]  = dirichlet (Tsat);
TG[right] = dirichlet (Twall);
TL[left]  = dirichlet (Tsat);
TL[right] = dirichlet (Twall);

/**
### Problem Data

We declare the maximum and minimum levels of refinement,
the $\lambda$ parameter, and the initial thickness of
the vapor layer.
*/

int maxlevel, minlevel = 3;
double lambdaval = 0.06779249298045148;
double  delta0 = 322.5e-6;
double tshift, teff;


int main (void) {
  /**
  We set the material properties of the
  fluids. */

  rho1 = 958., rho2 = 0.6;
  mu1 = 2.82e-4, mu2 = 1.23e-5;
  lambda1 = 0.68, lambda2 = 0.025;
  cp1 = 4216., cp2 = 2080.;
  dhev = 2.256e6,

  /**
  The initial temperature and the interface
  temperature are set to the same value. */

  Tsat = 373., Twall = 383.;
  TL0 = Tsat, TG0 = Tsat; TIntVal = Tsat;

  /**
  We change the dimension of the domain
  and the surface tension coefficient. */

  L0 = 10e-3;
  f.sigma = 0.059;

  /**
  We define a list with the maximum time
  steps and the maximum levels of refinement. */

  double dtlist[] = {0.005, 0.001, 0.0005, 0.00005};
  int mllist[] = {5, 6, 7, 8};

  /**
  We run the simulation for different maximum
  levels of refinement. */

  for (int sim = 0; sim < 1; sim++) {
    maxlevel = mllist[sim];
    DT = dtlist[sim];
    init_grid (1 << maxlevel);
    run();
  }
}

/**
We initialize the volume fraction field and the temperature
in the gas and in liquid phase. */

event init (i = 0) {
  fraction (f, -(x - 0.0096775));

  foreach() {
    TL[] = f[]*TL0;
    TG[] = (1. - f[])*TG0;
    T[] = f[]*TL0 + (1. - f[])*TG0;
  }
  boundary({T});
  boundary({TL,TG});

  /**
  At simulation time equal to zero, the thickness of the
  vapor layer is not zero. Therefore, we compute a time
  *shifting factor* (just for post-processing purposes). */

  double effective_height;
  effective_height = (sq(L0) - statsf(f).sum)/L0;
  tshift = sq(effective_height/2./lambdaval)*rho2*cp2/lambda2;
}

/**
We refine the interface and the region where the
temperature field changes. */

#if TREE
event adapt (i++) {
  adapt_wavelet_leave_interface ({T}, {f},
      (double[]){1.e-3,1.e-3}, maxlevel, minlevel, 1);
}
#endif

/**
## Post-Processing

The following lines of code are for post-processing purposes.
*/

/**
### Exact Solution

We write a function that computes the exact
solution to the thickness of the vapor layer.
*/

double exact (double time) {
  return 2.*lambdaval*sqrt(lambda2/rho2/cp2*time);
}

/**
### Output Files

We write the thickness of the vapor layer and the analytic
solution on a file.
*/

event output (t += 0.1) {
  double effective_height = 0.;
  foreach(reduction(+:effective_height))
    effective_height += (1. - f[])*dv();
  effective_height /= L0;

  double relerr = fabs (exact(t+tshift) - effective_height) / exact(t+tshift);

  char name[80];
  sprintf (name, "OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, "w");

  fprintf (fp, "%g %g %g %g\n", t+tshift, effective_height, exact (t+tshift), relerr);
  fflush (fp);
}

/**
### Movie

We write the animation with the evolution of the
temperature field and the gas-liquid interface.
*/

event movie (t += 0.1; t <= 10) {
  clear();
  view (tx = -0.5, ty = -0.5);
  box();
  draw_vof ("f", lw = 1.5);
  squares ("T", min = Tsat, max = Twall, linear = true);
  vectors ("u", scale = 1.e-1, lc = {1.,1.,1.});
  save ("movie.mp4");
}

/**
## Results

~~~gnuplot Evolution of the vapor layer thickness
set xlabel "t[s]"
set ylabel "Vapor Layer Thickness [m]"
set key left top
set size square
set grid

plot "OutputData-5" every 10 u 1:3 w p ps 2 t "Analytic", \
     "OutputData-5" u 1:2 w l lw 2 t "LEVEL 5"
~~~

## References

~~~bib
@article{malan2021geometric,
  title={A geometric VOF method for interface resolved phase change and conservative thermal energy advection},
  author={Malan, LC and Malan, Arnaud G and Zaleski, St{\'e}phane and Rousseau, PG},
  journal={Journal of Computational Physics},
  volume={426},
  pages={109920},
  year={2021},
  publisher={Elsevier}
}
~~~
*/
