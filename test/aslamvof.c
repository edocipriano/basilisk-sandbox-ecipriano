/**
# Aslam Extrapolations

This test case is similar to [aslam.c](aslam.c). However,
instead of starting from the level set field, we start from
a vof field, which is converted into level set, and then
Aslam extrapolations are performed.
*/

#include "utils.h"
#include "aslam.h"
#include "view.h"
#include "BOYD/src/LS_funcs/LS_reinit.h"

/**
We define a function that writes a picture with the
map and isolines of the extended scalar fields.
*/

void write_picture (char* name) {
  clear();
  isoline ("levelset", val = 0., lw = 2.);
  isoline("u", n = 30);
  squares ("u", spread = -1);
  box();
  save (name);
}

void write_levelset (void) {
  clear();
  isoline ("levelset", val = 0., lw = 2.);
  squares ("levelset", spread = -1);
  box();
  save ("levelset.png");
}

/**
We define a function that converts the vof fraction
to the level set field. */

void vof2ls (scalar f, scalar levelset) {
  double deltamin = L0/(1 << grid->maxdepth);
  foreach()
    levelset[] = -(2.*f[] - 1.)*deltamin*0.75;
#if TREE
  restriction({levelset});
#endif
  LS_reinit (levelset, dt = 0.5*L0/(1 << grid->maxdepth),
      it_max = 0.5*(1 << grid->maxdepth));
}

#define ufunc(x,y)(x*y)
# define circle(x,y,R)(sq(R) - sq(x) - sq(y))

/**
We declare the level set field *levelset*, and
the field to extrapolate *u*. */

scalar levelset[], u[], f[];

int main (void) {

  /**
  We set the domain geometry and we initialize
  the grid. */

  size (2.*pi);
  origin (-pi,-pi);
  init_grid (1 << 8);
  double R0 = 2.;

  /**
  We initialize the volume fraction field. */

  fraction (f, (sq(R0) - sq(x) - sq(y)));

  /**
  We reconstruct the levelset field. */

  vof2ls (f, levelset);
  write_levelset();

  /**
  We initialize the function *u* to be extrapolated
  and we call the *constant_extrapolation()* function.
  The extrapolations are perfomed using a $\Delta t$
  of 0.01 (it must be small enough to guarantee the
  stability of the explicit in time discretization),
  and using a total number of time steps equal to 300,
  in order to obtain a steady-state solution. */

  foreach()
    u[] = ufunc(x,y)*f[];
  write_picture ("initial.png");

  double dtmin = 0.5*L0/(1 << grid->maxdepth);

  constant_extrapolation (u, levelset, dt=dtmin, n=300, c=f);
  write_picture ("constant.png");
  fprintf (stderr, "constant = %g\n", statsf(u).sum);

  /**
  We re-initialize the function *u* and we apply the
  *linear_extrapolation()*. */

  foreach()
    u[] = ufunc(x,y)*f[];
  linear_extrapolation (u, levelset, dt=dtmin, n=300, c=f);
  write_picture ("linear.png");
  fprintf (stderr, "linear = %g\n", statsf(u).sum);
}

/**
## Results

![Level Set](aslam/levelset.png)

![Initial field](aslam/initial.png)

![Constant extrapolation](aslam/constant.png)

![Linear extrapolation](aslam/linear.png)

## References

~~~bib
@article{aslam2004partial,
  title={A partial differential equation approach to multidimensional extrapolation},
  author={Aslam, Tariq D},
  journal={Journal of Computational Physics},
  volume={193},
  number={1},
  pages={349--355},
  year={2004},
  publisher={Elsevier}
}
~~~
*/
