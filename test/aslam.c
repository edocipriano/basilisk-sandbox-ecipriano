/**
# Aslam Extrapolations

We try to use the constant and linear Aslam extrapolations
([Aslam 2003](#aslam2004partial)), defined in [aslam.h](../src/aslam.h).
The problem is characterized by a squared domain, with
dimensions $(-\pi,\pi)\times(-\pi,\pi)$, and by a level
set function $\phi = \sqrt{x^2 + y^2} - 2$. The field *u*
exist just in a region of the domain, and must be extrapolated:

$$
u =
\begin{cases}
  0 & \text{if } \phi > 0,\\
  \cos(x)\sin(x) & \text{if } \phi \leq 0.
\end{cases}
$$
*/

#include "utils.h"
#include "aslam.h"
#include "view.h"

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

/**
We declare the level set field *levelset*, and
the field to extrapolate *u*. */

scalar levelset[], u[];

int main (void) {

  /**
  We set the domain geometry and we initialize
  the grid. */

  size (2.*pi);
  origin (-pi,-pi);
  init_grid (1 << 8);
  double R0 = 2.;

  /**
  We initialize the level set function. */

  foreach()
    levelset[] = sqrt (sq (x) + sq (y)) - R0;

  /**
  We initialize the function *u* to be extrapolated
  and we call the *constant_extrapolation()* function.
  The extrapolations are perfomed using a $\Delta t$
  of 0.01 (it must be small enough to guarantee the
  stability of the explicit in time discretization),
  and using a total number of time steps equal to 300,
  in order to obtain a steady-state solution. */

  foreach()
    u[] = (levelset[] <= 0.) ? cos(x)*sin(y) : 0.;
  write_picture ("initial.png");

  constant_extrapolation (u, levelset, dt=0.01, n=300);
  write_picture ("constant.png");

  /**
  We re-initialize the function *u* and we apply the
  *linear_extrapolation()*. */

  foreach()
    u[] = (levelset[] <= 0.) ? cos(x)*sin(y) : 0.;
  linear_extrapolation (u, levelset, dt=0.01, n=300);
  write_picture ("linear.png");
}

/**
## Results

![Initial field](aslam/initial.png)(width="800" height="600")
![Constant extrapolation](aslam/constant.png)(width="800" height="600")
![Linear extrapolation](aslam/linear.png)(width="800" height="600")

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

