/**
# Centripetal Force

Implementation of the method proposed by [Saufi et al. 2019](#saufi2019dropletsmoke)
that allows to suspend liquid droplets on a specific point of the domain.
The idea is to to use an artificial suspending force to maintain
the center of the droplet approximately at the same coordinate,
also when the gravity force is present.

This method defines a centripetal potential as:
$$
\xi = \xi_0 \dfrac{1}{\mid {\mathbf{x} - \mathbf{x}_c}\mid}
$$
The suspending force is computed from the centripetal potential and
added to the acceleration terms in the momentum equation:
$$
\mathbf{f}_m = \rho \alpha \nabla \xi
$$
*/

#define CENTRIPETAL

/**
## Fields

We initialize the centripetal potential field *phicentripetal*.
*/

scalar phicentripetal[];

/**
## User Data

We need to define the coordinate of the point where the center of
the droplet should remain *p*, and a parameter that controls the
intensity of the suspending force *eps*.
*/

struct SuspendingForceModel {
  coord p;
  double eps;
};

struct SuspendingForceModel sfm = {
  .p = {0., 0., 0.},
  .eps = 1.25e-4,
};

/**
## Implementation

The droplet should be initialized at the position where it should
remain fixed. Therefore, the centripetal potential is computed
at the first time step.
*/

event init (i = 0) {
  foreach() {
    double magdistance = 0.;
    magdistance = sq(x - sfm.p.x) + sq(y - sfm.p.y);
#if dimension > 2
    magdistance += sq(z - sfm.p.z);
#endif
    magdistance = sqrt (magdistance);
    phicentripetal[] = 1./magdistance;
  }
  boundary({phicentripetal});
}

/**
We overload the acceleration event computing the suspending
force as a function of the gradient of the centripetal potential.
*/

event acceleration (i++) {
  face vector phiacc = a;
  foreach_face ()
    phiacc.x[] += sfm.eps*face_value(f,0)*face_gradient_x(phicentripetal,0);
}

/**
## References

~~~bib
@article{saufi2019dropletsmoke++,
  title={DropletSMOKE++: a comprehensive multiphase CFD framework for the evaporation of multidimensional fuel droplets},
  author={Saufi, Abd Essamade and Frassoldati, Alessio and Faravelli, T and Cuoci, A},
  journal={International Journal of Heat and Mass Transfer},
  volume={131},
  pages={836--853},
  year={2019},
  publisher={Elsevier}
}
~~~
*/
