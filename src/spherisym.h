/**
# Spherically-Symmetric coordinates

This file is adapted from ([spherisym.h](/sandbox/ysaade/allMach/src/spherisym.h)).
It defines the metrics that allows a spherical coordinates system to be
simulated using a 1D grid. The total volume of the sphere can be obtained
by multipliying statsf(f).sum by a factor $4\pi.$ */

#define SPHERISYM 1
#define dimension 1

event metric (i = 0) {

  /**
  By default *cm* is a constant scalar field. To make it variable, we
  need to allocate a new field. We also move it at the begining of the
  list of variables: this is important to ensure the metric is defined
  before other fields. */

  if (is_constant(cm)) {
    scalar * l = list_copy (all);
    cm = new scalar;
    free (all);
    all = list_concat ({cm}, l);
    free (l);
  }

  /**
  The volume/area of a cell is proportional to $r^2$ (i.e. $x*x$). We need
  to set boundary conditions at the top and bottom so that *cm* is
  interpolated properly when refining/coarsening the mesh. */

  scalar cmv = cm;
  foreach()
    cmv[] = x*x;
  cm[left] = dirichlet(x*x);
  cm[right] = dirichlet(x*x);

  /**
  We do the same for the length scale factors. The "length" of faces
  on the axis of revolution is zero ($x=r=0$ on the axis). To avoid
  division by zero we set it to epsilon (note that mathematically the
  limit is well posed). */

  if (is_constant(fm.x)) {
    scalar * l = list_copy (all);
    fm = new face vector;
    free (all);
    all = list_concat ((scalar *){fm}, l);
    free (l);
  }
  face vector fmv = fm;
  foreach_face()
    fmv.x[] = max(x*x, 1e-20);
  
  boundary ({cm, fm});
}
