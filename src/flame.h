/**
# Post-Processing of Flame Properties

This module is conceived to automatically compute and print interesting flame
properties, such as the mixture fraction, scalar dissipation rate, integral
quantities and any other flame property.

We assume that we are using OpenSMOKE to manage kinetics, thermodynamics, and
transport properties. */

#include "OpenSMOKE_Interface.h"

scalar zmix[], zsto[], chi[];
scalar flameind[];

FILE * fpflame;

event init (i = 0) {
  char name[80];
  sprintf (name, "Flame-%d", grid->maxdepth);
  if (pid() == 0)
    fpflame = fopen (name, "w");
}

event flame (i++) {

  /**
  We compute the mixture fraction using Bilger's formula. */

  double z_frac_fuel[NGS], z_frac_stoi[NGS];
  for (int jj=0; jj<NGS; jj++) {
    scalar YGInt = YGIntList[jj];
    z_frac_fuel[jj] = avg_interface (YGInt, f, tol=0.1);
    z_frac_stoi[jj] = 0.;
  }

  foreach() {
    double z_frac_mix[NGS];
    for (int jj=0; jj<NGS; jj++) {
      scalar Y = YList[jj];
      z_frac_mix[jj] = Y[];
    }
    zmix[] = OpenSMOKE_GetMixtureFractionFromMassFractions (z_frac_mix,
        z_frac_fuel, gas_start);

    zsto[] = OpenSMOKE_GetMixtureFractionFromMassFractions (z_frac_stoi,
        z_frac_fuel, gas_start);
  }

  vector gzmix[];
  gradients ({zmix}, {gzmix});

  scalar gzmixmag[];
  foreach() {
    double mag = 0.;
    foreach_dimension()
      mag += sq (gzmix.x[]);
    gzmixmag[] = sqrt (mag);
  }

  /**
  We calculate the scalar dissipation rate given the mixture fraction:

  $$
    \chi = 2\alpha_v|\nabla Z|^2
  $$
  */

  foreach() {
    double lambda2vh = lambda2, rho2vh = rho2, cp2vh = cp2;
#ifdef VARPROP
    lambda2vh = lambda2v[], rho2vh = rho2v[], cp2vh = cp2v[];
#endif
    chi[] = (rho2vh*cp2vh > 0.) ?
      2.*lambda2vh/rho2vh/cp2vh*sq (gzmixmag[]) : 0.;
  }

  /**
  We compute the flame interface field. */

  scalar zdiff[];
  foreach()
    zdiff[] = zmix[] - zsto[];

  face vector flameinterface[];
  foreach_face()
    flameinterface.x[] = (zdiff[]*zdiff[-1] < 0.) ? 1. : 0.;

  /**
  We compute the distance of the flame to the interface. The flame interface
  is the zero level set of the `zdiff` field. Therefore, we calculate the
  coordinates of the flame using definitions which are similar to those of
  the level set method. */

  face vector flamelambda[];
  foreach_face()
    flamelambda.x[] = (flameinterface.x[] == 1.) ? (zdiff[-1] - 0.)/(zdiff[-1] - zdiff[]) : 0.;

  /**
  We compute the flame position using `flamelambda`. */

  vector xp[];
  foreach() {
    coord o = {x, y, z};
    foreach_dimension()
      xp.x[] = o.x;
  }

  face vector flamepos[];
  foreach_face()
    flamepos.x[] = (flameinterface.x[] == 1.) ? xp.x[-1] + flamelambda.x[]*Delta : 0.;

  /**
  We calculate the horizontal flame diameter `Dx`, the vertical flame
  diameter `Dy`, and the effective diameter. */

  double Dx = 2.*statsf(flamepos.y).max;
  double Dy = statsf(flamepos.x).max - statsf(flamepos.x).min;
  if (!X0)
    Dy = 2.*statsf(flamepos.x).max;
  double De = 0.5*(Dx + Dy);

  /**
  Update flame temperature. */

  double nfaces = 0;
  double Tflame = 0.;
  foreach_face(reduction(+:nfaces) reduction(+:Tflame)) {
    if (flameinterface.x[]) {
      Tflame += TG[-1]*(1. - flamelambda.x[]) + TG[]*flamelambda.x[];
      nfaces++;
    }
  }
  Tflame *= (nfaces > 0) ? 1./nfaces : 0.;

  /**
  Interpolate the scalar dissipation rate on the flame. */

  foreach() {
    int ind = 0;
    foreach_dimension() {
      if (flameinterface.x[1] || flameinterface.x[])
        ind++;
    }
    flameind[] = (ind) ? 1. : 0.;
  }

  extern double d_over_d02, D0;
  double diam = sqrt (d_over_d02)*D0;
  double statTmax = statsf(T).max;

#if 0
  /**
  Flame diameter using the temperature peak. */

  Array * arrtemp = array_new();
  for (double y = 0.; y < L0; y += 0.5*L0/(double)N) {
    double val = interpolate (T, 0., y);
    val = (val == nodata) ? 0. : val;
    array_append (arrtemp, &val, sizeof(double));
  }
  double * temps = (double *)arrtemp->p;

  @if _MPI
  int size = arrtemp->len/sizeof(double);
  MPI_Allreduce (MPI_IN_PLACE, temps, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  @endif

  int index = 0;
  double Tmax = 0., ymax = 0.;
  for (double y = 0.; y < L0; y += 0.5*L0/(double)N) {
    if (temps[index] > Tmax) {
      Tmax = temps[index];
      ymax = y;
    }
    index++;
  }
  double DTmax = 2.*ymax;
  array_free (arrtemp);
#else
  /**
  The following calculation of the temperature peak is
  correct only for 1/4 droplet in microgravity conditions. */

  double Tmax = statTmax, radius = 0.;
  foreach(reduction(max:radius)) {
    if (T[] == Tmax)
      radius = sqrt (sq(x) + sq(y));
  }
  double DTmax = 2.*radius;

#endif

  /**
  We write the post-processing data on the Flame file. */

  if (pid() == 0) {
    fprintf (fpflame, "%g %g %g %g %g %g %g %g %g %g %g %g %g\n",
        t, t/sq(D0*1e3), Dx, Dy, De,
        Dx/diam, Dy/diam, De/diam, Tflame, statTmax, Tmax, DTmax, DTmax/diam);
  }
}


