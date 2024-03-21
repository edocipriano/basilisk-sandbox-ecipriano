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
    if (f[] < F_ERR) {
      double z_frac_mix[NGS];
      for (int jj=0; jj<NGS; jj++) {
        scalar YG = YGList[jj];
        z_frac_mix[jj] = YG[];
      }
      zmix[] = OpenSMOKE_GetMixtureFractionFromMassFractions (z_frac_mix,
          z_frac_fuel, gas_start);
    }
    else {
      zmix[] = 1.;
    }
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

  face vector flameinterface[];
  foreach_face()
    flameinterface.x[] = ((zmix[-1] - zsto[-1])*(zmix[] - zsto[]) < 0.) ? 1. : 0.;

  /**
  We compute the distance of the flame to the interface. These definitions are
  similar to those used in the context of the Ghost Fluid Method applied to
  the VOF field. */

  face vector flamelambda[];
  foreach_face() {
    if (flameinterface.x[])
      flamelambda.x[] = (zmix[-1] - zsto[-1])/(zmix[-1] - zsto[]);
    else
      flamelambda.x[] = 0.;
  }

  /**
  We compute the flame position using `flamelambda`. */

  face vector flamepos[];
  foreach_face() {
    coord xc = {x, y, z};
    double xo = xc.x - 0.5*Delta;
    double xn = xc.x + 0.5*Delta;
    flamepos.x[] = xo + flamelambda.x[]*(xn - xo);
  }

  /**
  We calculate the horizontal flame diameter `hflame`, the lower flame diameter
  `lflame` and the upper flame diameter `uflame`. */

  double hflame = 0.;
  foreach_face(y, reduction(max:hflame)) {
    hflame = max (hflame, flamepos.y[]);
  }

  double lflame = HUGE, uflame = 0.;
  foreach_face(x, reduction(min:lflame) reduction(max:uflame)) {
    lflame = min (lflame, flamepos.x[]);
    uflame = max (uflame, flamepos.x[]);
  }

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
    foreach_dimension()
      ind = (flameinterface.x[1] || flameinterface.x[]) ? ind+1 : 0.;
    flameind[] = (ind) ? chi[] : 0.;
  }

  /**
  We write the post-processing data on the Flame file. */

  extern double d_over_d02, D0;
  double diam = sqrt (d_over_d02)*2.*D0;

  fprintf (fpflame, "%g %g %g %g %g %g %g %g %g\n",
      t, t/sq(D0*1e3), hflame, lflame, uflame,
      hflame/diam, lflame/diam, uflame/diam, Tflame);

}


