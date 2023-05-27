#define VARPROP
#define FILTERED

#include "grid/multigrid.h"
#include "navier-stokes/centered-evaporation.h"
#define ufext uf
#include "two-phase-varprop.h"
#include "opensmoke-properties.h"
#include "tension.h"
#include "evaporation.h"
#include "thermal.h"
#include "view.h"

scalar drhodt[];

int maxlevel = 6;
double R0 = 0.5e-3;
double lambda1, lambda2, cp1, cp2, dhev, TIntVal;
double TL0, TG0, volume0;

int main (void) {
  kinfolder = "materials/n-heptane";

  rho1 = 10., rho2 = 1.;
  mu1 = 1.e-4, mu2 = 1.e-5;
  lambda1 = 0.1, lambda2 = 0.01;
  cp1 = 2000., cp2 = 1000.; dhev = 0.;
  TIntVal = 400., TL0 = 400., TG0 = 300.;

  frhocp1.inverse = false;
  frhocp2.inverse = true;

  f.sigma = 0.03;
  f.tracers = {TL,TG,frhocp1,frhocp2};

  size (3.*R0);
  init_grid (1 << maxlevel);

  //TOLERANCE = 1.e-2;
  run();
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y));

event init (i = 0) {
  fraction (f, -circle(x,y,R0));
  volume0 = 0.;
  foreach(reduction(+:volume0))
    volume0 += (1. - f[])*dv();
  foreach() {
    TL[] = f[]*TL0;
    TG[] = (1. - f[])*TG0;
    T[] = TL[] + TG[];
  }
}

event logfile (i++) {
  double volume = 0.;
  foreach(reduction(+:volume))
    volume += (1. - f[])*dv();
  fprintf (stdout, "%g %g\n", t, volume/volume0);
  fflush (stdout);
}

event movie (i += 50) {
  clear();
  view (tx = -0.5, ty = -0.5);
  draw_vof ("f", lw = 1.5);
  squares ("T", min = min (TL0,TG0), max = max (TL0,TG0),
      map=blue_white_red);
  save ("movie.mp4");
}

event snapshots (t += 0.01; t = end) {
  char name[80];
  sprintf (name, "snapshot-%g", t);
  dump (name);
}

event stop (t = 0.01);

