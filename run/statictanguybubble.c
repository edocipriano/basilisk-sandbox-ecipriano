#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "view.h"

double X0, Y0;
double R0 = 1.e-3;
int maxlevel = 7;
double tEnd = 0.1;

int main (void) {
  rho1 = 958., rho2 = 0.59;
  mu1 = 2.82e-4, mu2 = 1.23e-6;
  f.sigma = 0.059;

  L0 = 8.*R0;
  X0 = 0.5*L0, Y0 = 0.;
  origin (0., 0.);
  size (L0);
  init_grid (1 << maxlevel);
  run();
}

#define circle(x,y,R) (sq(x - 0.5*L0) + sq(y - 0.) - sq(R))

event init (i = 0) {
  fraction (f, circle(x,y,R0));
}

event movie (t += 0.1) {
  clear ();
  view (tx = -0.5, ty = -0.5);
  draw_vof ("f");
  squares ("f", min = 0., max = 1., linear = false);
  save ("t.mp4");
}

event snapshots (t += 0.1) {
  char name[80];
  sprintf (name, "snapshot-%g", t);
  dump (name);
}

event stop (t = tEnd) {
  return 0;
}
