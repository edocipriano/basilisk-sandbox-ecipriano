/**
# Solid Fiber Mode
*/

#define THERMAL_FIBER 1

scalar TS[], TSInt[], lambdas[], hs[], lambdafluid[], qS[], mapcells[];
double TSmin = HUGE, TSmax = -HUGE;

typedef enum FiberMaterials {
  quartz, SiC
} fiber_materials;

struct FiberModel {
  double T0;                // Intial temperature
  double rho, cp;           // Physical properties (SI units)
  double diam;              // Fiber diameter
  fiber_materials material; // Material of the solid fiber
  // Variable properties
  double (* hv) (double, double);   // Heat transfer coefficient
  double (* lambdav) (double);      // Thermal conductivity
};

struct FiberModel fiber = {
  .T0 = 300.,
  .rho = 1.,
  .cp = 1.,
  .diam = 120.e-6,
  .lambdav = NULL,
};

double h_nu_law_stagnant (double lambda, double diam) {
  return 0.36*lambda/diam;
}

double lambda_quartz (double T) {
  double lambda = 0.;
  if (T < 300.)
    lambda = 2.14749 - 298.76 / 300. +
    20.72e3 / 300. / 300. - 0.54e6 / 300. / 300. / 300.;
  else if (T >= 300. && T <= 1700.)
    lambda = 2.14749 - 298.76 / T +
    20.72e3 / T / T - 0.54e6 / T / T / T;
  else
    lambda = 2.14749 - 298.76 / 1700. +
    20.72e3 / 1700. / 1700. - 0.54e6 / 1700. / 1700. / 1700.;
  return lambda;
}

double lambda_SiC (double T) {
  double lambda = 0.;
  if (T < 300.)
    lambda = 2.516e-14 * 300 * 300 * 300 * 300 * 300
      - 1.186e-10 * 300.*300.*300.*300.
      + 2.206e-7 *300.*300.*300.
      - 0.0002 * 300 * 300 + 0.08796 * 300 - 9.1062;
  else if (T >= 300. && T <= 1300.)
    lambda = 2.516e-14 * T*T*T*T*T - 1.186e-10 * T*T*T*T
      + 2.206e-7 *T*T*T - 0.0002 * T*T + 0.08796 * T - 9.1062;
  else
    lambda = 2.516e-14 * 1300.*1300.*1300.*1300.*1300.
      - 1.186e-10 * 1300.*1300.*1300.*1300.
      + 2.206e-7 *1300.*1300.*1300. - 0.0002 * 1300.*1300.
      + 0.08796 * 1300. - 9.1062;
  return lambda;
}

event init (i = 0) {
  /**
  We set defaults properties depending on the fiber
  material. */

  switch (fiber.material) {
    case quartz:
      fiber.rho = 2220.;
      fiber.cp  = 760.;
      fiber.lambdav = lambda_quartz;
      fiber.hv = h_nu_law_stagnant;
      break;
    case SiC:
      fiber.rho = 2740.;
      fiber.cp  = 670.;
      fiber.lambdav = lambda_SiC;
      fiber.hv = h_nu_law_stagnant;
      break;
    default:
      fprintf (ferr, "WARNING: unknown fiber material.\n");
  }
#if AXI
  //fiber.diam = 2.*Y0;
  fiber.diam = 0.15e-3;
#endif

  foreach()
    TS[] = T[];
    //TS[] = fiber.T0;

  foreach() {
    qS[] = 0.;
    mapcells[] = 0.;
  }
}

event properties (i++) {
  foreach() {
    lambdas[] = fiber.lambdav (TS[]);
    hs[] = fiber.hv (lambdas[], fiber.diam);
  }
}

event phasechange (i++) {
  foreach_boundary(bottom)
    qS[] = 4./fiber.diam*hs[]*(T[] - TS[])*mapcells[]*0.1;

  foreach() {
#ifdef VARPROP
    slT[] -= qS[]*cm[];
    sgT[] -= qS[]*cm[];
#else
    slT[] -= qS[]/rho1/cp1*cm[];
    sgT[] -= qS[]/rho2/cp2*cm[];
#endif
  }
}

event tracer_diffusion (i++) {
  scalar bcells[];
  foreach() {
    bcells[] = -1.;
    mapcells[] = 0.;
  }

  foreach_boundary(bottom) {
    bcells[] = 1.;
    mapcells[] = 1.;
  }

  face vector thetas[];
  foreach_face()
    thetas.x[] = (bcells[]*bcells[-1] > 0.) ? 1. : 0.;
    //thetas.x[] = (bcells[] == 1. && bcells[-1] == 1.) ? 1. : 0.;

  face vector lambdasf[];
  foreach_face()
    lambdasf.x[] = face_value (lambdas, 0)/fiber.rho/fiber.cp*thetas.x[];

  scalar rs[];
  foreach()
    rs[] = qS[]/fiber.rho/fiber.cp;

  diffusion (TS, dt, D=lambdasf, r=rs);

  foreach_boundary(bottom)
    TSInt[] = (lambdafluid[]*T[] + lambdas[]*TS[])/(lambdas[] + lambdafluid[]);

  foreach_boundary(bottom, reduction(max:TSmax) reduction(min:TSmin)) {
    TSmax = max (TSmax, TS[]);
    TSmin = min (TSmin, TS[]);
  }
}

