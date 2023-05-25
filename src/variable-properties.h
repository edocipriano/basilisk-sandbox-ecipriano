/**
# Variable Properties

Simulations involving two-phase flows with variable material
properties can be performed using this module, which defines
structures and functions that help to setup such cases.
*/

#define VARPROP

/**
## Fields Allocations

We allocate fields required by this module. All properties
are intialized as constant fields, and they are initialized
just by the solver that needs them. If a thermal solver is
used, there is no need to initialize a non-constant
diffusivity of the chemical species, for example.
*/

(const) scalar rho1v = zeroc, rho2v = zeroc;
(const) scalar mu1v = zeroc, mu2v = zeroc;
(const) face vector lambda1v = zerof, lambda2v = zerof;
(const) scalar cp1v = zeroc, cp2v = zeroc;

/**
## Thermodynamic State

We define a structure with variables that define the thermodynamic
state of the mixture: temperature *T*, pressure *P*, composition *x*.
*/

typedef struct {
  double T, P;
  double * x;
} ThermoState;

/**
## Thermodynamic Properties

We define a struct that contains a bunch of function pointers
for each property of interest for the gas phase.
*/

typedef struct {
  // Mixture properties
  double (* rhov) (void *);
  double (* muv) (void *);
  double (* lambdav) (void *);
  double (* cpv) (void *);
  // Species properties
  double (* pvap) (void *, int);
  double (* dhev) (void *, int);
  double (* diff) (void *, int);
} ThermoProps;

