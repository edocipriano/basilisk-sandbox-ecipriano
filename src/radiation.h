/**
# Radiation Models

Radiation models for combustion simulations. Currently,
just the optically-thin model is implemented.
This implementation is based on the information reported
in [https://tnfworkshop.org/radiation/](https://tnfworkshop.org/radiation/)
*/

#define RADIATION

/**
We define the Stefan-Boltzmann constant $[W/m^2/K^4]$. */

#define STEFAN_BOLTZMANN  5.669e-8

/**
## Optically-Thin Model

This is the simplest model to implement, it computes the
heat loss for radiation according to the following equation:
$$
\nabla\cdot \mathbf{q}_{rad} =
-4\sigma\sum\limits_{i=1}^{NGS}a_{p,i} P x_i (T^4 - T_b^4)
$$
where $\sigma$ is the Stefan-Boltzmann constant, $a_{p,i}$
is the Planck mean absorption coefficient of species $i$,
$P x_i$ is the partial pressure of the chemical species,
$T$ is the local flame temperature, while $T_b$ is a baseline
temperature (approximated as 300K).

The following properties are required by the radiation
model:

* *T*: temperature of the system [K]
* *P*: pressure of the system [Pa]
* *xCO2*: mole fraction of CO2 [moli/molt]
* *xH2O*: mole fraction of H2O [moli/molt]

Only CO2 and H2O are considered, since their contribution to
the radiation is the gratest, if compared with the other
species.
*/

typedef struct {
  double T, P;
  double xCO2, xH2O;
} OpticallyThinProperties;

/**
We also declare a struct that gathers optically-thin settings. */

struct OpticalyThinModel {
  int indexCO2, indexH2O;
};

struct OpticalyThinModel otm = {
  .indexCO2 = -1,
  .indexH2O = -1,
};

/**
We compute the radiation according to the optically-thin model. */

double optically_thin (void * p) {

  /**
  We extract the data from the void pointer. */

  OpticallyThinProperties * otp;
  otp = (OpticallyThinProperties *)p;

  double T = otp->T, P = otp->P;
  double xCO2 = otp->xCO2, xH2O = otp->xH2O;

  /**
  The pressure value in SI (Pa) is converted in atm. */

  P *= 9.86923e-6;

  /**
  We compute the Plank-mean absorption coefficients using
  coefficients from [NIST](#grosshandler1993radcal). */

  double apCO2 = 18.741 - 121.310*(1000./T) + 273.500*pow (1000./T, 2.)
    - 194.050*pow (1000./T, 3.) + 56.310*pow (1000./T, 4.) - 5.8169*pow (1000./T, 5.);

  double apH2O = - 0.23093 - 1.12390*(1000./T) + 9.41530*pow (1000./T, 2.)
    - 2.99880*pow (1000./T, 3.) + 0.51382*pow (1000./T, 4.) - 1.86840e-5*pow (1000./T, 5.);

  double sum_pa = (apCO2*xCO2 + apH2O*xH2O)*P;

  /**
  We return the flux of radiant energy. */

  return 4.*STEFAN_BOLTZMANN*sum_pa*( pow (T, 4.) - pow (300., 4.) );
}

/**
## No Radiation

Implementation of a dummy no_radiation model, which returns
a 0 divq contribution from the radiation.
*/

double no_radiation (void * p) {
  return 0;
}

/**
## Main radiation function

We declare a function pointer to a generic radiation
model that computes the term $\nabla\cdot q_{rad}$
for the energy equation.
*/

double (* divq_rad) (void * p) = no_radiation;

/**
## Init Event

We overload the init event in order to setup properties
for the radiation models at the beginning of the simulation.
*/

event init (i = 0) {
  if (divq_rad == optically_thin) {

    for (int jj=0; jj<NGS; jj++) {
      if (strcmp (gas_species[jj], "CO2") == 0)
        otm.indexCO2 = jj; 
      if (strcmp (gas_species[jj], "H2O") == 0)
        otm.indexH2O = jj; 
    }

    // Species not found
    assert (otm.indexCO2 != -1);
    assert (otm.indexH2O != -1);
  }
}

/**
## References

~~~bib
@article{grosshandler1993radcal,
  title={Radcal: A narrow band model for radiation},
  author={Grosshandler, William L},
  journal={Calculations in a Combustion Environment, NIST Technical Note},
  volume={1402},
  year={1993}
}
~~~

*/

