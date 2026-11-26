/**
# Boiling

Temperature-gradient-driven phase change of pure liquids. The interface is at
the saturation temperature *TIntVal*, and the vaporization rate follows from
the interface energy balance
$$
  \dot{m} = \dfrac{\lambda_l\nabla T_l\cdot\mathbf{n}_\Gamma
    - \lambda_g\nabla T_g\cdot\mathbf{n}_\Gamma}{\Delta h_{ev}}
$$
where $\mathbf{n}_\Gamma$ is the interface normal pointing towards the liquid,
and $\dot{m}$ is positive when the liquid loses mass: the heat conducted from
the superheated phase towards the interface is used for the phase change. This
module includes [phasechange.h](phasechange.h), computes the vaporization
rate and the interface heat fluxes, and sets the phase change model
parameters (`pcm`) for boiling.

With `IMPLICIT`, the saturation temperature is imposed implicitly in the
diffusion (see *plic_flux()* in [gradients.h](gradients.h)). */

#ifndef IMPLICIT
# define IMPLICIT 0
#endif

#include "phasechange.h"
#include "gradients.h"

scalar mEvap[], * mEvapList = {mEvap};
scalar TInt[];

/**
## Phase change

The volume and face fractions of each phase are restricted to the fluid
region with embedded boundaries. */

void phasechange_rate (void) {
  phase_tracers_to_scalars (liq, f, tol = F_ERR);
  phase_tracers_to_scalars (gas, f, tol = F_ERR);

  foreach()
    TInt[] = TIntVal;

  scalar fl[], fg[];
  scalar_clone (fl, f);
  scalar_clone (fg, f);

  foreach() {
    fl[] = f[];
    fg[] = 1. - f[];
#if EMBED
    fl[] *= cs[];
    fg[] *= cs[];
#endif
  }
  fg.inverse = true;

  face vector fsl[], fsg[];
  face_fraction (fl, fsl);
  face_fraction (fg, fsg);

  //// fixme as:
  //// face vector fsl[], fsg[];
  //// face_fraction (f, fsl, inverse = false);
  //// face_fraction (f, fsg, inverse = true);

  /**
  The function *plic_gradient()* returns the derivative along the normal
  entering each phase, i.e. $\nabla T_l\cdot\mathbf{n}_\Gamma$ for the liquid
  and $-\nabla T_g\cdot\mathbf{n}_\Gamma$ for the gas: the two contributions
  are therefore summed. The rate is set to zero where the enthalpy of
  evaporation is not positive. */

  foreach_interfacial (f, F_ERR) {
    scalar TL = liq->T, TG = gas->T;
    bool vof = vof_gradients;
#if IMPLICIT
    double ltrgrad = plic_gradient_implicit (point, TL, fl, fsl, TIntVal, vof,
        NULL);
    double gtrgrad = plic_gradient_implicit (point, TG, fg, fsg, TIntVal, vof,
        NULL);
#else
    double ltrgrad = plic_gradient (point, TL, fl, fsl, TIntVal, vof, NULL);
    double gtrgrad = plic_gradient (point, TG, fg, fsg, TIntVal, vof, NULL);
#endif

    scalar lambdal = liq->lambda, lambdag = gas->lambda;
    scalar deltahev = liq->dhev;

    mEvap[] += (deltahev[] > 0.) ? lambdal[]*ltrgrad/deltahev[] : 0.;
    mEvap[] += (deltahev[] > 0.) ? lambdag[]*gtrgrad/deltahev[] : 0.;
  }

  /**
  Without `IMPLICIT`, the interface conductive heat fluxes enter the
  temperature equation of each phase as explicit source terms. */

#if !IMPLICIT
  foreach_interfacial_plic (f, F_ERR) {
    scalar TL = liq->T, TG = gas->T;
    bool vof = vof_gradients;
    double ltrgrad = plic_gradient (point, TL, fl, fsl, TIntVal, vof, NULL);
    double gtrgrad = plic_gradient (point, TG, fg, fsg, TIntVal, vof, NULL);

    scalar slT = liq->STexp, sgT = gas->STexp;
    scalar lambdal = liq->lambda, lambdag = gas->lambda;

    slT[] += lambdal[]*ltrgrad*dirac;
    sgT[] += lambdag[]*gtrgrad*dirac;
  }
#endif

  phase_scalars_to_tracers (liq, f);
  phase_scalars_to_tracers (gas, f);
}

/**
## Order of integration

The boiling problem is difficult from the evaporation, since the interfacial temperature
is often known a-priori. This allows using a different time integration scheme. Below, we
overwrite the default time integration of [phasechange.h](phasechange.h). */

#if IMPLICIT
macro phasechange_vof() {
  phasechange_transport();
}

macro phasechange_tracer_diffusion() {
  phasechange_reset();
  phasechange_diffusion();
  phasechange_rate();
  phasechange_expansion();
  phasechange_regression();
  phasechange_properties();
  phasechange_divergence();
  phasechange_reset_sources();
}
#endif

/**
## Phase change model parameters

We set the `pcm` parameters for boiling: the expansion term is shifted to the
gas, `pcm.boiling` selects the boiling treatment of the expansion term, and
the interface regression velocity uses the liquid density. */

// fixme: check compatibility with restore
event init (i = 0) {
  pcm.shifting = SHIFT_TO_GAS;
  pcm.boiling = true;
  pcm.byrhogas = false;

  /**
  The interface temperatures are the interfacial values of the temperatures of
  the phases (see *phasechange_tracers()* in [regression.h](regression.h)). */

  TInt.nodump = true;
  foreach()
    TInt[] = TIntVal;
  if (!liq->isothermal) {
    scalar TL = liq->T;
    TL.interfacial = TInt;
  }
  if (!gas->isothermal) {
    scalar TG = gas->T;
    TG.interfacial = TInt;
  }
}
