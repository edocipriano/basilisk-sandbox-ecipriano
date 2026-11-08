/**
# Evaporation

Diffusion-driven phase change of multicomponent mixtures in non-isothermal
conditions, following [Cipriano et al., 2023](#cipriano2023), and [Cipriano
et al., 2024](#cipriano2024) for variable properties. This module
includes [phasechange.h](phasechange.h), computes the vaporization rates, and
sets the phase change model parameters (`pcm`) for evaporation. In each
interfacial cell, it solves the interface jump condition, assuming
thermodynamic equilibrium at the interface:
$$
  \begin{cases}
    \dot{m}_i =
    \dot{m}\hat\omega_{i,l}
    + \rho_l\mathcal{D}_{i,l}\nabla\omega_{i,l}\cdot\mathbf{n}_\Gamma
    =
    \dot{m}\hat\omega_{i,g}
    + \rho_g\mathcal{D}_{i,g}\nabla\omega_{i,g}\cdot\mathbf{n}_\Gamma
    \\
    \sum_i \dot{m}_i\Delta h_{ev,i} =
    \lambda_l\nabla T_l\cdot\mathbf{n}_\Gamma
    - \lambda_g\nabla T_g\cdot\mathbf{n}_\Gamma
    \\
    \hat{x}_{i,g} = K_{eq,i}(\hat{T})\,\hat{x}_{i,l}
  \end{cases}
$$
where the hat denotes interface values, and $\mathbf{n}_\Gamma$ is the
interface normal pointing towards the liquid. The equilibrium is written here
for evaporation; in general, it relates the receiver to the carrier mole
fractions (see below). The solution gives the vaporization rate
$\dot{m}_i$ of each species, the interface mass fractions
$\hat\omega_{i,l}$ and $\hat\omega_{i,g}$, and the interface temperature
$\hat{T}$, which provide the source terms of the species and temperature
equations solved by [phasechange.h](phasechange.h).

## Carrier and receiver

The two phases play different roles in the phase change:

* the *carrier* contains only species that cross the interface;
* the *receiver* contains the same species, plus *inert* species that never
  cross the interface.

In droplet evaporation, the liquid fuel is the carrier, and the gas (fuel
vapor plus air, the inert) is the receiver. In bubble dissolution, or bubble
growth during electrolysis, the gas is the carrier, and the liquid (dissolved
gas plus water, the inert) is the receiver. The mass flux goes from the
carrier to the receiver: it is negative if species move the other way, as in
a growing bubble.

The carrier is therefore the phase with fewer species (the liquid if `NLS`
and `NGS` are equal), and each carrier species must also be a receiver
species. The following variables describe the two sets of species:

* *NCS*: number of carrier species (all the species of the carrier phase);
* *NRS*: number of receiver species;
* *NIS*: number of inert species, *NRS - NCS*;
* *CSI*: carrier species indices, i.e. the receiver index of each carrier
  species;
* *ISI*: inert species indices, i.e. the receiver index of each inert species.

The roles are chosen once, in *liquid_carrier*. The species part of the jump
condition is computed by *species_jump_condition()*, which takes the carrier
and the receiver as arguments: it is called with the liquid as the carrier for
evaporation, and with the gas as the carrier for bubble dissolution.

The vaporization rates `mEvapList` are indexed like the receiver species. The
sign convention is the one in [phasechange.h](phasechange.h): `mEvapTot` and
`mEvapList` are positive when the liquid loses mass. */

#include "phasechange.h"
#include "gradients.h"
#include "fsolve.h"

Phase * liq_int, * gas_int;

scalar * mEvapList = NULL;

bool liquid_carrier = true;
int NCS = 0, NRS = 0, NIS = 0;
int * CSI = NULL, * ISI = NULL;

double * YIntVals;

/**
## Interface energy balance */

typedef struct {
  Phase * liq, * gas, * carrier;
  scalar fl, fg;
  face vector fsl, fsg;
  coord c;
} EnergyBalanceData;

double divq_rad_int (double TInti, double Tbulk = 300., double alphacorr = 0) {
  return alphacorr*5.670373e-8*(pow (Tbulk, 4.) - pow (TInti, 4.));
}

void energy_balance (const double * xdata, double * fdata, void * params) {
  EnergyBalanceData * data = params;

  Phase * liq = data->liq, * gas = data->gas, * c = data->carrier;
  scalar fl = data->fl, fg = data->fg;
  face vector fsl = data->fsl, fsg = data->fsg;

  scalar TL = liq->T, TG = gas->T;

  foreach_point (data->c.x, data->c.y, data->c.z, serial) {
    double TInti = xdata[0];

    double gtrgrad = plic_gradient (point, TG, fg, fsg, TInti, true, NULL);
    double ltrgrad = plic_gradient (point, TL, fl, fsl, TInti, true, NULL);

    scalar lambda1 = liq->lambda, lambda2 = gas->lambda;

    double vapheat = 0.;
    for (size_t i = 0; i < c->n; i++) {
      scalar mEvap = mEvapList[CSI[i]];
      scalar dhev = c->dhevList[i];
      vapheat -= mEvap[]*dhev[];
    }

    fdata[0] = vapheat
      - divq_rad_int (TInti, TG0, pcm.emissivity)
      + lambda1[]*ltrgrad
      + lambda2[]*gtrgrad
      ;
  }
}

#if USE_GSL
int energy_balance_gsl (const gsl_vector * x, void * params, gsl_vector * f) {
  double * xdata = x->data, * fdata = f->data;
  energy_balance (xdata, fdata, params);
  return GSL_SUCCESS;
}
#endif

/**
## Interface equilibrium

The function *antoine()* returns the equilibrium constant of the i-th carrier
species, i.e. the ratio between its receiver-side and carrier-side interface
mole fractions. For evaporation it is Raoult's law, $p^{sat}_i(T)/P$; for
bubble dissolution it is Henry's law, $P/H_i(T)$. The default implementation
uses the *antoine* attribute of the carrier mass fraction, if set. */

double antoine_default (double T, double P, int i) {
  if (YIntVal)
    return YIntVal;
  else {
    Phase * c = liquid_carrier ? liq : gas;
    scalar YC = c->YList[i];
    if (YC.antoine)
      return YC.antoine (T, P);
    else
      return YIntVals[i];
  }
}

double (* antoine) (double, double, int) = &antoine_default;

/**
## Species jump condition

Normal gradient, at the interface, of the receiver species *k*, from mass or
mole fractions. */

double receiver_gradient (Point point, Phase * r, Phase * r_int, int k,
    scalar fr, face vector fsr)
{
  double grad = 0.;
  if (pcm.molar_diffusion) {
    scalar X = r->XList[k], XInt = r_int->XList[k], MWInt = r_int->MW;
    grad = plic_gradient (point, X, fr, fsr, XInt[], true, NULL);
    grad *= (MWInt[] > 0) ? r_int->MWs[k]/MWInt[] : 0.;
  }
  else {
    scalar Y = r->YList[k], YInt = r_int->YList[k];
    grad = plic_gradient (point, Y, fr, fsr, YInt[], true, NULL);
  }
  return grad;
}

/**

The function below solves the species part of the jump condition, from the
carrier *c* to the receiver *r*, with the decoupled procedure of [Cipriano et
al., 2023](#cipriano2023), Sec. 3.5: the carrier is assumed to evolve more
slowly than the receiver (e.g. the liquid during evaporation), so its interface
values are taken from the surrounding cells, while the receiver side follows
from equilibrium and mass balances. *c_int* and *r_int* store the interface
values, *fr* and *fsr* are the volume and face fractions of the receiver, and
*sign* converts the carrier-to-receiver mass flux into vaporization rates (+1
if the carrier is the liquid). The interface temperature must be already set
in *c_int*. */

void species_jump_condition (Phase * c, Phase * r, Phase * c_int,
    Phase * r_int, scalar fr, face vector fsr, double sign)
{

  /**
  The carrier interface mass fractions are the volume-fraction-weighted
  averages over the 3x3 stencil (Eq. 45). The receiver values only provide
  the inert composition, rescaled below. */

  foreach_interfacial (f, F_ERR) {
    foreach_species_in (c_int)
      Y[] = avg_neighbor (point, c->YList[i], f);
    foreach_species_in (r_int) {
      scalar YR = r->YList[i];
      Y[] = YR[];
    }
  }

  // Update mole fractions and moles
  memcpy (c_int->MWs, c->MWs, c->n*sizeof (double));
  memcpy (r_int->MWs, r->MWs, r->n*sizeof (double));
  phase_update_mw_moles (c_int);
  phase_update_mw_moles (r_int);

  /**
  The receiver interface mole fractions follow from the phase equilibrium
  (Eq. 47)
  $$
  \hat{x}_{i,r} = K_{eq,i}(\hat{T})\,\hat{x}_{i,c}
  $$
  with $K_{eq,i}$ given by *antoine()*. They are converted to mass fractions
  by lumping the inert species into a single pseudo-species (the last element
  of the arrays below), which closes the mole fractions to one. The inert
  species are then rescaled, preserving their relative composition. */

  scalar TInt = c_int->T;
  foreach_interfacial (f, F_ERR, serial) {
    double yrcorr[NCS + 1], xrcorr[NCS + 1], mwcorr[NCS + 1];
    double yinert[NIS > 0 ? NIS : 1], mwinert[NIS > 0 ? NIS : 1];

    // Store the sum of the inert species
    double sum_yinert_old = 0., sum_xinert_old = 0.;
    for (int j = 0; j < NIS; j++) {
      scalar YRInt = r_int->YList[ISI[j]];
      scalar XRInt = r_int->XList[ISI[j]];
      sum_yinert_old += YRInt[];
      sum_xinert_old += XRInt[];
    }

    // Thermodynamic equilibrium
    for (int i = 0; i < NCS; i++) {
      scalar XCInt = c_int->XList[i];
      scalar XRInt = r_int->XList[CSI[i]];
      XRInt[] = min (antoine (TInt[], Pref, i), 0.98)*XCInt[];
    }

    // Molecular weight of the inert species mixture
    for (int j = 0; j < NIS; j++) {
      scalar YRInt = r_int->YList[ISI[j]];
      yinert[j] = YRInt[];
      mwinert[j] = r_int->MWs[ISI[j]];
    }
    correctfrac (yinert, NIS);
    double mwinertmix = mass2mw (yinert, mwinert, NIS);

    // Convert interface mole fractions to mass fractions
    double sum_xr = 0.;
    for (int i = 0; i < NCS; i++) {
      scalar XRInt = r_int->XList[CSI[i]];
      sum_xr += XRInt[];
      xrcorr[i] = XRInt[];
      mwcorr[i] = c_int->MWs[i];
    }
    xrcorr[NCS] = 1. - sum_xr;
    mwcorr[NCS] = mwinertmix;
    mole2massfrac (yrcorr, xrcorr, mwcorr, NCS + 1);

    for (int i = 0; i < NCS; i++) {
      scalar YRInt = r_int->YList[CSI[i]];
      YRInt[] = yrcorr[i];
    }

    // Rescale the inert species
    for (int j = 0; j < NIS; j++) {
      scalar YRInt = r_int->YList[ISI[j]];
      scalar XRInt = r_int->XList[ISI[j]];
      YRInt[] *= yrcorr[NCS] / (sum_yinert_old + 1e-10);
      XRInt[] *= xrcorr[NCS] / (sum_xinert_old + 1e-10);
    }
  }

  /**
  Since the carrier diffusive fluxes sum to zero, summing the species
  balances over the carrier species gives an explicit expression for the
  total mass flux from the carrier to the receiver (Eqs. 49, 50)
  $$
  \dot{m} = \frac{\sum_{i=1}^{NCS} j_{i,r}}
    {1 - \sum_{i=1}^{NCS} \hat{\omega}_{i,r}}
  \quad
  \dot{m}_i = \dot{m}\,\hat{\omega}_{i,r} + j_{i,r}
  $$
  where $j_{i,r}$ is the diffusive flux of the receiver species along the
  normal $\mathbf{n}_r$ entering the receiver. It is defined by two options
  ([Cipriano et al., 2024](#cipriano2024), Eqs. 26, 27):

  * `pcm.molar_diffusion`: the diffusivities are molar based (e.g. from the
  Hirschfelder and Curtiss approximation), so Fick's law is written in terms of
  mole fractions, $j^F_{i,r} = -\rho_r\mathcal{D}_{i,r} (M_{w,i}/M_{w,r})\,
  \partial x_{i,r}/\partial\mathbf{n}_r$, with $M_{w,r}$ the molecular weight
  of the interface mixture. Otherwise, $j^F_{i,r} = -\rho_r\mathcal{D}_{i,r}\,
  \partial\omega_{i,r}/\partial\mathbf{n}_r$.
  * `pcm.fick_corrected`: Fick's law is exact only for binary mixtures, or for
  equal diffusivities. The correction
  $j_{i,r} = j^F_{i,r} - \hat{\omega}_{i,r}\sum_{k=1}^{NRS} j^F_{k,r}$
  makes the diffusive fluxes sum to zero. Otherwise, $j_{i,r} = j^F_{i,r}$.
  */

  foreach_interfacial (f, F_ERR) {
    scalar rhoR = r->rho;

    double jRtot = 0.;
    if (pcm.fick_corrected) {
      for (int k = 0; k < NRS; k++) {
        scalar DR = r->DList[k];
        double rtrgrad = receiver_gradient (point, r, r_int, k, fr, fsr);
        jRtot += -rhoR[]*DR[]*rtrgrad;
      }
    }

    double sum_jR = 0., sum_YRInt = 0.;
    for (int i = 0; i < NCS; i++) {
      scalar YRInt = r_int->YList[CSI[i]];
      scalar DR = r->DList[CSI[i]];
      double rtrgrad = receiver_gradient (point, r, r_int, CSI[i], fr, fsr);
      sum_jR += -rhoR[]*DR[]*rtrgrad - YRInt[]*jRtot;
      sum_YRInt += YRInt[];
    }
    double mtot = sum_jR / min (1. - sum_YRInt, 0.99);
    mEvapTot[] = sign*mtot;

    for (int i = 0; i < NCS; i++) {
      scalar mEvap = mEvapList[CSI[i]];
      scalar YRInt = r_int->YList[CSI[i]];
      scalar DR = r->DList[CSI[i]];
      double rtrgrad = receiver_gradient (point, r, r_int, CSI[i], fr, fsr);
      mEvap[] = sign*(mtot*YRInt[] - rhoR[]*DR[]*rtrgrad - YRInt[]*jRtot);
    }
  }

  /**
  The species equations get the explicit source $\dot{m}_i\delta_\Gamma$ and
  the implicit source $\dot{m}\,\omega_i\delta_\Gamma$, which accounts for
  the volume change (Eq. 43; both explicit with `EXPLICIT_ONLY`). Since
  vaporization rates are positive when the liquid loses mass, the sign of the
  sources is *sign* for the carrier and *-sign* for the receiver. */

  foreach_interfacial_plic (f, F_ERR) {
    for (int k = 0; k < NRS; k++) {
      scalar mEvap = mEvapList[k], YInt = r_int->YList[k];
      scalar SYexp = r->SYexpList[k], SYimp = r->SYimpList[k];
      if (pcm.diffusion == EXPLICIT_ONLY)
        SYexp[] += -sign*(mEvap[] - mEvapTot[]*YInt[])*dirac;
      else {
        SYexp[] += -sign*mEvap[]*dirac;
        SYimp[] += sign*mEvapTot[]*dirac;
      }
    }

    for (int i = 0; i < NCS; i++) {
      scalar mEvap = mEvapList[CSI[i]], YInt = c_int->YList[i];
      scalar SYexp = c->SYexpList[i], SYimp = c->SYimpList[i];
      if (pcm.diffusion == EXPLICIT_ONLY)
        SYexp[] += sign*(mEvap[] - mEvapTot[]*YInt[])*dirac;
      else {
        SYexp[] += sign*mEvap[]*dirac;
        SYimp[] += -sign*mEvapTot[]*dirac;
      }
    }
  }
}

/**
## Setup */

event defaults (i = 0) {
  liq_int = new_phase_minimal ("LInt", NLS, false, liq_species);
  gas_int = new_phase_minimal ("GInt", NGS, true, gas_species);

  liquid_carrier = (NLS <= NGS);
  Phase * c_int = liquid_carrier ? liq_int : gas_int;
  Phase * r_int = liquid_carrier ? gas_int : liq_int;

  NCS = c_int->n;
  NRS = r_int->n;
  NIS = NRS - NCS;

  foreach_species_in (r_int) {
    scalar mEvap = new scalar;
    char name[80];
    sprintf (name, "mEvap%s", r_int->species[i]);
    free (mEvap.name);
    mEvap.name = strdup (name);
    mEvap.nodump = true;
    mEvapList = list_add (mEvapList, mEvap);
  }

  // Carrier species indices *CSI* (aborts if missing from the receiver)
  CSI = (int *)malloc (NCS*sizeof (int));
  for (int i = 0; i < NCS; i++)
    CSI[i] = phase_species_index (r_int, c_int->species[i]);

  // Inert species indices *ISI*: the receiver species not in *CSI*
  bool carrier_species[NRS];
  for (int k = 0; k < NRS; k++)
    carrier_species[k] = false;
  for (int i = 0; i < NCS; i++)
    carrier_species[CSI[i]] = true;

  ISI = (int *)malloc (NIS*sizeof (int));
  int nis = 0;
  for (int k = 0; k < NRS; k++)
    if (!carrier_species[k])
      ISI[nis++] = k;
  assert (nis == NIS);

  // Default equilibrium constants of the carrier species
  YIntVals = (double *)calloc (NCS, sizeof (double));

  // We dump only non-null vaporization rates
  for (int i = 0; i < NCS; i++) {
    scalar mEvap = mEvapList[CSI[i]];
    mEvap.nodump = false;
  }
}

event cleanup (t = end) {
  delete_phase (liq_int), liq_int = NULL;
  delete_phase (gas_int), gas_int = NULL;
  delete (mEvapList), free (mEvapList), mEvapList = NULL;
  free (CSI), CSI = NULL;
  free (ISI), ISI = NULL;
  free (YIntVals);
}

event reset_sources (i++) {
  foreach() {
    foreach_scalar_in (liq_int) {
      T[] = 0.;
      foreach_species_in (liq_int) {
        Y[] = 0.;
        X[] = 0.;
      }
    }
    foreach_scalar_in (gas_int) {
      T[] = 0.;
      foreach_species_in (gas_int) {
        Y[] = 0.;
        X[] = 0.;
      }
    }
  }
}

/**
## Phase change */

event phasechange (i++) {
  phase_tracers_to_scalars (liq, f, tol = F_ERR);
  phase_tracers_to_scalars (gas, f, tol = F_ERR);

  scalar fl[], fg[];
  foreach() {
    fl[] = f[];
    fg[] = 1. - f[];
  }

  face vector fsl[], fsg[];
  face_fraction (fl, fsl);
  face_fraction (fg, fsg);

  /**
  The interface temperature is first estimated from the liquid temperature
  (or set to *TIntVal* if `pcm.isothermal_interface`), using the
  volume-fraction-weighted average over the 3x3 stencil. */

  scalar TLInt = liq_int->T, TGInt = gas_int->T;
  foreach_interfacial (f, F_ERR) {
    TLInt[] = TIntVal;
    TGInt[] = TIntVal;
    if (!pcm.isothermal_interface) {
      TLInt[] = avg_neighbor (point, liq->T, f);
      TGInt[] = TLInt[];
    }
  }

  if (liquid_carrier)
    species_jump_condition (liq, gas, liq_int, gas_int, fg, fsg, 1.);
  else
    species_jump_condition (gas, liq, gas_int, liq_int, fl, fsl, -1.);

  /**
  Unless the interface is isothermal, the interface temperature is refined
  by solving the energy balance ([Cipriano et al., 2023](#cipriano2023),
  Eq. 51)
  $$
  \sum_{i=1}^{NCS} \dot{m}_i\Delta h_{ev,i} =
    \lambda_l\nabla T_l\cdot\mathbf{n}_\Gamma
    - \lambda_g\nabla T_g\cdot\mathbf{n}_\Gamma
  $$
  with $\mathbf{n}_\Gamma$ pointing towards the liquid, the vaporization
  rates frozen, and $\Delta h_{ev,i}$ taken from the carrier species. In
  *energy_balance()*, each gradient is computed along the normal entering its
  own phase, hence the plus sign. */

  if (!pcm.isothermal && !pcm.isothermal_interface) {
    EnergyBalanceData data;
    data.liq = liq, data.gas = gas;
    data.carrier = liquid_carrier ? liq : gas;
    data.fl = fl, data.fg = fg;
    data.fsl = fsl, data.fsg = fsg;

    foreach_interfacial (f, F_ERR, serial) {
      Array * arrUnk = array_new();
      double vali = TLInt[];
      array_append (arrUnk, &vali, sizeof (double));

      coord o = {x, y, z};
      foreach_dimension()
        data.c.x = o.x;

#if USE_GSL
      fsolve (energy_balance_gsl, arrUnk, &data);
#else
      fprintf (ferr,
          "src/evaporation.h:%d: error: missing root-finding method\n",
          LINENO), fflush (ferr);
      exit(1);
#endif

      double * unk = (double *)arrUnk->p;
      TLInt[] = unk[0];
      TGInt[] = TLInt[];
      array_free (arrUnk);
    }
  }

  /**
  The temperature equations get the interface conductive fluxes ([Cipriano
  et al., 2023](#cipriano2023), Eq. 44). */

  foreach_interfacial_plic (f, F_ERR) {
    scalar TL = liq->T, TG = gas->T;
    double ltrgrad = plic_gradient (point, TL, fl, fsl, TLInt[], true, NULL);
    double gtrgrad = plic_gradient (point, TG, fg, fsg, TGInt[], true, NULL);

    scalar slT = liq->STexp, sgT = gas->STexp;
    scalar lambdal = liq->lambda, lambdag = gas->lambda;

    slT[] += lambdal[]*ltrgrad*dirac;
    sgT[] += lambdag[]*gtrgrad*dirac;
  }

  // We add the heat from the mass diffusion process
  if (pcm.mass_diffusion_enthalpy) {
    phase_add_heat_species_diffusion (liq, f, pcm.molar_diffusion, F_ERR);
    phase_add_heat_species_diffusion (gas, f, pcm.molar_diffusion, F_ERR);
  }

  phase_scalars_to_tracers (liq, f);
  phase_scalars_to_tracers (gas, f);
}

/**
## Phase change model parameters

We set the `pcm` parameters for evaporation: the expansion term is shifted to
the liquid, and the species source terms are split into explicit and implicit
contributions. With a single velocity field (`nv == 1`), the interface
regression velocity uses the gas density. */

// fixme: check compatibility with restore
event init (i = 0) {
  pcm.shifting = SHIFT_TO_LIQUID;
  pcm.diffusion = EXPLICIT_IMPLICIT;
  if (nv == 1)
    pcm.byrhogas = true;
}

/**
## References

~~~bib
@hal{cipriano2023, hal-04252476}

@hal{cipriano2024, hal-04580590}
~~~
*/
