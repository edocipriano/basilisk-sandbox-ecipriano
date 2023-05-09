/**
# Multicomponent Phase Change Model

This phase change model describes the evaporation of non-isothermal
systems with multiple chemical species. The interface jump condition
is computed coupling the species and temperature equation in every
interfacial cell, assuming that thermodynamic equilibrium conditions
occur at the gas-liquid interface. The following balance is solved:

$$
  \begin{cases}
    \dot{m}_i = 
    \dot{m}\hat \omega_{i,l} - \rho_l \mathcal{D}_{i,l} \left.\dfrac{\partial
    \omega_{i,l}}{\partial\mathbf{n}_\Gamma}\right\vert_l 
    =
    \dot{m}\hat\omega_{i,g} - \rho_g \mathcal{D}_{i,g} \left.\dfrac{\partial
    \omega_{i,g}}{\partial\mathbf{n}_\Gamma}\right\vert_g
    \\
    \sum_{j=0}^{NS} \dot{m}_i \Delta h_{ev,i} = \lambda_l \left.\dfrac{
    \partial T_l}{\partial \mathbf{n}_\Gamma}\right\vert_l + \lambda_g
    \left.\dfrac{\partial T_g}{\partial \mathbf{n}_\Gamma}\right\vert_g
    \\
    \hat x_{i,g} = k_{eq,i}(\hat T) \hat x_{i,l}
    \\
  \end{cases}
$$

obtaining the vaporization rate $\dot{m}_i$ for each chemical species,
the liquid interface mass fractions $\hat{\omega}_{i,l}$, the gas
interface mass fractions $\hat{\omega}_{i,g}$, and the interface temperature
$\hat{T}$.
model
*/

#define PCM_MULTICOMPONENT

/**
We include the definitions of the interface gradients for the jump condition,
the face fraction for the diffusion equation, and thermodynamic functions
for the equilibrium at the gas-liquid interface.
*/

#include "intgrad.h"
#include "fracface.h"
#include "diffusion.h"
#include "thermodynamics.h"

/**
## Model Options

The following compiler directives can be defined to modify some options
of this code:

* *NGS*: number of chemical species in gas-phase
* *NLS*: number of chemical species in liquid-phase
* *SOLVE_TEMPERATURE*: non-isothermal system: solve the temperature equation
* *USE_CLAPEYRON*: Clausius-Clapeyron equation for the thermodynamic equilibrium
* *USE_ANTOINE*: Antoine equation for the thermodynamic equilibrium (see [thermodynamics.h](thermodynamics.h))

If *USE_SUNDIALS* is defined, we can solve the jump condition as a non-linear
system of equations. The calculation of the interface temperature from the
vaporization rate also relies of this keywork:

* *undefined*: temperature not solved and decoupled solution of jump conditions
* *USE_SUNDIALS = 0*: decoupled solution of jump condition for the chemical species.
The interface temperature is obtained from the non-linear algebraic equations solver.
Required when *SOLVE_TEMPERATURE* is defined.
* *USE_SUNDIALS = 1*: coupdecoupled solution of jump condition and interface temperature
to obtain fist guess values which are then refined by the non-linear algebraic equations
solver
*/

#ifdef USE_SUNDIALS
# include "ijc-coupled.h"
#endif

/**
## User Data

The following data should be provided by the user:

* *gas_species*: list of gas phase species
* *liq_species*: list of liquid phase species
* *inert_species*: name of the inert species
* *gas_start*: initial mass fractions of gas species
* *liq_start*: initial mass fractions of liquid species
* *inDmix1*: diffusivities of species in liquid phase
* *inDmix2*: diffusivities of species in gas phase
* *inKeq*: thermodynamic equilibrium constant
* *Tboil*: boiling temperature of liquid species
(only if *USE_CLAPEYRON*)
*/

extern char* gas_species[NGS];
extern char* liq_species[NLS];
extern char* inert_species[1];

extern double gas_start[NGS];
extern double liq_start[NLS];

extern double inDmix1[NLS];
extern double inDmix2[NGS];

extern double inKeq[NLS];
#ifdef USE_CLAPEYRON
extern double Tboil[NLS];
#endif

/**
By default, the molecular weights are set to 1 for each
species. This can be modified in an init event. The default
thermodynamic reference pressure of the system is 1atm.
*/

double inMW[NGS];
double Pref = 101325.;

/**
If *SOLVE_TEMPERATURE* is defined, additional data must be
defined:

* *lambda1*: thermal conductivity in liquid phase
* *lambda2*: thermal conductivity in gas phase
* *dhev*: enthalpy of evaporation
* *cp1*: specific heat for the liquid phase
* *cp2*: specific heat for the gas phase
* *TG0, TL0*: inital temperatures in gas and liquid phase
*/

#ifdef SOLVE_TEMPERATURE
extern double lambda1, lambda2, dhev, cp1, cp2;
extern double TL0, TG0;

//scalar T[], TL[], TG[], TInt[];
scalar T[], TInt[];
scalar TL, TG;
scalar sgT[], slT[], sgTimp[], slTimp[];
face vector lambda1f[], lambda2f[];
#endif

/**
## Fields

Since we don't know a-priori how many species
we want to simulate, we initialize empty lists
containing the vaporization rate for every species,
the mass fractions, the interface properties and
phase change source terms for the transport equations.
*/

scalar * mEvapList = NULL;    // [NGS]
scalar * YList     = NULL;    // [NGS]
scalar * YLList    = NULL;    // [NLS]
scalar * YGList    = NULL;    // [NGS]
scalar * YLIntList = NULL;    // [NLS]
scalar * YGIntList = NULL;    // [NGS]
scalar * slexpList = NULL;    // [NLS]
scalar * slimpList = NULL;    // [NLS]
scalar * sgexpList = NULL;    // [NGS]
scalar * sgimpList = NULL;    // [NGS]

/**
We declare useful fields used for loops over chemical
species:

* *LSI*: vector with Liquid Species Indices [NLS]
* *GOSI*: vector with Gas-Only Species Indices [NGOS]
* *NGOS*: number of Gas-Only Species Indices
* *inertIndex*: index/position of the inert species.
*/

int * LSI = NULL;
int * GOSI = NULL;
int NGOS;
int inertIndex;

/**
We initilize other useful fields. */

bool success;
bool init_fields;

scalar fG[], fL[], fuT[];
face vector fsL[], fsG[];
scalar f0[];

scalar divu[], fold[];

/**
## Defaults

In the defaults event, the lists containing the species
are initialized from the lists of species names in gas
and liquid phase. The tracers are assigned in this event.
*/

event defaults (i = 0)
{
  /**
  Fill Liquid Species Indices *LSI* array. */

  Array * arrLSI = array_new();
  for (int ii=0; ii<NGS; ii++) {
    for (int jj=0; jj<NLS; jj++) {
      if (strcmp(gas_species[ii], liq_species[jj]) == 0) {
        int idx = ii; array_append (arrLSI, &idx, sizeof(int));
      }
    }
  }
  LSI = (int *) array_shrink (arrLSI);

  /**
  Fill Gas Only Species Indices *GOSI* array. */

  Array * arrGOSI = array_new();
  for (int ii=0; ii<NGS; ii++) {
    bool thisSpeciesIsAlsoLiquid = false;
    for (int jj=0; jj<NLS; jj++) {
      if (strcmp(liq_species[jj], gas_species[ii]) == 0) {
        thisSpeciesIsAlsoLiquid = true;
      }
    }
    if (!thisSpeciesIsAlsoLiquid) {
      int idx = ii; array_append (arrGOSI, &idx, sizeof(int));
    }
  }
  NGOS = arrGOSI->len/sizeof(int);
  GOSI = (int *) array_shrink (arrGOSI);
  assert (NGOS == (NGS - NLS));

  /**
  Fill index of inert species among the gas species. */

  inertIndex = -1;
  for (int ii=0; ii<NGS; ii++) {
    if (strcmp(gas_species[ii], inert_species[0]) == 0)
      inertIndex = ii;
  }

  /**
  We create species fields from the list of species names. */

  YList     = NULL;
  YLList    = NULL;
  YGList    = NULL;
  YGIntList = NULL;
  YLIntList = NULL;
  mEvapList = NULL;
  slexpList = NULL;
  slimpList = NULL;
  sgexpList = NULL;
  sgimpList = NULL;

  for (int jj=0; jj<NLS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "%s_L", liq_species[jj]);
    a.name = strdup (name);
    YLList = list_append (YLList, a);
  }
  reset (YLList, 0.);

  for (int jj=0; jj<NGS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "%s_G", gas_species[jj]);
    a.name = strdup (name);
    YGList = list_append (YGList, a);
  }
  reset (YGList, 0.);

  for (int jj=0; jj<NLS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "%s_LInt", liq_species[jj]);
    a.name = strdup (name);
    YLIntList = list_append (YLIntList, a);
  }
  reset (YLIntList, 0.);

  for (int jj=0; jj<NGS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "%s_GInt", gas_species[jj]);
    a.name = strdup (name);
    YGIntList = list_append (YGIntList, a);
  }
  reset (YGIntList, 0.);

  for (int jj=0; jj<NGS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "%s", gas_species[jj]);
    a.name = strdup (name);
    YList = list_append (YList, a);
  }
  reset (YList, 0.);

  for (int jj=0; jj<NGS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "mEvap_%s", gas_species[jj]);
    a.name = strdup (name);
    mEvapList = list_append (mEvapList, a);
  }
  reset (mEvapList, 0.);

  for (int jj=0; jj<NLS; jj++) {
    scalar a = new scalar;
    scalar b = new scalar;
    free (a.name);
    free (b.name);
    char aname[20];
    char bname[20];
    sprintf (aname, "slexp_%s", liq_species[jj]);
    sprintf (bname, "slimp_%s", liq_species[jj]);
    a.name = strdup (aname);
    b.name = strdup (bname);
    a.nodump = true;
    b.nodump = true;
    slexpList = list_append (slexpList, a);
    slimpList = list_append (slimpList, b);
  }
  reset (slexpList, 0.);
  reset (slimpList, 0.);

  for (int jj=0; jj<NGS; jj++) {
    scalar a = new scalar;
    scalar b = new scalar;
    free (a.name);
    free (b.name);
    char aname[20];
    char bname[20];
    sprintf (aname, "sgexp_%s", gas_species[jj]);
    sprintf (bname, "sgimp_%s", gas_species[jj]);
    a.name = strdup (aname);
    b.name = strdup (bname);
    a.nodump = true;
    b.nodump = true;
    sgexpList = list_append (sgexpList, a);
    sgimpList = list_append (sgimpList, b);
  }
  reset (sgexpList, 0.);
  reset (sgimpList, 0.);

  fL.nodump = true;
  fG.nodump = true;

  for (scalar s in YGList)
    s.inverse = true;

  for (scalar s in YLList)
    s.inverse = false;

  /**
  We restore to NULL the tracers associated with fu
  and fuext in order to re-run simulations correctly. */

  fuext.tracers = NULL;
  fu.tracers    = NULL;

#ifdef CONSISTENTPHASE1
  fuext.tracers = list_concat (fuext.tracers, YLList);
#else
  fu.tracers = list_concat (fu.tracers, YLList);
#endif
#ifdef CONSISTENTPHASE2
  fuext.tracers = list_concat (fuext.tracers, YGList);
#else
  fu.tracers = list_concat (fu.tracers, YGList);
#endif

  /**
  On adaptive meshes, tracers need to use linear interpolation (rather
  than the default bilinear interpolation) to ensure conservation when
  refining cells. */

#if TREE
    for (scalar s in YGList) {
#if EMBED
      s.refine = s.prolongation = refine_embed_linear;
#else
      s.refine  = refine_linear;
#endif
      s.restriction = restriction_volume_average;
      s.dirty = true; // boundary conditions need to be updated
    }
#endif

#if TREE
    for (scalar s in YLList) {
#if EMBED
      s.refine = s.prolongation = refine_embed_linear;
#else
      s.refine  = refine_linear;
#endif
      s.restriction = restriction_volume_average;
      s.dirty = true; // boundary conditions need to be updated
    }
#endif

#ifdef SOLVE_TEMPERATURE

  TL = new scalar;
  TG = new scalar;

  TL.inverse = false;
  TG.inverse = true;

  sgT.nodump    = true;
  slT.nodump    = true;
  sgTimp.nodump = true;
  slTimp.nodump = true;

#ifdef CONSISTENTPHASE1
  fuext.tracers = list_append (fuext.tracers, TL);
#else
  fu.tracers = list_append (fu.tracers, TL);
#endif
#ifdef CONSISTENTPHASE2
  fuext.tracers = list_append (fuext.tracers, TG);
#else
  fu.tracers = list_append (fu.tracers, TG);
#endif
  TL.refine  = refine_linear;
  TL.restriction = restriction_volume_average;
  TL.dirty = true; // boundary conditions need to be updated
  TG.refine  = refine_linear;
  TG.restriction = restriction_volume_average;
  TG.dirty = true; // boundary conditions need to be updated
#endif

  /**
  Set default multicomponent properties. */

  for (int jj=0; jj<NGS; jj++)
    inMW[jj] = 1.;

  init_fields = true;
}

/**
## Init

In the init event, We initialize the temperature and
chemical species and set to zero additional fields. */

event init (i = 0)
{
  foreach() {
    for (int jj=0; jj<NLS; jj++) {
      scalar s = YLList[jj];
      if (init_fields)
        s[] = f[]*liq_start[jj];
      else
        s[] = 0.;
    }
    for (int jj=0; jj<NGS; jj++) {
      scalar s = YGList[jj];
      if (init_fields)
        s[] = (1. - f[])*gas_start[jj];
      else
        s[] = 0.;
    }
    for (int jj=0; jj<NLS; jj++) {
      scalar s = YLIntList[jj];
      s[] = 0.;
    }
    for (int jj=0; jj<NGS; jj++) {
      scalar s = YGIntList[jj];
      s[] = 0.;
    }
    for (int jj=0; jj<NGS; jj++) {
      scalar s = mEvapList[jj];
      s[] = 0.;
    }
    for (int jj=0; jj<NGOS; jj++) {
      scalar s  = YList[GOSI[jj]];
      scalar sg = YGList[GOSI[jj]];
      s[] = sg[];
    }
    for (int jj=0; jj<NLS; jj++) {
      scalar s  = YList[LSI[jj]];
      scalar sg = YGList[LSI[jj]];
      scalar sl = YLList[jj];
      s[] = sg[] + sl[];
    }
    for (int jj=0; jj<NLS; jj++) {
      scalar slexp = slexpList[jj];
      scalar slimp = slimpList[jj];
      slexp[] = 0.;
      slimp[] = 0.;
    }
    for (int jj=0; jj<NGS; jj++) {
      scalar sgexp = sgexpList[jj];
      scalar sgimp = sgimpList[jj];
      sgexp[] = 0.;
      sgimp[] = 0.;
    }
    for (int jj=0; jj<NGS; jj++) {
      scalar mEvap = mEvapList[jj];
      mEvap[] = 0.;
    }
  }
  boundary(YList);
  boundary(YLList);
  boundary(YGList);
  boundary(YGIntList);
  boundary(YLIntList);
  boundary(mEvapList);
  boundary(slexpList);
  boundary(slimpList);
  boundary(sgexpList);
  boundary(sgimpList);

#ifdef SOLVE_TEMPERATURE
  foreach() {
    TL[] = TL0*f[];
    TG[] = TG0*(1. - f[]);
    T[]  = TL[] + TG[];
  }
  boundary({T,TL,TG});
#endif
  
}

/**
## Cleanup

We deallocate the various lists from the memory. */

event cleanup (t = end)
{
  delete (mEvapList), free (mEvapList), mEvapList = NULL;
  delete (YList), free (YList), YList = NULL;
  delete (YLList), free (YLList), YLList = NULL;
  delete (YGList), free (YGList), YGList = NULL;
  delete (YLIntList), free (YLIntList), YLIntList = NULL;
  delete (YGIntList), free (YGIntList), YGIntList = NULL;
  delete (slexpList), free (slexpList), slexpList = NULL;
  delete (slimpList), free (slimpList), slimpList = NULL;
  delete (sgexpList), free (sgexpList), sgexpList = NULL;
  delete (sgimpList), free (sgimpList), sgimpList = NULL;
  delete (fu.tracers), free (fu.tracers), fu.tracers = NULL;
  delete (fuext.tracers), free (fuext.tracers), fuext.tracers = NULL;
  free (LSI);
  free (GOSI);
#ifdef SOLVE_TEMPERATURE
  delete ({TL,TG});
#endif
}

/**
## Phase Change

In the *phasechange* event, the vaporization rate is computed
and the diffusion step for the mass fraction field (in liquid
and gas phase) is solved. */

event phasechange (i++)
{
  /**
  First we compute the value of the non volume-averaged
  temperature fields. This procedure allows a better
  calculation of the gradients close to the interface. */

  foreach() {
    f[] = clamp (f[], 0., 1.);
    f[] = (f[] > F_ERR) ? f[] : 0.;
    f0[] = f[];
    fL[] = f[]; fG[] = 1. - f[];

#ifdef SOLVE_TEMPERATURE
    TL[] = f[] > F_ERR ? TL[]/f[] : 0.;
    TG[] = ((1. - f[]) > F_ERR) ? TG[]/(1. - f[]) : 0.;
#endif

    /**
    We recover value of mass fractions from their tracer form. */

    for (int jj=0; jj<NLS; jj++) {
      scalar YL = YLList[jj];
      YL[] = f[] > F_ERR ? YL[]/f[] : 0.;
    }

    for (int jj=0; jj<NGS; jj++) {
      scalar YG = YGList[jj];
        YG[] = ((1. - f[]) > F_ERR) ? YG[]/(1. - f[]) : 0.;
    }
  }
  boundary({f,fL,fG});
  boundary(YGList);
  boundary(YLList);

  /**
  We compute the value of volume fraction *f* on the
  cell-faces using a geometric approach (necessary
  for interface gradients and diffusion equations). */

  face_fraction (fL, fsL);
  face_fraction (fG, fsG);

  /**
  We assign the interface temperature value. */

#ifdef SOLVE_TEMPERATURE
  foreach() {
    TInt[] = 0.;
    if (f[] > F_ERR && f[] < 1.-F_ERR)
      TInt[] = avg_neighbor (point, TL, f);
  }
  boundary({TInt});
#endif

  /**
  We compute the molecular weight of the gas-only
  species mixture. */

  scalar MWGmix[];
  foreach() {
    MWGmix[] = 0.;
    if (f[] > F_ERR && f[] < 1.-F_ERR) {
      double yGinorm[NGOS];
      double sumYGi = 0.;
      for (int jj=0; jj<NGOS; jj++) {
        scalar YG = YGList[GOSI[jj]];
        yGinorm[jj] = YG[];
        sumYGi += yGinorm[jj];
      }
      for (int jj=0; jj<NGOS; jj++) {
        yGinorm[jj] /= (sumYGi + 1.e-10);
      }
      for (int jj=0; jj<NGOS; jj++) {
        MWGmix[] += yGinorm[jj] / inMW[GOSI[jj]];
      }
      MWGmix[] = 1. / (MWGmix[] + 1.e-10);
    }
  }
  boundary({MWGmix});

  /**
  We compute total vaporization flowrate. */

  scalar sumYGInt[];
  foreach() {
    mEvapTot[] = 0.;
    double sum_jG = 0.;
    double sum_YGInt = 0.;
    sumYGInt[] = 0.;

    /**
    We reset to zero mEvap for every species, and we set to zero
    the interface mass fraction fields. */

    for (int jj=0; jj<NGS; jj++) {
      scalar mEvap = mEvapList[jj];
      scalar YGInt = YGIntList[jj];
      mEvap[] = 0.;
      YGInt[] = 0.;
    }
    for (int jj=0; jj<NLS; jj++) {
      scalar YLInt = YLIntList[jj];
      YLInt[] = 0.;
    }

    if (f[] > F_ERR && f[] < 1.-F_ERR) {

      /**
      We create fields to store the local mole fractions, the
      conversion is required by the thermodynamic equilibrium. */

      double XGIntConv[NLS+1], YGIntConv[NLS+1];
      double XLIntConv[NLS], YLIntConv[NLS];
      double inMWG[NLS+1], inMWL[NLS];

      for (int jj=0; jj<NLS; jj++) {
        inMWL[jj] = inMW[LSI[jj]];
        inMWG[jj] = inMW[LSI[jj]];
        XGIntConv[jj] = 0.;
        YGIntConv[jj] = 0.;
        XLIntConv[jj] = 0.;
        YLIntConv[jj] = 0.;
      }

      /**
      We convert *YLInt* to mole fractions *XLInt*. */

      for (int jj=0; jj<NLS; jj++) {
        scalar YL = YLList[jj];
        YLIntConv[jj] = avg_neighbor (point, YL, f);
      }
      mass2molefrac (XLIntConv, YLIntConv, inMWL, NLS);

      /**
      We compute *XGInt* from the thermodynamic equilibrium.
      Different equilibrium options are available: constant
      thermodynamic equilibrium constant, Clausius-Clapeyron
      relation, and Antoine equation. */

      double sumXGi = 0.;
      for (int jj=0; jj<NLS; jj++) {
        XGIntConv[jj] = inKeq[jj]*XLIntConv[jj];
#ifdef USE_CLAPEYRON
        XGIntConv[jj] = clapeyron (min (TInt[], Tboil[jj]-1.), Tboil[jj], dhev, inMW[jj])*XLIntConv[jj];
#endif
#ifdef USE_ANTOINE
        scalar YL = YLList[jj];
        XGIntConv[jj] = min (YL.antoine (TInt[], Pref), 0.98)*XLIntConv[jj];
#endif
        sumXGi += XGIntConv[jj];
      }
      XGIntConv[NLS] = 1. - sumXGi;
      inMWG[NLS] = MWGmix[];
      mole2massfrac (YGIntConv, XGIntConv, inMWG, NLS+1);

      /**
      We compute the total vaporization rate from the sum
      of the interface mass balances over all the chemical
      species in liquid phase. */

      for (int jj=0; jj<NLS; jj++) {
        scalar YGInt = YGIntList[LSI[jj]];
        scalar YLInt = YLIntList[jj];
        scalar YL    = YLList[jj];
        scalar YG    = YGList[LSI[jj]];

        YLInt[] = avg_neighbor (point, YL, f);
        YGInt[] = YGIntConv[jj];

        double gtrgrad = ebmgrad (point, YG, fL, fG, fsL, fsG, true, YGInt[], &success);
        sum_jG += -rho2*inDmix2[LSI[jj]]*gtrgrad;
        sum_YGInt += YGInt[];
        sumYGInt[] += YGInt[];
      }
#ifdef DIFFUSIVE
      mEvapTot[] = sum_jG;
#else
      mEvapTot[] = sum_jG/min(1. - sum_YGInt, 0.99);
#endif
    }
  }

  /**
  From the knowledge of the total vaporization rate, we
  compute the vaporization rate for each species. */

  foreach() {
    if (f[] > F_ERR && f[] < 1.-F_ERR) {

      for (int jj=0; jj<NLS; jj++) {
        scalar mEvap = mEvapList[LSI[jj]];
        scalar YG    = YGList[LSI[jj]];
        scalar YGInt = YGIntList[LSI[jj]];

        double gtrgrad = ebmgrad (point, YG, fL, fG, fsL, fsG, true, YGInt[], &success);
#ifdef DIFFUSIVE
        mEvap[] = - rho2*inDmix2[LSI[jj]]*gtrgrad;
#else
        mEvap[] = mEvapTot[]*YGInt[] - rho2*inDmix2[LSI[jj]]*gtrgrad;
#endif
      }
    }

    /**
    We adjust the interface mass fractions of the gas-only
    species in the system. */

    double yGinorm[NGOS];
    double sumYGi = 0.;
    if (f[] > F_ERR && f[] < 1.-F_ERR) {
      for (int jj=0; jj<NGOS; jj++) {
        scalar YG = YGList[GOSI[jj]];
        yGinorm[jj] = YG[];
        sumYGi += yGinorm[jj];
      }
    }

    for (int jj=0; jj<NGOS; jj++) {
      scalar mEvap = mEvapList[GOSI[jj]];
      scalar YGInt = YGIntList[GOSI[jj]];
      //scalar YG    = YGList[GOSI[jj]];

      mEvap[] = 0.; YGInt[] = 0.; YGInt[] = 0.;

      if (f[] > F_ERR && f[] < 1.-F_ERR) {
        YGInt[] = (1. - sumYGInt[])*yGinorm[jj]/(sumYGi + 1.e-10);
      }
    }
  }

  /**
  If *SOLVE_TEMPERATURE*, the value of interface temperature
  is obtained from the total vaporization rate just computed
  and the interface temperature gradients. This result in an
  algebraic equation which is numerically zero-ed. */

#ifdef USE_SUNDIALS

#ifdef SOLVE_TEMPERATURE
  ijc_CoupledTemperature();
#endif

  /**
  Finally, we can refine the first guess values from the decoupled
  solution from the fully-coupled interface jump condition. Solved
  as a non-linear system of equations. */

#if USE_SUNDIALS > 0
  ijc_CoupledNls();
#endif

#endif

  /**
  The source terms for the diffusion equation of the species
  mass fractions in gas an liquid phase are computed here. */

  foreach() {

    for (int jj=0; jj<NLS; jj++) {
      scalar slexp = slexpList[jj];
      scalar slimp = slimpList[jj];

      slexp[] = 0.;
      slimp[] = 0.;
    }
    for (int jj=0; jj<NGS; jj++) {
      scalar sgexp = sgexpList[jj];
      scalar sgimp = sgimpList[jj];

      sgexp[] = 0.;
      sgimp[] = 0.;
    }

    if (f[] > F_ERR && f[] < 1.-F_ERR) {
      coord n = facet_normal (point, fL, fsL), p;
      double alpha = plane_alpha (fL[], n);
      double area = plane_area_center (n, alpha, &p);
      normalize (&n);

      for (int jj=0; jj<NGS; jj++) {
        scalar sgexp = sgexpList[jj];
        scalar sgimp = sgimpList[jj];
        scalar mEvap = mEvapList[jj];

#ifdef AXI
        sgexp[] = -mEvap[]/rho2*area*(y + p.y*Delta)/(Delta*y)*cm[];
        sgimp[] = +mEvapTot[]/rho2*area*(y + p.y*Delta)/(Delta*y)*cm[];
#else
        sgexp[] = -mEvap[]/rho2*area/Delta*cm[];
        sgimp[] = +mEvapTot[]/rho2*area/Delta*cm[];
#endif
      }

      for (int jj=0; jj<NLS; jj++) {
        scalar slexp = slexpList[jj];
        scalar slimp = slimpList[jj];
        scalar mEvap = mEvapList[LSI[jj]];

#ifdef AXI
        slexp[] = +mEvap[]/rho1*area*(y + p.y*Delta)/(Delta*y)*cm[];
        slimp[] = -mEvapTot[]/rho1*area*(y + p.y*Delta)/(Delta*y)*cm[];
#else
        slexp[] = +mEvap[]/rho1*area/Delta*cm[];
        slimp[] = -mEvapTot[]/rho1*area/Delta*cm[];
#endif
      }
    }
  }
  boundary(slexpList);
  boundary(sgexpList);
  boundary(slimpList);
  boundary(sgimpList);

  /**
  We restore the tracer form of the liquid and gas-phase
  mass fraction fields. */

  foreach() {
    for (scalar YL in YLList)
      YL[] *= f[]*(f[] > F_ERR);
    for (scalar YG in YGList) 
      YG[] *= (1. - f[])*((1. - f[]) > F_ERR);
#ifdef SOLVE_TEMPERATURE
    TL[] *= f[]*(f[] > F_ERR);
    TG[] *= (1. - f[])*((1. - f[]) > F_ERR);
#endif
  }
  boundary(YLList);
  boundary(YGList);
#ifdef SOLVE_TEMPERATURE
  boundary({TL,TG});
#endif
}

/**
## Tracer Advection

We let the volume fractions *fu* and *fuext* to
advect the fields YL and YG, as implemented in
the tracer_advection event of [evaporation.h](evaporation.h)
*/

event tracer_advection (i++);

/**
## Tracer Diffusion

We solve the diffusion equations for species and temperature
accounting for the phase change contributions. */

event tracer_diffusion (i++)
{
  /**
  We remove the fractions of f and mass fractions
  lower than F_ERR and we reconstruct the non-volume
  averaged form of the mass fraction fields, in order
  to improve the discretization of the face gradients
  in the diffusion equation. */

  foreach() {
    f[] = clamp (f[], 0., 1.);
    f[] = (f[] > F_ERR) ? f[] : 0.;

    for (scalar YL in YLList)
      YL[] = (fuext[] > F_ERR) ? YL[]/fuext[] : 0.;
    for (scalar YG in YGList)
      YG[] = ((1. - fu[]) > F_ERR) ? YG[]/(1. - fu[]) : 0.;

    fL[] = f[]; fG[] = 1. - f[];

#ifdef SOLVE_TEMPERATURE
    TL[] = (fuext[] > F_ERR) ? TL[]/fuext[] : 0.;
    TG[] = ((1. - fu[]) > F_ERR) ? TG[]/(1. - fu[]) : 0.;
#endif
  }
  boundary({fL,fG});
  boundary(YLList);
  boundary(YGList);
#ifdef SOLVE_TEMPERATURE
  boundary({TL,TG});
#endif

  /**
  We compute the value of volume fraction *f* on the
  cell-faces using a geometric approach (necessary
  for interface gradients and diffusion equations). */

  face_fraction (fL, fsL);
  face_fraction (fG, fsG);

  /**
  We solve the diffusion equations, confined by means of
  the face fraction fields *fsL* and *fsG*. */

  scalar theta1[], theta2[];

#if TREE
  theta1.refine = theta1.prolongation = fraction_refine;
  theta2.refine = theta2.prolongation = fraction_refine;
  //theta1.dirty = true;
  //theta2.dirty = true;
#endif

#ifdef MODIFIED_DIFFUSION
  vector pcs1[], pcs2[];

  foreach() {
    if (fL[] > F_ERR && fL[] < 1.-F_ERR) {
      coord m = mycs (point, fL); 
      double alpha = plane_alpha (fL[], m); 
      coord prel;
      plane_area_center (m, alpha, &prel);
      coord pc; 
      plane_center (m, alpha, fL[], &pc);

      pcs1.x[] = x + pc.x*Delta;
      pcs1.y[] = y + pc.y*Delta;
      pcs1.z[] = z + pc.z*Delta;
    }
    else {
      pcs1.x[] = x;
      pcs1.y[] = y;
      pcs1.z[] = z;
    }

    if (fG[] > F_ERR && fG[] < 1.-F_ERR) {
      coord m = mycs (point, fG); 
      double alpha = plane_alpha (fG[], m); 
      coord prel;
      plane_area_center (m, alpha, &prel);
      coord pc; 
      plane_center (m, alpha, fG[], &pc);

      pcs2.x[] = x + pc.x*Delta;
      pcs2.y[] = y + pc.y*Delta;
      pcs2.z[] = z + pc.z*Delta;
    }
    else {
      pcs2.x[] = x;
      pcs2.y[] = y;
      pcs2.z[] = z;
    }
  }

  face vector corrdist1[], corrdist2[];
  foreach_face() {
    corrdist1.x[] = Delta/fabs (pcs1.x[] - pcs1.x[-1]);
    corrdist2.x[] = Delta/fabs (pcs2.x[] - pcs2.x[-1]);

    // check: Palmore version, I don't think it's correct

    //double magdist1 = 0., magdist2 = 0.;
    //foreach_dimension() {
    //  magdist1 += sq (pcs1.x[] - pcs1.x[-1]);
    //  magdist2 += sq (pcs2.x[] - pcs2.x[-1]);
    //}
    //magdist1 = sqrt (magdist1);
    //magdist2 = sqrt (magdist2);

    //corrdist1.x[] = Delta*(pcs1.x[] - pcs1.x[-1])/magdist1;
    //corrdist2.x[] = Delta*(pcs2.x[] - pcs2.x[-1])/magdist2;
  }
  boundary({corrdist2});
#endif

  for (int jj=0; jj<NLS; jj++) {
    face vector Dmix1f[];
    foreach_face() {
      Dmix1f.x[] = inDmix1[jj]*fsL.x[]*fm.x[];
#ifdef MODIFIED_DIFFUSION
      Dmix1f.x[] *= corrdist1.x[];
#endif
    }
    boundary((scalar *){Dmix1f});

    foreach()
      theta1[] = cm[]*max(fL[], 1.e-3);
      //theta1[] = cm[]*max(fL[], T_ERR);
    boundary({theta1});

    scalar YL = YLList[jj];
    scalar slexp = slexpList[jj];
#ifndef USE_DALPHADT
    scalar slimp = slimpList[jj];
#else
    scalar slimp[];
    foreach() {
      slimp[] = -(f[] - f0[])/dt;
#ifdef AXI
      slimp[] *= y;
#endif
    }
    boundary({slimp});
#endif

    diffusion (YL, dt, D=Dmix1f, r=slexp, beta=slimp, theta=theta1);
  }

  for (int jj=0; jj<NGS; jj++) {

    face vector Dmix2f[];
    foreach_face() {
      Dmix2f.x[] = inDmix2[jj]*fsG.x[]*fm.x[];
#ifdef MODIFIED_DIFFUSION
      Dmix2f.x[] *= corrdist2.x[];
#endif
    }
    boundary((scalar *){Dmix2f});

    foreach()
      theta2[] = cm[]*max(fG[], F_ERR);
    boundary({theta2});

    scalar YG = YGList[jj];
    scalar sgexp = sgexpList[jj];
#ifndef USE_DALPHADT
    scalar sgimp = sgimpList[jj];
#else
    scalar sgimp[];
    foreach() {
      sgimp[] = -(f0[] - f[])/dt;
#ifdef AXI
      sgimp[] *= y;
#endif
    }
    boundary({sgimp});
#endif

    diffusion (YG, dt, D=Dmix2f, r=sgexp, beta=sgimp, theta=theta2);
  }

#ifdef SOLVE_TEMPERATURE

  foreach_face() {
    lambda1f.x[] = lambda1/rho1/cp1*fsL.x[]*fm.x[];
    lambda2f.x[] = lambda2/rho2/cp2*fsG.x[]*fm.x[];
#ifdef MODIFIED_DIFFUSION
    lambda1f.x[] *= corrdist1.x[];  
    lambda2f.x[] *= corrdist2.x[];
#endif
  }
  boundary((scalar*){lambda1f, lambda2f});

  /**
  Compute source terms for temperature equations. */

  foreach() {
    sgT[] = 0.; slT[] = 0.;
    sgTimp[] = 0.; slTimp[] = 0.;
    theta1[] = cm[]*max(f[], F_ERR);
    theta2[] = cm[]*max(1. - f[], F_ERR);
    if (f0[] > F_ERR && f0[] < 1.-F_ERR) {
      coord n = facet_normal (point, fL, fsL), p;
      double alpha = plane_alpha (fL[], n);
      double area = plane_area_center (n, alpha, &p);
      normalize (&n);

      double bc = TInt[];
      double gtrgrad = ebmgrad (point, TG, fL, fG, fsL, fsG, true, bc, &success);
      double ltrgrad = ebmgrad (point, TL, fL, fG, fsL, fsG, false, bc, &success);

      double lheatflux = lambda1*ltrgrad;
      double gheatflux = lambda2*gtrgrad;

#ifdef AXI
      slT[] = lheatflux/rho1/cp1*area*(y + p.y*Delta)/(Delta*y)*cm[];
      sgT[] = gheatflux/rho2/cp2*area*(y + p.y*Delta)/(Delta*y)*cm[];
      slTimp[] = mEvapTot[]/rho1*area*(y + p.y*Delta)/(Delta*y)*cm[];
      sgTimp[] = mEvapTot[]/rho2*area*(y + p.y*Delta)/(Delta*y)*cm[];
#else
      slT[] = lheatflux/rho1/cp1*area/Delta*cm[];
      sgT[] = gheatflux/rho2/cp2*area/Delta*cm[];
      slTimp[] = mEvapTot[]/rho1*area/Delta*cm[];
      sgTimp[] = mEvapTot[]/rho2*area/Delta*cm[];
#endif
    }
  }
  boundary({slT,sgT});
  boundary({slTimp, sgTimp});

  foreach() {
    theta1[] = cm[]*max(fL[], F_ERR);
    theta2[] = cm[]*max(fG[], F_ERR);
  }
  boundary({theta1,theta2});

  /**
  Solve diffusion equations for temperature. */

  diffusion (TL, dt, D=lambda1f, r=slT, theta=theta1);
  diffusion (TG, dt, D=lambda2f, r=sgT, theta=theta2);

#endif

  foreach() {
    /**
    We reconstruct the mass fractions from the total mass
    in every cell. */

    double totmassliq = 0.;
    for (scalar YL in YLList)
      totmassliq += rho1*YL[]*fL[]*dv();
    for (scalar YL in YLList)
      YL[] = (totmassliq > 0.) ? rho1*YL[]*fL[]*dv() / totmassliq : 0.;

    double totmassgas = 0.;
    for (scalar YG in YGList)
      totmassgas += rho2*YG[]*fG[]*dv();
    for (scalar YG in YGList)
      YG[] = (totmassgas > 0.) ? rho2*YG[]*fG[]*dv() / totmassgas : 0.;

    for (scalar YL in YLList)
      YL[] *= f[];
    for (scalar YG in YGList)
      YG[] *= (1. - f[]);
#ifdef SOLVE_TEMPERATURE
    TL[] *= f[];
    TG[] *= (1. - f[]);
#endif
  }
  boundary(YLList);
  boundary(YGList);
#ifdef SOLVE_TEMPERATURE
  boundary({TL,TG});
#endif

  /**
  We reconstruct the volume-averaged mass fractions and
  the one-field mass fraction field. */
  
  foreach() {
    for (scalar YL in YLList) {
      YL[] = clamp (YL[], 0., 1.);
      YL[] = (YL[] > F_ERR) ? YL[] : 0.;
    }
    for (scalar YG in YGList) {
      YG[] = clamp (YG[], 0., 1.);
      YG[] = (YG[] > F_ERR) ? YG[] : 0.;
    }
    for (int jj=0; jj<NGOS; jj++) {
      scalar Y  = YList[GOSI[jj]];
      scalar YG = YGList[GOSI[jj]];

      Y[] = YG[];
    }
    for (int jj=0; jj<NLS; jj++) {
      scalar Y  = YList[LSI[jj]];
      scalar YG = YGList[LSI[jj]];
      scalar YL = YLList[jj];

      Y[] = YL[] + YG[];
    }
  }
  boundary(YList);
  boundary(YGList);

#ifdef SOLVE_TEMPERATURE
  foreach()
    T[] = TL[] + TG[];
  boundary({T,TL,TG});
#endif
}

