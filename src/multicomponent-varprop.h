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

If *USE_GSL* is defined, we can solve the jump condition as a non-linear
system of equations. The calculation of the interface temperature from the
vaporization rate also relies of this keyword:

* *undefined*: temperature not solved, only decoupled solution of jump conditions
* *USE_GSL = 0*: decoupled solution of jump condition for the chemical species.
The interface temperature is obtained from the non-linear algebraic equations solver.
Required when *SOLVE_TEMPERATURE* is defined.
* *USE_GSL = 1*: decoupled solution of jump condition and interface temperature
to obtain fist guess values which are then refined by the non-linear algebraic equations
solver.
*/

#ifdef USE_GSL
# include "ijc-coupled-gsl.h"
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
#ifdef VARPROP
scalar * Dmix1List = NULL;    // [NLS]
scalar * Dmix2List = NULL;    // [NGS]
scalar * dhevList  = NULL;    // [NLS]
#endif

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
bool update_properties = true;
bool first_iter = true;

scalar fG[], fL[], fuT[];
face vector fsL[], fsG[];
scalar f0[];

scalar divu[], fold[];

/**
Variable properties stuff. */

#ifdef VARPROP
scalar dummy[];
scalar frho1r[], frho2r[];
scalar frho1[], frho2[];
scalar frhocp1r[], frhocp2r[];
scalar frhocp1[], frhocp2[];
scalar frhocp1r0[], frhocp2r0[];
scalar rho1v0[], rho2v0[];
#endif

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
#ifdef VARPROP
  Dmix1List = NULL;
  Dmix2List = NULL;
  dhevList  = NULL;
#endif

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

#ifdef VARPROP
  for (int jj=0; jj<NLS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "Dmix1_%s", liq_species[jj]);
    a.name = strdup (name);
    //a.nodump = true;
    Dmix1List = list_append (Dmix1List, a);
  }
  for (int jj=0; jj<NGS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "Dmix2_%s", gas_species[jj]);
    a.name = strdup (name);
    //a.nodump = true;
    Dmix2List = list_append (Dmix2List, a);
  }
  reset (Dmix1List, 0.);
  reset (Dmix2List, 0.);

  for (int jj=0; jj<NLS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "dhev_%s", liq_species[jj]);
    a.name = strdup (name);
    //a.nodump = true;
    dhevList = list_append (dhevList, a);
  }
  reset (dhevList, 0.);
#endif

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

#ifdef VARPROP
  frho1r.inverse = false;
  frho2r.inverse = true;
  frhocp1r.inverse = false;
  frhocp2r.inverse = true;

  scalar * forig = list_copy (f.tracers);
  free (f.tracers);
  f.tracers = list_concat (forig, {frho1r,frho2r,frhocp1r,frhocp2r});
  free (forig);

  frho1.inverse = false;
  frho2.inverse = true;
  frhocp1.inverse = false;
  frhocp2.inverse = true;

  fuext.tracers = list_concat (fuext.tracers, {frho1});
  fuext.tracers = list_concat (fuext.tracers, {frhocp1});
  fu.tracers = list_concat (fu.tracers, {frho2});
  fu.tracers = list_concat (fu.tracers, {frhocp2});
#endif

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

#if TREE
#if EMBED
  TL.refine = TL.prolongation = refine_embed_linear;
  TG.refine = TG.prolongation = refine_embed_linear;
#else
  TL.refine  = refine_linear;
  TG.refine  = refine_linear;
#endif
  TL.restriction = restriction_volume_average;
  TL.dirty = true; // boundary conditions need to be updated
  TG.restriction = restriction_volume_average;
  TG.dirty = true; // boundary conditions need to be updated
#endif

#endif

  /**
  Set default multicomponent properties. */

  for (int jj=0; jj<NGS; jj++)
    inMW[jj] = 1.;

  init_fields = true;
  update_properties = true;
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
  //boundary(YList);
  //boundary(YLList);
  //boundary(YGList);
  //boundary(YGIntList);
  //boundary(YLIntList);
  //boundary(mEvapList);
  //boundary(slexpList);
  //boundary(slimpList);
  //boundary(sgexpList);
  //boundary(sgimpList);

#ifdef SOLVE_TEMPERATURE
  foreach() {
    TL[] = TL0*f[];
    TG[] = TG0*(1. - f[]);
    T[]  = TL[] + TG[];
  }
  //boundary({T,TL,TG});
#endif

#ifdef VARPROP
  if (is_constant (rho2v)) {
    scalar * l = list_copy (all);
    rho2v = new scalar;
    free (all);
    all = list_concat ({rho2v}, l);
    free (l);
  }
  if (is_constant (mu2v)) {
    scalar * l = list_copy (all);
    mu2v = new scalar;
    free (all);
    all = list_concat ({mu2v}, l);
    free (l);
  }
  if (is_constant (cp2v)) {
    scalar * l = list_copy (all);
    cp2v = new scalar;
    free (all);
    all = list_concat ({cp2v}, l);
    free (l);
  }
  if (is_constant (lambda2v)) {
    scalar * l = list_copy (all);
    lambda2v = new scalar;
    free (all);
    all = list_concat ({lambda2v}, l);
    free (l);
  }
  if (is_constant (rho1v)) {
    scalar * l = list_copy (all);
    rho1v = new scalar;
    free (all);
    all = list_concat ({rho1v}, l);
    free (l);
  }
  if (is_constant (mu1v)) {
    scalar * l = list_copy (all);
    mu1v = new scalar;
    free (all);
    all = list_concat ({mu1v}, l);
    free (l);
  }
  if (is_constant (cp1v)) {
    scalar * l = list_copy (all);
    cp1v = new scalar;
    free (all);
    all = list_concat ({cp1v}, l);
    free (l);
  }
  if (is_constant (lambda1v)) {
    scalar * l = list_copy (all);
    lambda1v = new scalar;
    free (all);
    all = list_concat ({lambda1v}, l);
    free (l);
  }
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
## Update Properties

Update non-constant properties. */

event properties (i++)
{
  if (update_properties) {

    if (!first_iter) {
      foreach() {
        frhocp1r0[] = f[]*rho1v[]*cp1v[];
        frhocp2r0[] = (1. - f[])*rho2v[]*cp2v[];
      }
    }

    /**
    Update non-constant properties. */

    double yliq[NLS], ygas[NGS];
    double xliq[NLS], xgas[NGS];
    double mwl[NLS];
    for (int jj=0; jj<NLS; jj++)
      mwl[jj] = inMW[LSI[jj]];

    scalar betaexp1[], betaexp2[];
    ts1.P = Pref;
    ts2.P = Pref;

    ts1.x = xliq;
    ts2.x = xgas;

    double T_PROP = 0.1;

    double mass_t = 0.;
    foreach(reduction(+:mass_t))
      mass_t += f[]*rho1v[]*dv();
#if AXI
    mass_t *= 2.*pi;
#endif

    foreach() {
      if (f[] > T_PROP) {

        double * T1h = &ts1.T;
        double * x1h = ts1.x;
        *T1h = TL[]/f[];

#ifdef FORCEINITPROP
        *T1h = TL0;
        for (int jj=0; jj<NLS; jj++)
          x1h[jj] = liq_start[jj];
#endif

        foreach_elem (YLList, jj) {
          scalar YL = YLList[jj];
          yliq[jj] = (NLS == 1) ? 1. : YL[]/f[];
          mwl[jj] = inMW[LSI[jj]];
        }
        mass2molefrac (x1h, yliq, mwl, NLS);

        //rho1v0[] = rho1v[];
        rho1v[] = rho1;
        mu1v[] = mu1;
        cp1v[] = cp1;
        lambda1v[] = lambda1;
        rho1v[] = tp1.rhov (&ts1);
        mu1v[] = tp1.muv (&ts1);
        cp1v[] = tp1.cpv (&ts1);
        lambda1v[] = tp1.lambdav (&ts1);
        betaexp1[] = liqprop_thermal_expansion (&tp1, &ts1);

        // We want the liquid phase diffusivity
        // to be weighted on the mass fractions
        double x1hbkp[NLS];
        foreach_elem (YLList, jj) {
          x1hbkp[jj] = x1h[jj];
          x1h[jj] = yliq[jj];
        }

        foreach_elem (YLList, jj) {
          // Enthalpy of evaporation
          scalar dhevjj = dhevList[jj];
          dhevjj[] = dhev;
          dhevjj[] = tp1.dhev (&ts1, jj);

          // Liquid phase diffusivity
          scalar Dmix1v = Dmix1List[jj];
          Dmix1v[] = inDmix1[jj];
          Dmix1v[] = tp1.diff (&ts1, jj);
        }

        // We recover the real mole fraction values
        foreach_elem (YLList, jj) {
          x1h[jj] = x1hbkp[jj];
        }

      }
      else {
        //rho1v0[] = 0.;
        rho1v[] = 0.;
        mu1v[] = 0.;
        cp1v[] = 0.;
        cp1v[] = 0.;
        lambda1v[] = 0.;
        betaexp1[] = 0.;

        foreach_elem (Dmix1List, jj) {
          scalar Dmix1v = Dmix1List[jj];
          scalar dhevjj = dhevList[jj];
          Dmix1v[] = 0.;
          dhevjj[] = 0.;
        }
      }
      if ((1. - f[]) > T_PROP) {

        double * T2h = &ts2.T;
        double * x2h = ts2.x;
        *T2h = TG[]/(1. - f[]);

#ifdef FORCEINITPROP
          *T2h = TG0;
          for (int jj=0; jj<NGS; jj++)
            x2h[jj] = gas_start[jj];
#endif

        foreach_elem (YGList, jj) {
          scalar YG = YGList[jj];
          ygas[jj] = YG[]/(1. - f[]);
        }
        mass2molefrac (x2h, ygas, inMW, NGS);

        rho2v0[] = rho2v[];
        rho2v[] = rho2;
        mu2v[] = mu2;
        cp2v[] = cp2;
        lambda2v[] = lambda2;
        rho2v[] = tp2.rhov (&ts2);
        mu2v[] = tp2.muv (&ts2);
        cp2v[] = tp2.cpv (&ts2);
        lambda2v[] = tp2.lambdav (&ts2);
        betaexp2[] = gasprop_thermal_expansion (&ts2);

        foreach_elem (Dmix2List, jj) {
          scalar Dmix2v = Dmix2List[jj];
          Dmix2v[] = inDmix2[jj];
          Dmix2v[] = tp2.diff (&ts2, jj);
        }
      }
      else {
        rho2v0[] = 0.;
        rho2v[] = 0.;
        mu2v[] = 0.;
        cp2v[] = 0.;
        betaexp2[] = 0.;

        foreach_elem (Dmix2List, jj) {
          scalar Dmix2v = Dmix2List[jj];
          Dmix2v[] = 0.;
        }
      }

      frho1r[] = f[]*rho1v[];
      frho2r[] = (1. - f[])*rho2v[];
      frhocp1r[] = f[]*rho1v[]*cp1v[];
      frhocp2r[] = (1. - f[])*rho2v[]*cp2v[];
      frho1[] = frho1r[];
      frho2[] = frho2r[];
      frhocp1[] = frhocp1r[];
      frhocp2[] = frhocp2r[];
    }

    foreach() {
      if (f[] <= T_PROP) {
        double rho1vgh = 0.;
        double mu1vgh = 0.;
        double cp1vgh = 0.;
        double lambda1vgh = 0.;
        double dhevgh[NLS];
        for (int jj=0; jj<NLS; jj++)
          dhevgh[jj] = 0.;

        int counter = 0;
        foreach_neighbor(1) {
          if (f[] > T_PROP) {
            counter++;
            rho1vgh += rho1v[];
            mu1vgh += mu1v[];
            cp1vgh += cp1v[];
            lambda1vgh += lambda1v[];

            for (int jj=0; jj<NLS; jj++) {
              scalar dhevjj = dhevList[jj];
              dhevgh[jj] += dhevjj[];
            }
          }
        }
        //rho1v0[] = rho1v[];
        rho1v[] = (counter != 0) ? rho1vgh/counter : 0.;
        mu1v[] = (counter != 0) ? mu1vgh/counter : 0.;
        cp1v[] = (counter != 0) ? cp1vgh/counter : 0.;
        lambda1v[] = (counter != 0) ? lambda1vgh/counter : 0.;

        for (int jj=0; jj<NLS; jj++) {
          scalar dhevjj = dhevList[jj];
          dhevjj[] = (counter != 0) ? dhevgh[jj]/counter : 0.;
        }
      }

      if ((1. - f[]) <= T_PROP) {
        double rho2vgh = 0.;
        double mu2vgh = 0.;
        double cp2vgh = 0.;
        double lambda2vgh = 0.;
        double Dmix2vgh[NGS];
        for (int jj=0; jj<NGS; jj++)
          Dmix2vgh[jj] = 0.;

        int counter = 0;
        foreach_neighbor(1) {
          if ((1. - f[]) > T_PROP) {
            counter++;
            rho2vgh += rho2v[];
            mu2vgh += mu2v[];
            cp2vgh += cp2v[];
            lambda2vgh += lambda2v[];

            for (int jj=0; jj<NGS; jj++) {
              scalar Dmix2jj = Dmix2List[jj];
              Dmix2vgh[jj] += Dmix2jj[];
            }
          }
        }
        rho2v0[] = rho2v[];
        rho2v[] = (counter != 0) ? rho2vgh/counter : 0.;
        mu2v[] = (counter != 0) ? mu2vgh/counter : 0.;
        cp2v[] = (counter != 0) ? cp2vgh/counter : 0.;
        lambda2v[] = (counter != 0) ? lambda2vgh/counter : 0.;

        for (int jj=0; jj<NGS; jj++) {
          scalar Dmix2jj = Dmix2List[jj];
          Dmix2jj[] = (counter != 0) ? Dmix2vgh[jj]/counter : 0.;
        }
      }

      frho1r[] = f[]*rho1v[];
      frhocp1r[] = f[]*rho1v[]*cp1v[];
      frho2r[] = (1. - f[])*rho2v[];
      frhocp2r[] = (1. - f[])*rho2v[]*cp2v[];
      frho1[] = frho1r[];
      frho2[] = frho2r[];
      frhocp1[] = frhocp1r[];
      frhocp2[] = frhocp2r[];
    }

  //  if (!first_iter) {
  //    double mass_tp1 = 0.;
  //    foreach(reduction(+:mass_tp1))
  //      mass_tp1 += f[]*rho1v[]*dv();
  //#if AXI
  //    mass_tp1 *= 2.*pi;
  //#endif
  //    double dmassdt = (mass_tp1 - mass_t)/dt;
  //
  //    double segment = 0.;
  //    foreach(reduction(+:segment)) {
  //      if (f[] > F_ERR && f[] < 1.-F_ERR) {
  //        coord n = mycs (point, f), p;
  //        double alpha = plane_alpha (f[], n);
  //        double area = plane_area_center (n, alpha, &p);
  //#if AXI
  //        segment += area*(y + p.y*Delta);
  //#else
  //        segment += area;
  //#endif
  //      }
  //    }
  //    foreach() {
  //      mDilation[] = 0.;
  //      if (f[] > F_ERR && f[] < 1.-F_ERR) {
  //        mDilation[] = dmassdt/segment;
  //      }
  //    }
  //  }


    // Compute lagrangian derivative of density
    scalar rhovt[], cpvt[], betavt[], lambdavt[];
    foreach() {
      rhovt[] = aavg (f[], rho1v[], rho2v[]);
      betavt[] = aavg (f[], betaexp1[], betaexp2[]);
      lambdavt[] = aavg (f[], lambda1v[], lambda2v[]);
      cpvt[] = aavg (f[], cp1v[], cp2v[]);

      TL[] = (f[] > F_ERR) ? TL[]/f[] : 0.;
      TG[] = (1. - f[] > F_ERR) ? TG[]/(1. - f[]) : 0.;
    }

    face_fraction (f, fsL);

    face vector lambdagT[];
    foreach_face() {
      double lambdavf = 0.5*(lambdavt[] + lambdavt[-1]);
      lambdagT.x[] = fm.x[]*lambdavf*face_gradient_x (T, 0); /// !<<
    }
    //face vector lambdagT[], lambdalT[];
    //foreach_face() {
    //  double lambda1vf = 0.5*(lambda1v[] + lambda1v[-1]);
    //  lambdalT.x[] = fm.x[]*fsL.x[]*lambda1vf*face_gradient_x (TL, 0);
    //  double lambda2vf = 0.5*(lambda2v[] + lambda2v[-1]);
    //  lambdagT.x[] = fm.x[]*(1. - fsL.x[])*lambda2vf*face_gradient_x (TG, 0);
    //}

    //face vector rhovflux[];
    //tracer_fluxes (frho1r, uf, rhovflux, dt, zeroc);

    foreach() {
      double laplT = 0.;
      foreach_dimension()
        laplT += (lambdagT.x[1] - lambdagT.x[]);
      laplT /= Delta;
      //double laplT1 = 0., laplT2 = 0.;
      //foreach_dimension() {
      //  laplT1 += (lambdalT.x[1] - lambdalT.x[]);
      //  laplT2 += (lambdagT.x[1] - lambdagT.x[]);
      //}
      //laplT1 /= Delta;
      //laplT2 /= Delta;

      double drho1dt = (f[] > F_ERR) ?
        -betaexp1[]/(rho1v[]*cp1v[])*laplT : 0.;

      //double drho2dt = ((1. - f[]) > F_ERR) ?
      double drho2dt = ((1. - f[]) > 1.-F_ERR) ?
        -1./(rho2v[]*cp2v[]*T[])*laplT : 0.;

      //double drho1dt = (f[] > F_ERR) ?
      //  betaexp1[]/(rho1v[]*cp1v[])*laplT1 : 0.;

      //double drho2dt = ((1. - f[]) > F_ERR) ?
      //  -betaexp2[]/(rho2v[]*cp2v[])*laplT2 : 0.;

      drhodt[] = aavg (f[], drho1dt, 0.);

      TL[] *= f[];
      TG[] *= (1. - f[]);

      //// New calculation using explicit density
      //double temporal1 = (rho1v[] - rho1v0[])/dt;
      //double temporal2 = (rho2v[] - rho2v0[])/dt;
      //double temporal = aavg (f[], temporal1, temporal2);

      //double div_rhov = 0.;
      //foreach_dimension()
      //  div_rhov += (rhovflux.x[1] - rhovflux.x[]);
      //div_rhov /= (Delta*cm[]);

      //double div = 0.;
      //foreach_dimension()
      //  div += (uf.x[1] - uf.x[]);
      //div /= Delta;

      ////drhodt[] = 1./rhovt[]*(temporal + div_rhov - rhovt[]*div);
      ////dummy[] = (f[] > 1.e-3) ? 1./rho1v[]*(temporal + div_rhov - rho1v[]*div) : 0.;
      //dummy[] = (f[] > F_ERR) ? f[]*(1./rho1v[]*((rho1v[] - rho1v0[])/dt*cm[] + div_rhov - rho1v[]*div)) : 0.;
      //drhodt[] = -dummy[];
    }
  }
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
        double dhevvh = dhev;
# ifdef VARPROP
        scalar dhevjj = dhevList[jj];
        dhevvh = dhevjj[];
# endif
        XGIntConv[jj] = clapeyron (min (TInt[], Tboil[jj]-1.), Tboil[jj], dhevvh, inMW[jj])*XLIntConv[jj];
#endif
#ifdef USE_ANTOINE
        scalar YL = YLList[jj];
        XGIntConv[jj] = min (YL.antoine (TInt[], Pref), 0.98)*XLIntConv[jj];
        //XGIntConv[jj] = min (opensmoke_antoine (TInt[], Pref, jj), 0.98)*XLIntConv[jj];
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

        double rho2vh = rho2;
        double inDmix2vh = inDmix2[LSI[jj]];
#ifdef VARPROP
        rho2vh = rho2v[];
        //scalar Dmix2v = Dmix2List[jj];
        scalar Dmix2v = Dmix2List[LSI[jj]];
        inDmix2vh = Dmix2v[];
#endif
        double gtrgrad = ebmgrad (point, YG, fL, fG, fsL, fsG, true, YGInt[], &success);
        sum_jG += -rho2vh*inDmix2vh*gtrgrad;
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

        double rho2vh = rho2;
        double inDmix2vh = inDmix2[LSI[jj]];
#ifdef VARPROP
        rho2vh = rho2v[];
        //scalar Dmix2v = Dmix2List[jj];
        scalar Dmix2v = Dmix2List[LSI[jj]];
        inDmix2vh = Dmix2v[];
#endif
        double gtrgrad = ebmgrad (point, YG, fL, fG, fsL, fsG, true, YGInt[], &success);
#ifdef DIFFUSIVE
        mEvap[] = - rho2vh*inDmix2vh*gtrgrad;
#else
        mEvap[] = mEvapTot[]*YGInt[] - rho2vh*inDmix2vh*gtrgrad;
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

#ifdef USE_GSL

# ifdef SOLVE_TEMPERATURE
  ijc_CoupledTemperature();
# endif

  /**
  Finally, we can refine the first guess values from the decoupled
  solution from the fully-coupled interface jump condition. Solved
  as a non-linear system of equations. */

# if USE_GSL > 0
  ijc_CoupledNls();
# endif

#endif

  /**
  Check if boiling conditions occur and correct the
  evaporation rate switching to the boiling regime. */

#ifdef POSSIBLE_BOILING
  if (NLS == 1) {
    scalar mBoil[];
    foreach() {
      scalar dhev0 = dhevList[0];
      scalar YL    = YLList[0];
      scalar YLInt = YLIntList[0];
      scalar YGInt = YGIntList[LSI[0]];
      scalar mEvap = mEvapList[LSI[0]];

      if (f[] > F_ERR && f[] < 1.-F_ERR) {
        TInt[] = 373.;
        double gtrgrad = ebmgrad (point, TG, fL, fG, fsL, fsG, true, TInt[], false);
        mBoil[] = lambda2v[]*gtrgrad/dhev0[];
        double Keq = YL.antoine (TInt[], Pref);
        double sigmoid = max ( 2.*(1./(1. + exp (-20.*(Keq-0.9)))-0.5), 0.);
        mEvap[] = mEvap[]*(1.-sigmoid) + mBoil[]*sigmoid;
        YLInt[] = 1.;
        YGInt[] = YGInt[]*(1.-sigmoid) + YLInt[]*sigmoid;
      }
    }
  }
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

        double rho2vh = rho2;
#ifdef VARPROP
        rho2vh = rho2v[];
#endif

#ifdef AXI
        //sgexp[] = -mEvap[]/rho2vh*area*(y + p.y*Delta)/(Delta*y)*cm[];
        //sgimp[] = +mEvapTot[]/rho2vh*area*(y + p.y*Delta)/(Delta*y)*cm[];
        sgexp[] = -mEvap[]*area*(y + p.y*Delta)/(Delta*y)*cm[];
        sgimp[] = +mEvapTot[]*area*(y + p.y*Delta)/(Delta*y)*cm[];
#else
        //sgexp[] = -mEvap[]/rho2vh*area/Delta*cm[];
        //sgimp[] = +mEvapTot[]/rho2vh*area/Delta*cm[];
        sgexp[] = -mEvap[]*area/Delta*cm[];
        sgimp[] = +mEvapTot[]*area/Delta*cm[];
#endif
      }

      for (int jj=0; jj<NLS; jj++) {
        scalar slexp = slexpList[jj];
        scalar slimp = slimpList[jj];
        scalar mEvap = mEvapList[LSI[jj]];

        double rho1vh = rho1;
#ifdef VARPROP
        rho1vh = rho1v[];
#endif

#ifdef AXI
        //slexp[] = +mEvap[]/rho1vh*area*(y + p.y*Delta)/(Delta*y)*cm[];
        //slimp[] = -mEvapTot[]/rho1vh*area*(y + p.y*Delta)/(Delta*y)*cm[];
        slexp[] = +mEvap[]*area*(y + p.y*Delta)/(Delta*y)*cm[];
        slimp[] = -mEvapTot[]*area*(y + p.y*Delta)/(Delta*y)*cm[];
#else
        //slexp[] = +mEvap[]/rho1vh*area/Delta*cm[];
        //slimp[] = -mEvapTot[]/rho1vh*area/Delta*cm[];
        slexp[] = +mEvap[]*area/Delta*cm[];
        slimp[] = -mEvapTot[]*area/Delta*cm[];
#endif
      }
    }
  }

#ifdef SOLVE_TEMPERATURE
  /**
  Compute source terms for temperature equations. */

  foreach() {
    sgT[] = 0.; slT[] = 0.;
    sgTimp[] = 0.; slTimp[] = 0.;
    if (f[] > F_ERR && f[] < 1.-F_ERR) { // <<!!
      coord n = facet_normal (point, fL, fsL), p;
      double alpha = plane_alpha (fL[], n);
      double area = plane_area_center (n, alpha, &p);
      normalize (&n);

      double bc = TInt[];
      double gtrgrad = ebmgrad (point, TG, fL, fG, fsL, fsG, true, bc, &success);
      double ltrgrad = ebmgrad (point, TL, fL, fG, fsL, fsG, false, bc, &success);

      double rho1vh = rho1;
      double rho2vh = rho2;
      double cp1vh = cp1;
      double cp2vh = cp2;
      double lambda1vh = lambda1;
      double lambda2vh = lambda2;
#ifdef VARPROP
      rho1vh = rho1v[];
      rho2vh = rho2v[];
      cp1vh = cp1v[];
      cp2vh = cp2v[];
      lambda1vh = lambda1v[];
      lambda2vh = lambda2v[];
#endif

      double lheatflux = lambda1vh*ltrgrad;
      double gheatflux = lambda2vh*gtrgrad;

#ifdef AXI
      slT[] = lheatflux/rho1vh/cp1vh*area*(y + p.y*Delta)/(Delta*y)*cm[];
      sgT[] = gheatflux/rho2vh/cp2vh*area*(y + p.y*Delta)/(Delta*y)*cm[];
      slTimp[] = -mEvapTot[]/rho1vh/cp1vh*area*(y + p.y*Delta)/(Delta*y)*cm[];
      sgTimp[] = mEvapTot[]/rho2vh/cp2vh*area*(y + p.y*Delta)/(Delta*y)*cm[];
      //slT[] = lheatflux*area*(y + p.y*Delta)/(Delta*y)*cm[];
      //sgT[] = gheatflux*area*(y + p.y*Delta)/(Delta*y)*cm[];
      //slTimp[] = -mEvapTot[]*area*(y + p.y*Delta)/(Delta*y)*cm[];
      //sgTimp[] = mEvapTot[]*area*(y + p.y*Delta)/(Delta*y)*cm[];
#else
      slT[] = lheatflux/rho1vh/cp1vh*area/Delta*cm[];
      sgT[] = gheatflux/rho2vh/cp2vh*area/Delta*cm[];
      slTimp[] = -mEvapTot[]/rho1vh/cp1vh*area/Delta*cm[];
      sgTimp[] = mEvapTot[]/rho2vh/cp2vh*area/Delta*cm[];
      //slT[] = lheatflux*area/Delta*cm[];
      //sgT[] = gheatflux*area/Delta*cm[];
      //slTimp[] = -mEvapTot[]*area/Delta*cm[];
      //sgTimp[] = mEvapTot[]*area/Delta*cm[];
#endif
    }
  }
#endif

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
    T[] = TL[] + TG[];
#endif
  }
}

/**
## Tracer Advection

We let the volume fractions *fu* and *fuext* to
advect the fields YL and YG, as implemented in
the tracer_advection event of [evaporation.h](evaporation.h)
*/

event tracer_advection (i++)
{
  foreach() {
    for (scalar YL in YLList) {
      YL[] = (f[] > F_ERR) ? YL[]/f[] : 0.;
      YL[] = f[]*rho1v[]*YL[];
    }
    for (scalar YG in YGList) {
      YG[] = ((1. - f[]) > F_ERR) ? YG[]/(1. - f[]) : 0.;
      YG[] = (1. - f[])*rho2v[]*YG[];
    }
    //TL[] = (f[] > F_ERR) ? TL[]/f[] : 0.;
    //TL[] = frhocp1[]*TL[];
    //TG[] = ((1. - f[]) > F_ERR) ? TG[]/(1. - f[]) : 0.;
    //TG[] = frhocp2[]*TG[];
  }
}

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
      //YL[] = (fuext[] > F_ERR) ? YL[]/fuext[] : 0.;
      YL[] = (fuext[] > F_ERR) ? YL[]/frho1[] : 0.;
    for (scalar YG in YGList)
      //YG[] = ((1. - fu[]) > F_ERR) ? YG[]/(1. - fu[]) : 0.;
      YG[] = ((1. - fu[]) > F_ERR) ? YG[]/frho2[] : 0.;

    //fL[] = f[]; fG[] = 1. - f[];

#ifdef VARPROP
    rho1v[] = (f[] > F_ERR) ? frho1r[]/f[] : 0.;
    rho2v[] = ((1. - f[]) > F_ERR) ? frho2r[]/(1. - f[]) : 0.;
    cp1v[] = (f[] > F_ERR) ? frhocp1r[]/(f[]*rho1v[]) : 0.;
    cp2v[] = ((1. - f[]) > F_ERR) ? frhocp2r[]/((1. - f[])*rho2v[]) : 0.;
#endif

#ifdef SOLVE_TEMPERATURE
    TL[] = (fuext[] > F_ERR) ? TL[]/fuext[] : 0.;
    TG[] = ((1. - fu[]) > F_ERR) ? TG[]/(1. - fu[]) : 0.;
    //TL[] = (fuext[] > F_ERR) ? TL[]/frhocp1[] : 0.;
    //TG[] = ((1. - fu[]) > F_ERR) ? TG[]/frhocp2[] : 0.;
#endif

    fL[] = f[]; fG[] = 1. - f[];
  }

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
  theta1.dirty = true;
  theta2.dirty = true;
#endif

  for (int jj=0; jj<NLS; jj++) {
    face vector Dmix1f[];
    foreach_face() {
      double Dmix1vh = inDmix1[jj];
      double rho1vh = rho1;
#ifdef VARPROP
      scalar Dmix1jj = Dmix1List[jj];
      Dmix1vh = 0.5*(Dmix1jj[] + Dmix1jj[-1]);
      rho1vh = 0.5*(rho1v[] + rho1v[-1]);
#endif
      Dmix1f.x[] = rho1vh*Dmix1vh*fsL.x[]*fm.x[];
    }

    foreach()
      //theta1[] = cm[]*max(fL[]*rho1v[], 1.e-3);
      theta1[] = cm[]*max(frho1r[], F_ERR);
      //theta1[] = cm[]*max(fL[], T_ERR);

    scalar YL = YLList[jj];
    scalar slexp = slexpList[jj];
    scalar slimp = slimpList[jj];

    diffusion (YL, dt, D=Dmix1f, r=slexp, beta=slimp, theta=theta1);
  }

  for (int jj=0; jj<NGS; jj++) {

    face vector Dmix2f[];
    foreach_face() {
      double Dmix2vh = inDmix2[jj];
      double rho2vh = rho2;
#ifdef VARPROP
      scalar Dmix2v = Dmix2List[jj];
      Dmix2vh = 0.5*(Dmix2v[] + Dmix2v[-1]);
      rho2vh = 0.5*(rho2v[] + rho2v[-1]);
#endif
      Dmix2f.x[] = rho2vh*Dmix2vh*fsG.x[]*fm.x[];
    }

    foreach()
      //theta2[] = cm[]*max(fG[]*rho2v[], F_ERR);
      theta2[] = cm[]*max(frho2r[], F_ERR);

    scalar YG = YGList[jj];
    scalar sgexp = sgexpList[jj];
    scalar sgimp = sgimpList[jj];

    diffusion (YG, dt, D=Dmix2f, r=sgexp, beta=sgimp, theta=theta2);
  }

#ifdef SOLVE_TEMPERATURE

  foreach_face() {
#ifdef VARPROP
    double alpha1l = (rho1v[] != 0.) ? lambda1v[]/rho1v[]/cp1v[] : 0.;
    double alpha2l = (rho2v[] != 0.) ? lambda2v[]/rho2v[]/cp2v[] : 0.;
    double alpha1r = (rho1v[-1] != 0.) ? lambda1v[-1]/rho1v[-1]/cp1v[-1] : 0.;
    double alpha2r = (rho2v[-1] != 0.) ? lambda2v[-1]/rho2v[-1]/cp2v[-1] : 0.;
    lambda1f.x[] = 0.5*(alpha1r + alpha1l)*fsL.x[]*fm.x[];
    lambda2f.x[] = 0.5*(alpha2r + alpha2l)*fsG.x[]*fm.x[];
    //lambda1f.x[] = 0.5*(lambda1v[] + lambda1v[-1])*fsL.x[]*fm.x[];
    //lambda2f.x[] = 0.5*(lambda2v[] + lambda1v[-1])*fsG.x[]*fm.x[];
#else
    lambda1f.x[] = lambda1/rho1/cp1*fsL.x[]*fm.x[];
    lambda2f.x[] = lambda2/rho2/cp2*fsG.x[]*fm.x[];
#endif
  }

  foreach() {
    theta1[] = cm[]*max(fL[], F_ERR);
    theta2[] = cm[]*max(fG[], F_ERR);
    //theta1[] = cm[]*max(frhocp1[], F_ERR);
    //theta2[] = cm[]*max(frhocp2[], F_ERR);
  }

  /**
  Solve diffusion equations for temperature. */

  diffusion (TL, dt, D=lambda1f, r=slT, theta=theta1);
  diffusion (TG, dt, D=lambda2f, r=sgT, theta=theta2);

#endif

  foreach() {
    /**
    We reconstruct the mass fractions from the total mass
    in every cell. */

    double rho1vh = rho1;
    double rho2vh = rho2;
#ifdef VARPROP
    rho1vh = rho1v[];
    rho2vh = rho2v[];
#endif

    double totmassliq = 0.;
    for (scalar YL in YLList)
      totmassliq += YL[];
    for (scalar YL in YLList)
      YL[] = (totmassliq > 0.) ? YL[] / totmassliq : 0.;

    double totmassgas = 0.;
    for (scalar YG in YGList)
      totmassgas += YG[];
    for (scalar YG in YGList)
      YG[] = (totmassgas > 0.) ? YG[] / totmassgas : 0.;

    for (scalar YL in YLList)
      YL[] *= f[];
    for (scalar YG in YGList)
      YG[] *= (1. - f[]);
#ifdef SOLVE_TEMPERATURE
    TL[] *= f[];
    TG[] *= (1. - f[]);
#endif
  }

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

#ifdef SOLVE_TEMPERATURE
  foreach()
    T[] = TL[] + TG[];
#endif
  update_properties = true;
}

event end_timestep (i++) {
  update_properties = false;
  first_iter = false;
}

