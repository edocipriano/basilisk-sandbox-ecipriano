/**
# Update Properties

Update the thermodynamic and transport properties for the
multicomponent phase change model, and compute the density
material derivative which is used as a source term for the
velocity divergence in order to describe low-Mach compressibility
effects. */

#ifndef T_PROP
# define T_PROP 0.1
#endif

#ifdef VARPROP

scalar DTDt1[], DTDt2[];
scalar * DYDt1 = NULL;    // [NLS]
scalar * DYDt2 = NULL;    // [NGS]

scalar Hcheck[];
scalar betaexp1[], betaexp2[];
scalar rho1vInt[], rho2vInt[];

/**
## *update_properties_constant()*: update variable property fields setting constant properties (for debug purposes)
*/

void update_properties_constant (void) {
  foreach() {
    rho1v[] = rho1;
    rho1v0[] = rho1;
    rho1vInt[] = rho1;
    mu1v[] = mu1;
    cp1v[] = cp1;
    lambda1v[] = lambda1;

    for (int jj=0; jj<NLS; jj++) {
      scalar dhevjj = dhevList[jj];
      dhevjj[] = dhev;

      scalar Dmix1v = Dmix1List[jj];
      Dmix1v[] = inDmix1[jj];

      scalar Cp1v = Cp1List[jj];
      Cp1v[] = cp1;
    }

    rho2v[] = rho2;
    rho2v0[] = rho2;
    mu2v[] = mu2;
    cp2v[] = cp2;
    lambda2v[] = lambda2;

    for (int jj=0; jj<NGS; jj++) {
      scalar Dmix2v = Dmix2List[jj];
      Dmix2v[] = inDmix2[jj];

      scalar Cp2v = Cp2List[jj];
      Cp2v[] = cp2;
    }
  }
}

/**
## *update_properties_initial()*: update variable properties with initial conditions
*/

void update_properties_initial (void) {

  foreach() {
    ThermoState ts1h, ts2h;
    ts1h.T = TL0;
    ts2h.T = TG0;
    ts1h.P = Pref;
    ts2h.P = Pref;
    ts1h.x = liq_start;
    ts2h.x = gas_start;

    rho1v[] = tp1.rhov (&ts1h);
    rho1vInt[] = rho1v[];
    rho1v0[] = rho1v[];
    mu1v[] = tp1.muv (&ts1h);
    cp1v[] = tp1.cpv (&ts1h);
    lambda1v[] = tp1.lambdav (&ts1h);

    for (int jj=0; jj<NLS; jj++) {
      scalar dhevjj = dhevList[jj];
      dhevjj[] = tp1.dhev (&ts1h, jj);

      scalar Dmix1v = Dmix1List[jj];
      Dmix1v[] = tp1.diff (&ts1h, jj);

      scalar Cp1v = Cp1List[jj];
      Cp1v[] = tp1.cps (&ts1h, jj);
    }

    rho2v[] = tp2.rhov (&ts2h);
    rho2v0[] = rho2v[];
    mu2v[] = tp2.muv (&ts2h);
    cp2v[] = tp2.cpv (&ts2h);
    lambda2v[] = tp2.lambdav (&ts2h);

    for (int jj=0; jj<NGS; jj++) {
      scalar Dmix2v = Dmix2List[jj];
      Dmix2v[] = tp2.diff (&ts2h, jj);

      scalar Cp2v = Cp2List[jj];
      Cp2v[] = tp2.cps (&ts2h, jj);
    }
  }
}


event defaults (i = 0)
{
  for (int jj=0; jj<NLS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "DYDt1_%s", liq_species[jj]);
    a.name = strdup (name);
    a.nodump = true;
    DYDt1 = list_append (DYDt1, a);
  }

  for (int jj=0; jj<NGS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "DYDt2_%s", gas_species[jj]);
    a.name = strdup (name);
    a.nodump = true;
    DYDt2 = list_append (DYDt2, a);
  }
}

double mLiq0 = 0.;

/**
## Initialization

We set update the variable property fields as a function
of the initial conditions. */

event init (i = 0)
{
  update_properties_initial();

  mLiq0 = 0.;
  foreach(reduction(+:mLiq0))
    mLiq0 += rho1v[]*f[]*dv();

  MW1mix.dirty = true;
  MW2mix.dirty = true;

#if TREE
  for (scalar s in {drhodt, drhodtext}) {
#if EMBED
    s.refine = s.prolongation = refine_embed_linear;
#else
    s.refine  = refine_linear;
#endif
    s.restriction = restriction_volume_average;
    s.dirty = true; // boundary conditions need to be updated
  }
#endif
}

event cleanup (t = end)
{
  delete (DYDt1), free (DYDt1), DYDt1 = NULL;
  delete (DYDt2), free (DYDt2), DYDt2 = NULL;
}

/**
The following event resets the material derivatives, so
that they can be modified by external modules by adding
source terms which contribute to the density variations.
*/

event reset_sources (i++)
{
  foreach() {
    DTDt1[] = 0.;
    DTDt2[] = 0.;

    for (int jj=0; jj<NLS; jj++) {
      scalar DYDt1jj = DYDt1[jj];
      DYDt1jj[] = 0.;
    }
    for (int jj=0; jj<NGS; jj++) {
      scalar DYDt2jj = DYDt2[jj];
      DYDt2jj[] = 0.;
    }
  }
}

/**
## *update_properties()*: update variable properties as a function of the thermodynamic state

Every phase property is calculated from expressions and correlations as
a function of the thermodynamic state of the mixture (in this case temperature,
thermodynamic pressure, and mole fractions):

$$
\phi_k = f(T_k, P , x_{i,k})
$$

These functions can be easily overwritten by assigning the `ThermoProp` 
functions in [variable-properties.h](/sandbox/ecipriano/src/variable-properties.h#thermodynamic-properties).
*/

void update_properties (void)
{
  foreach() {
    rho1v0[] = rho1v[];
    rho2v0[] = rho2v[];
  }

  double MW1[NLS], MW2[NGS];
  for (int jj=0; jj<NLS; jj++)
    MW1[jj] = inMW[LSI[jj]];
  for (int jj=0; jj<NGS; jj++)
    MW2[jj] = inMW[jj];

  foreach() {
    //MW1mix[] = 0.;
    //MW2mix[] = 0.;

    if (f[] > T_PROP) {
      double x1[NLS], y1[NLS];
      for (int jj=0; jj<NLS; jj++) {
        scalar YL = YLList[jj];
        y1[jj] = (NLS == 1) ? 1. : YL[];
      }
      correctfrac (y1, NLS);
      mass2molefrac (x1, y1, MW1, NLS);
      //MW1mix[] = mass2mw (y1, MW1, NLS);

      ThermoState ts1h;
      ts1h.T = TL[];
      ts1h.P = Pref;
      ts1h.x = x1;

      rho1v[] = tp1.rhov (&ts1h);
      mu1v[] = tp1.muv (&ts1h);
      cp1v[] = tp1.cpv (&ts1h);
      lambda1v[] = tp1.lambdav (&ts1h);
      betaexp1[] = liqprop_thermal_expansion (&tp1, &ts1h);

      for (int jj=0; jj<NLS; jj++) {
        // Enthalpy of evaporation
        scalar dhevjj = dhevList[jj];
        dhevjj[] = tp1.dhev (&ts1h, jj);

        // Liquid phase diffusivity
        scalar Dmix1v = Dmix1List[jj];
        Dmix1v[] = tp1.diff (&ts1h, jj);

        // Liquid phase species heat capacity
        scalar Cp1v = Cp1List[jj];
        Cp1v[] = tp1.cps (&ts1h, jj);
      }
    }

    if ((1. - f[]) > T_PROP) {
      double x2[NGS], y2[NGS];
      for (int jj=0; jj<NGS; jj++) {
        scalar YG = YGList[jj];
        y2[jj] = YG[];
      }
      correctfrac (y2, NGS);
      mass2molefrac (x2, y2, MW2, NGS);
      //MW2mix[] = mass2mw (y2, MW2, NGS);

      ThermoState ts2h;
      ts2h.T = TG[];
      ts2h.P = Pref;
      ts2h.x = x2;

      rho2v[] = tp2.rhov (&ts2h);
      mu2v[] = tp2.muv (&ts2h);
      cp2v[] = tp2.cpv (&ts2h);
      lambda2v[] = tp2.lambdav (&ts2h);
      betaexp2[] = gasprop_thermal_expansion (&ts2h);

      for (int jj=0; jj<NGS; jj++) {
        scalar Dmix2v = Dmix2List[jj];
        Dmix2v[] = tp2.diff (&ts2h, jj);

        scalar Cp2v = Cp2List[jj];
        Cp2v[] = tp2.cps (&ts2h, jj);
      }
    }
  }

  foreach() {
    if (f[] <= T_PROP) {
      double rho1vgh = 0.;
      double mu1vgh = 0.;
      double cp1vgh = 0.;
      double lambda1vgh = 0.;
      double dhevgh[NLS];
      double beta1vgh = 0.;
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
          beta1vgh += betaexp1[];

          for (int jj=0; jj<NLS; jj++) {
            scalar dhevjj = dhevList[jj];
            dhevgh[jj] += dhevjj[];
          }
        }
      }
      rho1v[] = (counter != 0) ? rho1vgh/counter : 0.;
      mu1v[] = (counter != 0) ? mu1vgh/counter : 0.;
      cp1v[] = (counter != 0) ? cp1vgh/counter : 0.;
      lambda1v[] = (counter != 0) ? lambda1vgh/counter : 0.;
      betaexp1[] = (counter != 0) ? beta1vgh/counter : 0.;

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
      rho2v[] = (counter != 0) ? rho2vgh/counter : 0.;
      mu2v[] = (counter != 0) ? mu2vgh/counter : 0.;
      cp2v[] = (counter != 0) ? cp2vgh/counter : 0.;
      lambda2v[] = (counter != 0) ? lambda2vgh/counter : 0.;

      for (int jj=0; jj<NGS; jj++) {
        scalar Dmix2jj = Dmix2List[jj];
        Dmix2jj[] = (counter != 0) ? Dmix2vgh[jj]/counter : 0.;
      }
    }
  }
}

/**
## *update_divergence()*: update density material derivative

According to the low-Mach approximation, the density material derivative
is computed by considering that temperature and composition gradients
dominate over thermodynamic pressure gradients. The final expression of
the velocity divergence reads:

$$
\nabla\cdot\mathbf{u}_k =
\left. \dfrac{1}{\rho}\dfrac{D\rho}{D t} \right\vert_k = 
\beta_k\dfrac{D T_k}{D t}
+ M_{w,k} \sum_{i=1}^{NS} \dfrac{1}{M_{w,i}}\dfrac{D\omega_{i,k}}{D t}
$$

which is used both for the liquid and for the gas phase. However, in
liquid phase we assume that the density changes due to the temperature
gradients dominate over composition effects. This is always true for
pure liquid droplets, but it is generally valid also for multicomponent
droplets. Relaxing this hypotesis is not complex, it requires the last
term on the RHS to be computed without the ideal gas approximation.
*/

void update_divergence (void) {

  /**
  We define the variables used to compute the lagrangian derivative
  on each level. */

  restriction ({T,TL,TG});
  restriction (YLList);
  restriction (YGList);
#ifdef MOLAR_DIFFUSION
  restriction (XLList);
  restriction (XGList);
#endif

  /**
  We calculate the Lagrangian derivative of the temperature fields. */

  face vector lambdagradT1[], lambdagradT2[];
  foreach_face() {
    lambdagradT1.x[] = face_value (lambda1v, 0)*face_gradient_x (TL, 0)*fm.x[]*fsL.x[];
    lambdagradT2.x[] = face_value (lambda2v, 0)*face_gradient_x (TG, 0)*fm.x[]*fsG.x[];
  }

  foreach() {
    foreach_dimension()
      DTDt1[] += (lambdagradT1.x[1] - lambdagradT1.x[])/Delta;
    DTDt1[] += slT[];

    foreach_dimension()
      DTDt2[] += (lambdagradT2.x[1] - lambdagradT2.x[])/Delta;
    DTDt2[] += sgT[];
  }

  /**
  We calculate the Lagrangian derivative for the chemical species mass
  fractions. */

  for (int jj=0; jj<NGS; jj++) {
    scalar YG = YGList[jj];
    scalar Dmix2v = Dmix2List[jj];
    scalar DYDt2j = DYDt2[jj];

    face vector rhoDmixYGjj[];
    foreach_face() {
      double rho2f = face_value (rho2v, 0);
      double Dmix2f = face_value (Dmix2v, 0);
      rhoDmixYGjj.x[] = rho2f*Dmix2f*face_gradient_x (YG, 0)*fm.x[]*fsG.x[];
    }

    scalar sgexp = sgexpList[jj];
    scalar sgimp = sgimpList[jj];
    scalar YGInt = YGIntList[jj];

    foreach() {
      foreach_dimension()
        DYDt2j[] += (rhoDmixYGjj.x[1] - rhoDmixYGjj.x[])/Delta;
      DYDt2j[] += (sgexp[] + sgimp[]*YGInt[]);
    }
  }

  /**
  We add diffusion correction contributions to the chemical species
  mass fraction derivatives. */

  face vector phicGtot[];
  foreach_face() {
    phicGtot.x[] = 0.;
#ifdef FICK_CORRECTED
    for (int jj=0; jj<NGS; jj++) {
      scalar Dmix2v = Dmix2List[jj];

      double rho2f = face_value (rho2v, 0);
      double Dmix2f = face_value (Dmix2v, 0);
# ifdef MOLAR_DIFFUSION
      double MW2mixf = face_value (MW2mix, 0);

      scalar XG = XGList[jj];
      phicGtot.x[] += (MW2mixf > 0.) ?
        rho2f*Dmix2f*inMW[jj]/MW2mixf*face_gradient_x (XG, 0)*fm.x[]*fsG.x[] : 0.;
# else
      scalar YG = YGList[jj];
      phicGtot.x[] += rho2f*Dmix2f*face_gradient_x (YG, 0)*fm.x[]*fsG.x[];
# endif // MOLAR_DIFFUSION
    }
#endif  // FICK_CORRECTED
  }

  for (int jj=0; jj<NGS; jj++) {
    face vector phicGjj[];
    foreach_face() {
      phicGjj.x[] = phicGtot.x[];
#ifdef MOLAR_DIFFUSION
      scalar Dmix2v = Dmix2List[jj];

      double rho2f = face_value (rho2v, 0);
      double Dmix2f = face_value (Dmix2v, 0);
      double MW2mixf = face_value (MW2mix, 0);

      phicGjj.x[] -= (MW2mixf > 0.) ?
        rho2f*Dmix2f/MW2mixf*face_gradient_x (MW2mix, 0)*fm.x[]*fsG.x[] : 0.;
#endif

      scalar YG = YGList[jj];
      phicGjj.x[] *= face_value (YG, 0);
    }

    scalar DYDt2j = DYDt2[jj];

    foreach()
      foreach_dimension()
        DYDt2j[] -= (phicGjj.x[1] - phicGjj.x[])/Delta;
  }

  /**
  We calculate the one-field divergence by volume-averaging the liquid and the
  gas-phase contributions. */

  foreach() {
    double divu1 = 0., divu2 = 0.;

    // Add liquid temperature contribution
    divu1 += (rho1v[]*cp1v[] > 0.) ?
      betaexp1[]/(rho1v[]*cp1v[])*DTDt1[] : 0.;

    // Add gas temperature contribution
    divu2 += (TG[]*rho2v[]*cp2v[] > 0.) ?
      1./(TG[]*rho2v[]*cp2v[])*DTDt2[] : 0.;

    // Add gas chemical species contribution
    double divu2species = 0.;
    for (int jj=0; jj<NGS; jj++) {
      scalar DYDt2j = DYDt2[jj];
      divu2species += 1./inMW[jj]*DYDt2j[];
    }
    divu2 += (rho2v[] > 0.) ? MW2mix[]/rho2v[]*divu2species : 0.;

    // Volume averaged contributions
    drhodt[] = divu1*f[] + divu2*(1. - f[]);
    drhodtext[] = divu1*f[];

    // Adjust sign for internal convention
    drhodt[] *= -1.;
    drhodtext[] *= -1.;
  }
}

void update_divergence_density (void) {
  vector grho1[], grho2[];
  gradients ({rho1v, rho2v}, {grho1, grho2});

  scalar DrhoDt1[], DrhoDt2[];
  foreach() {
    DrhoDt1[] = (rho1v[] - rho1v0[])/dt;
    DrhoDt2[] = (rho2v[] - rho2v0[])/dt;

    foreach_dimension() {
#ifdef VELOCITY_JUMP
      DrhoDt1[] += u1.x[]*grho1.x[];
      DrhoDt2[] += u2.x[]*grho2.x[];
#else
      DrhoDt1[] += uext.x[]*grho1.x[];
      DrhoDt2[] += u.x[]*grho2.x[];
#endif
    }

    DrhoDt1[] = DrhoDt1[]*cm[];
    DrhoDt2[] = DrhoDt2[]*cm[];

    double one_over_rho1 = (rho1v[] > 0.) ? 1./rho1v[] : 0.;
    double one_over_rho2 = (rho2v[] > 0.) ? 1./rho2v[] : 0.;

    if (iter > 1) {
      drhodt[] = (one_over_rho1*DrhoDt1[]*f[] + one_over_rho2*DrhoDt2[]*(1. - f[]));
      drhodtext[] = one_over_rho1*DrhoDt1[]*f[];
    }
  }
}

#endif
