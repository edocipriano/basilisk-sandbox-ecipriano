/**
## Update Properties

Update the thermodynamic properties for the multicomponent phase
change model, and compute the lagrangian derivative of the density,
which is used as a sorce term for the velocity divergence, to
describe low Mach compressibility effects. */

#ifndef T_PROP
# define T_PROP 0.1
#endif

scalar Hcheck[];
scalar betaexp1[], betaexp2[];
scalar rho1vInt[], rho2vInt[];

void update_properties_constant (void) {
  foreach() {
    ThermoState ts1h, ts2h;
    ts1h.T = TL0;
    ts2h.T = TG0;
    ts1h.P = Pref;
    ts2h.P = Pref;
    ts1h.x = liq_start;
    ts2h.x = gas_start;

    rho1v[] = rho1;
    rho1v0[] = rho1;
    rho1vInt[] = rho1;
    mu1v[] = mu1;
    cp1v[] = cp1;
    lambda1v[] = lambda1;

    foreach_elem (YLList, jj) {
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

    foreach_elem (Dmix2List, jj) {
      scalar Dmix2v = Dmix2List[jj];
      Dmix2v[] = inDmix2[jj];

      scalar Cp2v = Cp2List[jj];
      Cp2v[] = cp2;
    }
  }
  boundary ({rho1v,rho1v0,rho1vInt,mu1v,cp1v,lambda1v,
      rho2v, rho2v0, mu2v, cp2v, lambda2v});
  boundary (dhevList);
  boundary (Dmix1List);
  boundary (Cp1List);
  boundary (Dmix2List);
  boundary (Cp2List);
}

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

    foreach_elem (YLList, jj) {
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

    foreach_elem (Dmix2List, jj) {
      scalar Dmix2v = Dmix2List[jj];
      Dmix2v[] = tp2.diff (&ts2h, jj);

      scalar Cp2v = Cp2List[jj];
      Cp2v[] = tp2.cps (&ts2h, jj);
    }
  }
  boundary ({rho1v,rho1v0,rho1vInt,mu1v,cp1v,lambda1v,
      rho2v, rho2v0, mu2v, cp2v, lambda2v});
  boundary (dhevList);
  boundary (Dmix1List);
  boundary (Cp1List);
  boundary (Dmix2List);
  boundary (Cp2List);
}

extern double mLiq0;

event init (i = 0) {
  update_properties_initial();

  mLiq0 = 0.;
  foreach(reduction(+:mLiq0))
    mLiq0 += rho1v[]*f[]*dv();

  MW1mix.dirty = false;
  MW2mix.dirty = false;

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

void update_properties (void)
{
  double MW1[NLS], MW2[NGS];
  foreach_elem (YLList, jj)
    MW1[jj] = inMW[LSI[jj]];
  foreach_elem (YGList, jj)
    MW2[jj] = inMW[jj];

  foreach() {
    MW1mix[] = 0.;
    MW2mix[] = 0.;

    if (f[] > T_PROP) {
      double x1[NLS], y1[NLS];
      foreach_elem (YLList, jj) {
        scalar YL = YLList[jj];
        y1[jj] = (NLS == 1) ? 1. : YL[];
      }
      correctfrac (y1, NLS);
      mass2molefrac (x1, y1, MW1, NLS);
      MW1mix[] = mass2mw (y1, MW1, NLS);

      ThermoState ts1h;
      ts1h.T = TL[];
      ts1h.P = Pref;
      ts1h.x = x1;

      rho1v[] = tp1.rhov (&ts1h);
      mu1v[] = tp1.muv (&ts1h);
      cp1v[] = tp1.cpv (&ts1h);
      lambda1v[] = tp1.lambdav (&ts1h);
      betaexp1[] = liqprop_thermal_expansion (&tp1, &ts1h);

      foreach_elem (YLList, jj) {
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
    //else {
    //  //rho1v[] = 0.;
    //  rho1v0[] = 0.;
    //  mu1v[] = 0.;
    //  cp1v[] = 0.;
    //  cp1v[] = 0.;
    //  lambda1v[] = 0.;
    //  betaexp1[] = 0.;

    //  foreach_elem (Dmix1List, jj) {
    //    scalar Dmix1v = Dmix1List[jj];
    //    scalar dhevjj = dhevList[jj];
    //    scalar Cp1v   = Cp1List[jj];
    //    Dmix1v[] = 0.;
    //    dhevjj[] = 0.;
    //    Cp1v[]   = 0.;
    //  }
    //}

    if ((1. - f[]) > T_PROP) {
      double x2[NGS], y2[NGS];
      foreach_elem (YGList, jj) {
        scalar YG = YGList[jj];
        y2[jj] = YG[];
      }
      correctfrac (y2, NGS);
      mass2molefrac (x2, y2, MW2, NGS);
      MW2mix[] = mass2mw (y2, MW2, NGS);

      ThermoState ts2h;
      ts2h.T = TG[];
      ts2h.P = Pref;
      ts2h.x = x2;

      rho2v[] = tp2.rhov (&ts2h);
      mu2v[] = tp2.muv (&ts2h);
      cp2v[] = tp2.cpv (&ts2h);
      lambda2v[] = tp2.lambdav (&ts2h);
      betaexp2[] = gasprop_thermal_expansion (&ts2h);

      foreach_elem (Dmix2List, jj) {
        scalar Dmix2v = Dmix2List[jj];
        Dmix2v[] = tp2.diff (&ts2h, jj);

        scalar Cp2v = Cp2List[jj];
        Cp2v[] = tp2.cps (&ts2h, jj);
      }
    }
    //else {
    //  //rho2v[] = 0.;
    //  rho2v0[] = 0.;
    //  mu2v[] = 0.;
    //  cp2v[] = 0.;
    //  betaexp2[] = 0.;

    //  foreach_elem (Dmix2List, jj) {
    //    scalar Dmix2v = Dmix2List[jj];
    //    scalar Cp2v   = Cp2List[jj];
    //    Dmix2v[] = 0.;
    //    Cp2v[] = 0.;
    //  }
    //}
  }

  //// Update interface properties
  //foreach() {
  //  if (f[] > F_ERR && f[] < 1.-F_ERR) {

  //    // Liquid interface side
  //    double x1[NLS], y1[NLS];
  //    foreach_elem (YLIntList, jj) {
  //      scalar YLInt = YLIntList[jj];
  //      y1[jj] = YLInt[];
  //    }
  //    correctfrac (y1, NLS);
  //    mass2molefrac (x1, y1, MW1, NLS);

  //    ThermoState ts1h;
  //    ts1h.T = TInt[];
  //    ts1h.P = Pref;
  //    ts1h.x = x1;

  //    //rho1v[] = tp1.rhov (&ts1h);
  //    //mu1v[] = tp1.muv (&ts1h);
  //    //cp1v[] = tp1.cpv (&ts1h);
  //    //lambda1v[] = tp1.lambdav (&ts1h);
  //    //betaexp1[] = liqprop_thermal_expansion (&tp1, &ts1h);

  //    //foreach_elem (YLIntList, jj) {
  //    //  scalar dhevjj = dhevList[jj];
  //    //  dhevjj[] = tp1.dhev (&ts1h, jj);

  //    //  scalar Dmix1v = Dmix1List[jj];
  //    //  Dmix1v[] = tp1.diff (&ts1h, jj);
  //    //}

  //    // Gas interface side
  //    double x2[NGS], y2[NGS];
  //    foreach_elem (YGIntList, jj) {
  //      scalar YGInt = YGIntList[jj];
  //      y2[jj] = YGInt[];
  //    }
  //    correctfrac (y2, NGS);
  //    mass2molefrac (x2, y2, MW2, NGS);

  //    ThermoState ts2h;
  //    ts2h.T = TInt[];
  //    ts2h.P = Pref;
  //    ts2h.x = x2;

  //    rho2vInt[] = tp2.rhov (&ts2h);
  //    //mu2v[] = tp2.muv (&ts2h);
  //    //cp2v[] = tp2.cpv (&ts2h);
  //    //lambda2v[] = tp2.lambdav (&ts2h);
  //    //betaexp2[] = gasprop_thermal_expansion (&ts2h);

  //    //foreach_elem (YGIntList, jj) {
  //    //  scalar Dmix2v = Dmix2List[jj];
  //    //  Dmix2v[] = tp2.diff (&ts2h, jj);
  //    //}
  //  }
  //}

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
      rho1v0[] = rho1v[];
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
  }
  boundary ({rho1v,rho1v0,rho1vInt,mu1v,cp1v,lambda1v,
      rho2v, rho2v0, mu2v, cp2v, lambda2v});
  boundary (dhevList);
  boundary (Dmix1List);
  boundary (Cp1List);
  boundary (Dmix2List);
  boundary (Cp2List);

  for (int b = 0; b < nboundary; b++) {
    foreach_boundary(b) {

      // liquid phase
      double x1[NLS], y1[NLS];
      foreach_elem (YLList, jj) {
        scalar YL = YLList[jj];
        double YLf = 0.5*(YL[] + get_ghost (point, YL, b));
        y1[jj] = (NLS == 1.) ? 1. : YLf;
      }
      correctfrac (y1, NLS);
      mass2molefrac (x1, y1, MW1, NLS);
      double MW1mixf = mass2mw (y1, MW1, NLS);
      set_ghost (point, MW1mix, b, MW1mixf);

      // gas phase
      double x2[NGS], y2[NLS];
      foreach_elem (YGList, jj) {
        scalar YG = YGList[jj];
        double YGf = 0.5*(YG[] + get_ghost (point, YG, b));
        y2[jj] = YGf;
      }
      correctfrac (y2, NGS);
      mass2molefrac (x2, y2, MW2, NGS);
      double MW2mixf = mass2mw (y2, MW2, NGS);
      set_ghost (point, MW2mix, b, MW2mixf);
    }
  }
}

void update_divergence (void) {

  /**
  We define the variables used to compute the lagrangian derivative
  on each level. */

  restriction ({T,TL,TG});
  restriction (YLList);
  restriction (YGList);

  foreach() {

#ifdef CHEMISTRY
    double Qr = 0.;
    double massfracs[NGS], molefracs[NGS], ci[NGS], ri[NGS];
    foreach_elem (YGList, jj) {
      massfracs[jj] = 0.;
      molefracs[jj] = 0.;
      ci[jj] = 0.;
      ri[jj] = 0.;
    }

    // Set up chemical reactions contribution
    if ((1. - f[]) < F_ERR) {
      OpenSMOKE_GasProp_SetTemperature (TG[]);
      OpenSMOKE_GasProp_SetPressure (Pref);

      foreach_elem (YGList, jj) {
        scalar YG = YGList[jj];
        massfracs[jj] = YG[];
      }
      correctfrac (massfracs, NGS);
      mass2molefrac (molefracs, massfracs, inMW, NGS);

      double ctot = (TG[] > 0.) ? Pref/(R_GAS*1000.)/TG[] : 0.;
      foreach_elem (YGList, jj) {
        ci[jj] = ctot*molefracs[jj];
        ci[jj] = (ci[jj] < 0.) ? 0. : ci[jj];
        ri[jj] = 0.;
      }

      OpenSMOKE_GasProp_KineticConstants();
      OpenSMOKE_GasProp_ReactionRates (ci);     // [kmol/m3]
      OpenSMOKE_GasProp_FormationRates (ri);    // [kmol/m3/s]

      Qr = OpenSMOKE_GasProp_HeatRelease (ri);
    }
#endif

    // Compute chemical species contributions
    double laplYtot = 0.;
    foreach_elem (YGList, jj) {
      scalar YG = YGList[jj];
      scalar Dmix2v = Dmix2List[jj];
      double laplYjj = 0.;
      foreach_dimension() {
        double Dmixfr = 0.5*(Dmix2v[1] + Dmix2v[]);
        double Dmixfl = 0.5*(Dmix2v[] + Dmix2v[-1]);
        double rhofr  = 0.5*(rho2v[1] + rho2v[]);
        double rhofl  = 0.5*(rho2v[] + rho2v[-1]);
        laplYjj += (fm.x[1]*fsG.x[]*rhofr*Dmixfr*face_gradient_x (YG, 1) -
            fm.x[]*fsG.x[]*rhofl*Dmixfl*face_gradient_x (YG, 0));
      }
      laplYjj /= Delta;

      // Add interfacial contribution
      scalar sgexp = sgexpList[jj];
      scalar sgimp = sgimpList[jj];
      laplYjj += sgexp[];
      laplYjj += sgimp[]*YG[];

      // Add chemical reactions contribution
#ifdef CHEMISTRY
        //laplYjj += OpenSMOKE_MW(jj)*ri[jj]*cm[]*(1. - f[]);
        laplYjj += OpenSMOKE_MW(jj)*ri[jj]*(1. - f[]);
#endif

      // Multiply by the species molecular weight
      //laplYtot += 1./(inMW[jj]*Delta)*laplYjj;
      laplYtot += 1./inMW[jj]*laplYjj;
    }
    laplYtot *= MW2mix[]/rho2v[];

    // Compute temperature contribution
    double laplT1 = 0.;
    foreach_dimension() {
      double lambdafr = 0.5*(lambda1v[1] + lambda1v[]);
      double lambdafl = 0.5*(lambda1v[] + lambda1v[-1]);
      laplT1 += (fm.x[1]*fsL.x[1]*lambdafr*face_gradient_x (TL, 1) -
          fm.x[]*fsL.x[]*lambdafl*face_gradient_x (TL, 0));
    }
    laplT1 /= Delta;

    double laplT2 = 0.;
    foreach_dimension() {
      double lambdafr = 0.5*(lambda2v[1] + lambda2v[]);
      double lambdafl = 0.5*(lambda2v[] + lambda2v[-1]);
      laplT2 += (fm.x[1]*fsG.x[1]*lambdafr*face_gradient_x (TG, 1) -
          fm.x[]*fsG.x[]*lambdafl*face_gradient_x (TG, 0));
    }
    laplT2 /= Delta;

#ifdef CHEMISTRY
      //laplT2 += Qr*cm[]*(1. - f[]);
      laplT2 += Qr*(1. - f[]);
#endif

//    // Add temperature interface contribution
//    if (f[] > F_ERR && f[] < 1.-F_ERR) {
//      coord n = facet_normal (point, f, fsL), p;
//      double alpha = plane_alpha (f[], n);
//      double area = plane_area_center (n, alpha, &p);
//      normalize (&n);
//
//      double bc = TInt[];
//      double gtrgrad = ebmgrad (point, TG, fL, fG, fsL, fsG, true, bc, &success);
//      double ltrgrad = ebmgrad (point, TL, fL, fG, fsL, fsG, false, bc, &success);
//
//      double lheatflux = lambda1v[]*ltrgrad;
//      double gheatflux = lambda2v[]*gtrgrad;
//
//#ifdef AXI
//    laplT1 += lheatflux*area*(y + p.y*Delta)/(Delta*y)*cm[];
//    laplT2 += gheatflux*area*(y + p.y*Delta)/(Delta*y)*cm[];
//#else
//    laplT1 += lheatflux*area/Delta*cm[];
//    laplT2 += gheatflux*area/Delta*cm[];
//#endif
//    }
    laplT1 += slT[];
    laplT2 += sgT[];

    double drho1dt = 0.;
    double drho2dt = 0.;

    // Add liquid compressibility due to temperature
    drho1dt += (rho1v[]*cp1v[] > 0.) ?
      -betaexp1[]/(rho1v[]*cp1v[])*laplT1 : 0.;

    drho2dt += (TG[]*rho2v[]*cp2v[] > 0.) ?
      -1./(TG[]*rho2v[]*cp2v[])*laplT2 : 0.;

    // Add gas compressibility due to composition
    //drho2dt += ((1. - f[]) > F_ERR) ? -laplYtot : 0.;
    drho2dt += (f[] == 0.) ? -laplYtot : 0.;

    drhodt[] = drho1dt*f[] + drho2dt*(1. - f[]);
    drhodtext[] = drho1dt;
  }
  boundary ({drhodt, drhodtext});
}

