
#include "intgrad.h"
#include "fsolve-gsl.h"
#include "thermodynamics.h"

extern scalar * mEvapList;
extern scalar * YLList;
extern scalar * YGList;
extern scalar * YLIntList;
extern scalar * YGIntList;
#ifdef VARPROP
extern scalar * Dmix1List;
extern scalar * Dmix2List;
extern scalar * dhevList;
#endif
extern double inDmix1[NLS];
extern double inDmix2[NGS];
extern double inKeq[NLS];
extern double inMW[NGS];
extern double Pref;
extern scalar fL, fG;
extern face vector fsL, fsG;
extern int * LSI;
extern int * GOSI;
extern int inertIndex;
extern int NGOS;
extern double rho1, rho2;
#ifdef SOLVE_TEMPERATURE
extern double lambda1, lambda2, dhev, cp1, cp2, TL0, TG0;
extern scalar TInt, TL, TG, T;
#endif

#ifdef CLAPEYRON
extern double Tsat;
#endif
#ifdef USE_CLAPEYRON
extern double Tboil[NLS];
#endif

typedef struct {
  Point point;
} UserDataNls;

void Equations (const double * xdata, double * fdata, void * params) {
  UserDataNls * data = (UserDataNls *)params;
  Point point = data->point;

  double mEvapToti[NGS];
  double mEvapi[NLS];
  double YLInti[NLS];
  double YGInti[NGS];
  double XLInti[NLS];
  double XGInti[NGS];
  double jL[NLS];
  double jG[NGS];
  double mEvapSum;
  bool success;
  int count;
#ifdef SOLVE_TEMPERATURE
  double TInti;
  double gradTGn;
  double gradTLn;
#endif

  /**
  Rename unknowns for convenience. */

  count = 0;

  for (int jj=0; jj<NLS; jj++)
    mEvapi[jj] = xdata[count++];

  for (int jj=0; jj<NLS; jj++)
    YLInti[jj] = xdata[count++];

  for (int jj=0; jj<NGS; jj++)
    YGInti[jj] = xdata[count++];

#ifdef SOLVE_TEMPERATURE
  TInti = xdata[count++];
#endif

  /**
  Set to zero gas-only species vaporization mass flow-rates. */

  for (int jj=0; jj<NLS; jj++)
    mEvapToti[LSI[jj]] = mEvapi[jj];

  for (int jj=0; jj<NGOS; jj++)
    mEvapToti[GOSI[jj]] = 0.;

  /**
  Compute diffusive fluxes. */

  for (int jj=0; jj<NGS; jj++) {
    scalar YG = YGList[jj];
    double gtrgrad = ebmgrad (point, YG, fL, fG, fsL, fsG, true, YGInti[jj], &success);
    jG[jj] = -rho2*inDmix2[jj]*gtrgrad;
  }
  for (int jj=0; jj<NLS; jj++) {
    scalar YL = YLList[jj];
    double ltrgrad = ebmgrad (point, YL, fL, fG, fsL, fsG, false, YLInti[jj], &success);
    jL[jj] = rho1*inDmix1[jj]*ltrgrad;
  }

#ifdef SOLVE_TEMPERATURE
  gradTGn = ebmgrad (point, TG, fL, fG, fsL, fsG, true,  TInti, &success);
  gradTLn = ebmgrad (point, TL, fL, fG, fsL, fsG, false, TInti, &success);
#endif

  /**
  Convert mass-to-mole fractions. */

  {
    double inMWG[NGS];
    for (int jj=0; jj<NGS; jj++)
      inMWG[jj] = inMW[jj];

    double inMWL[NLS];
    for (int jj=0; jj<NLS; jj++)
      inMWL[jj] = inMW[LSI[jj]];

    mass2molefrac (XLInti, YLInti, inMWL, NLS);
    mass2molefrac (XGInti, YGInti, inMWG, NGS);
  }

  /**
  Compute sum of vaporization mass flow-rate. */

  mEvapSum = 0.;
  for (int jj=0; jj<NLS; jj++) {
    mEvapSum += mEvapi[jj];
  }

  //double sum_YGInt = 0.;
  //double sum_jG = 0.;
  //for (int jj=0; jj<NLS; jj++) {
  //  sum_jG += jG[LSI[jj]];
  //  sum_YGInt += YGInti[LSI[jj]];
  //}
  //mEvapSum = sum_jG/(1. - sum_YGInt);

  count = 0;

  /**
  Mass balance for species in liquid phase. */

  //for (int jj=0; jj<NLS-1; jj++) {
  ////for (int jj=0; jj<NLS; jj++) {

  //  fdata[count++] = mEvapi[jj]
  //                 - mEvapSum*YLInti[jj]
  //                 - jL[jj]
  //                 ;
  //}

  //double sum_YLInt = 0.;
  //for (int jj=0; jj<NLS; jj++)
  //  sum_YLInt += YLInti[jj];
  //fdata[count++] = 1. - sum_YLInt;

  // TODO: only for pure liquid
  if (NLS > 1) {
    for (int jj=0; jj<NLS; jj++) {

      fdata[count++] = mEvapi[jj]
                     - mEvapSum*YLInti[jj]
                     - jL[jj]
                     ;
    }
  }
  else {
    fdata[count++] = YLInti[0] - 1.;
  }

  /**
  Mass balance for species in gas phase and for gas-only species. */

  for (int jj=0; jj<NGS; jj++) {

    fdata[count++] = mEvapToti[jj]
                   - mEvapSum*YGInti[jj]
                   - jG[jj]
                   ;
  }
  /*
     OK per 1 specie
  for (int jj=0; jj<NLS; jj++) {

    fdata[count++] = mEvapi[jj]
                   - mEvapSum*YGInti[LSI[jj]]
                   - jG[LSI[jj]]
                   ;
  }

  fdata[count++] = 1. - sum_YGInt - YGInti[inertIndex];
  */

  /**
  Thermodynamic (VLE) equilibrium at the interface. */

#ifdef CLAPEYRON
  for (int jj=0; jj<NLS; jj++)
    inKeq[jj] = exp(-dhev/8.314*inMW[LSI[jj]]*(1./TInti - 1./Tsat))*YLInti[jj];
#endif

  double Keq[NLS];
  for (int jj=0; jj<NLS; jj++)
    Keq[jj] = inKeq[jj];

#ifdef USE_CLAPEYRON
  for (int jj=0; jj<NLS; jj++) {
    Keq[jj] = min (clapeyron ( min (TInti, Tboil[jj]-1), Tboil[jj], dhev, inMW[LSI[jj]]), 0.98);
  }
#endif
#ifdef USE_ANTOINE
  for (int jj=0; jj<NLS; jj++) {
    scalar YL = YLList[jj];
    Keq[jj] = min (YL.antoine (TInti, Pref), 0.98);
  }
#endif

  for (int jj=0; jj<NLS; jj++) {
    fdata[count++] = XLInti[jj]*Keq[jj] - XGInti[LSI[jj]];
  }

#ifdef SOLVE_TEMPERATURE

  /**
  Interface energy balance. */

  double vapheat = 0.;
  for (int jj=0; jj<NLS; jj++)
    vapheat -= dhev*mEvapi[jj];

  fdata[count++] = vapheat
                 + lambda1*gradTLn
                 + lambda2*gradTGn
                 ;

#endif
}

int EquationsGsl (const gsl_vector * x, void * params, gsl_vector * f) {
  double * xdata = x->data;
  double * fdata = f->data;

  Equations (xdata, fdata, params);
  return GSL_SUCCESS;
}

/**
A partire da *YLIntList*, *YGIntList*, *mEvapList* e delle proprietà
fisiche come densità e diffusività, bisogna calcolare la jump condition. */

void ijc_CoupledNls ()
{
  foreach() {
    if (f[] > F_ERR && f[] < 1.-F_ERR) {
      /**
      Reset variables. */

      Array * arrUnk = array_new();

      // mInt field
      for (int jj=0; jj<NLS; jj++) {
        scalar mEvapi = mEvapList[LSI[jj]];
        double vali = mEvapi[];
        array_append (arrUnk, &vali, sizeof(double));
      }

      // YLInt field
      for (int jj=0; jj<NLS; jj++) {
        //scalar YLi = YLList[jj];
        //double vali = YLi[]/max(f[], I_TOL);
        scalar YLInti = YLIntList[jj];
        double vali = YLInti[];
        array_append (arrUnk, &vali, sizeof(double));
      }

      // YGInt field
      for (int jj=0; jj<NGS; jj++) {
        //scalar YGi = YGList[jj];
        //double vali = YGi[]/max(1. - f[], I_TOL);
        scalar YGInti = YGIntList[jj];
        double vali = YGInti[];
        array_append (arrUnk, &vali, sizeof(double));
      }

#ifdef SOLVE_TEMPERATURE
      // Temperature field
      {
        double vali = TInt[];
        array_append (arrUnk, &vali, sizeof(double));
      }
#endif

      /**
      Gather parameters that must be passed to the
      non-linear system function. */

      UserDataNls data;
      data.point = point;

      /**
      Solve the non-linear system of equations. */

      fsolve (EquationsGsl, arrUnk, &data);

      /**
      Recover Nls solution. */
      {
        double * unk = (double *)arrUnk->p;
        //double * unk = (double *) array_shrink (arrUnk);
        int count = 0;

        // mInt field
        for (int jj=0; jj<NLS; jj++) {
          scalar mEvapi = mEvapList[LSI[jj]];
          mEvapi[] = unk[count++];
        }

        // YLInt field
        for (int jj=0; jj<NLS; jj++) {
          scalar YLInti = YLIntList[jj];
          YLInti[] = unk[count++];
        }

        // YGInt field
        for (int jj=0; jj<NGS; jj++) {
          scalar YGInti = YGIntList[jj];
          YGInti[] = unk[count++];
        }
#ifdef SOLVE_TEMPERATURE
        // Temperature field
        {
          TInt[] = unk[count++];
        }
#endif
      }
      array_free (arrUnk);
    }
  }
}

#ifdef SOLVE_TEMPERATURE

void EqTemperature (const double * xdata, double * fdata, void * params) {
  UserDataNls * data = (UserDataNls *)params;
  Point point = data->point;

  double TInti = xdata[0];
  bool success = false;

  double gtrgrad = ebmgrad (point, TG, fL, fG, fsL, fsG, true,  TInti, &success);
  double ltrgrad = ebmgrad (point, TL, fL, fG, fsL, fsG, false, TInti, &success);

  double vapheat = 0.;
  double mEvapi[NLS];
  for (int jj=0; jj<NLS; jj++) {
    scalar mEvapSi = mEvapList[LSI[jj]];
    mEvapi[jj] = mEvapSi[];
    double dhevvh = dhev;
#ifdef VARPROP
    scalar dhevjj = dhevList[jj];
    dhevvh = dhevjj[jj];
#endif
    vapheat -= mEvapi[jj]*dhevvh;
  }

  fdata[0] = vapheat
#ifdef VARPROP
      + lambda1v[]*ltrgrad
      + lambda2v[]*gtrgrad
#else
       + lambda1*ltrgrad
       + lambda2*gtrgrad
#endif
       ;
}

int EqTemperatureGsl (const gsl_vector * x, void * params, gsl_vector * f) {
  double * xdata = x->data;
  double * fdata = f->data;

  EqTemperature (xdata, fdata, params);
  return GSL_SUCCESS;
}


void ijc_CoupledTemperature ()
{
  foreach() {
    if (f[] > F_ERR && f[] < 1.-F_ERR) {
      Array * arrUnk = array_new();
      {
        double vali = TInt[];
        array_append (arrUnk, &vali, sizeof(double));
      }
      UserDataNls data;
      data.point = point;
#ifdef USE_GSL
      fsolve (EqTemperatureGsl, arrUnk, &data);
#endif
      {
        double * unk = (double *)arrUnk->p;
        TInt[] = unk[0];
      }
      array_free (arrUnk);
    }
  }
}

#endif // SOLVE_TEMPERATURE

