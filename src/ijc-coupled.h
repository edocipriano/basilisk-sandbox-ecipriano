
#ifdef PCM_MULTICOMPONENT

#include "curvature.h"
#include "intgrad.h"
#include "fsolve.h"
#include "thermodynamics.h"

#ifdef USE_SUNDIALS

//#include <kinsol/kinsol.h>
//#include <nvector/nvector_serial.h>
//#include <sunmatrix/sunmatrix_dense.h>
//#include <sunlinsol/sunlinsol_dense.h>
//#include <sundials/sundials_types.h>
//
//#define KIN_FTOL   1.e-5 // function tolerance 
//#define KIN_STOL   1.e-5 // step tolerance
//#define KIN_FTOL   1.e-6 // function tolerance 
//#define KIN_STOL   1.e-6 // step tolerance

#pragma autolink -L$SUNDIALS_LIBRARY_PATH/build/src/kinsol -lsundials_kinsol

extern scalar * mEvapList;
extern scalar * YLList;
extern scalar * YGList;
extern scalar * YLIntList;
extern scalar * YGIntList;
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
  int neq;
  Point point;
} *UserDataNls;

static int Equations (N_Vector usol, N_Vector fsol, void *user_data)
{
  realtype * udata = NULL;
  realtype * fdata = NULL;

  udata = N_VGetArrayPointer(usol);
  fdata = N_VGetArrayPointer(fsol);

  UserDataNls data;
  data = (UserDataNls)user_data;
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
    mEvapi[jj] = udata[count++];

  for (int jj=0; jj<NLS; jj++)
    YLInti[jj] = udata[count++];

  for (int jj=0; jj<NGS; jj++)
    YGInti[jj] = udata[count++];

#ifdef SOLVE_TEMPERATURE
  TInti = udata[count++];
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

  double sum_jG = 0.;
  double sum_YGInt = 0.;
  for (int jj=0; jj<NLS; jj++) {
    sum_jG += jG[LSI[jj]];
    sum_YGInt += YGInti[LSI[jj]];
  }

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

  //double P_Pa = 101325.;

#ifdef THERMOPROP
  double dhevp, fugacity_G, fugacity_L;
  thermoprop_nHeptane (TInti, P_Pa, 1., 0.1, &dhevp, &fugacity_L, &fugacity_G);
  inKeq[0] = min(thermodata.Keq, 0.99);
  inKeq[1] = 0.;
  //fprintf (stdout, "Keq = %f\n", inKeq[0]); fflush (stdout);
#endif

#ifdef CLAPEYRON
  for (int jj=0; jj<NLS; jj++)
    inKeq[jj] = exp(-dhev/8.314*inMW[LSI[jj]]*(1./TInti - 1./Tsat))*YLInti[jj];
#endif

  double Keq[NLS];
  for (int jj=0; jj<NLS; jj++)
    Keq[jj] = inKeq[jj];

#ifdef USE_CLAPEYRON
  for (int jj=0; jj<NLS; jj++) {
    //inKeq[jj] = min (clapeyron ( min (TInti, Tboil[jj]-1), Tboil[jj], dhev, inMW[LSI[jj]]), 0.98);
    Keq[jj] = min (clapeyron ( min (TInti, Tboil[jj]-1), Tboil[jj], dhev, inMW[LSI[jj]]), 0.98);
  }
#endif
#ifdef USE_ANTOINE
  for (int jj=0; jj<NLS; jj++) {
    scalar YL = YLList[jj];
    Keq[jj] = min (YL.antoine (TInti, Pref), 0.99);
  }
#endif

  for (int jj=0; jj<NLS; jj++) {
    fdata[count++] = XLInti[jj]*Keq[jj] - XGInti[LSI[jj]];
    //fdata[count++] = (P_Pa*YLInti[jj]*inKeq[jj] - P_Pa*YGInti[LSI[jj]]);
    //fdata[count++] = (P_Pa*XLInti[jj]*inKeq[jj] - P_Pa*XGInti[LSI[jj]]);
    //fdata[count++] = fugacity_G - fugacity_L;
  }

#ifdef SOLVE_TEMPERATURE

  /**
  Interface energy balance. */

  double vapheat = 0.;
  for (int jj=0; jj<NLS; jj++)
    //vapheat -= mEvapi[jj]*(dhev + (cp1 - cp2)*(Tboil[jj] - TInti));
    vapheat -= dhev*mEvapi[jj];

  fdata[count++] = vapheat
                 + lambda1*gradTLn
                 + lambda2*gradTGn
                 ;

#endif

  return 0;
}
/*
static int Equations2 (N_Vector usol, N_Vector fsol, void *user_data)
{
  realtype * udata = NULL;
  realtype * fdata = NULL;

  udata = N_VGetArrayPointer(usol);
  fdata = N_VGetArrayPointer(fsol);

  UserDataNls data;
  data = (UserDataNls)user_data;
  Point point = data->point;

  double mEvapToti[NGS];
  double mEvapi[NLS];
  double YLInti[NLS];
  double YGInti[NGS];
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
  */
  /**
  Rename unknowns for convenience. */
  /*
  count = 0;

  for (int jj=0; jj<NLS; jj++)
    mEvapi[jj] = udata[count++];

  for (int jj=0; jj<NLS; jj++)
    YLInti[jj] = udata[count++];

  for (int jj=0; jj<NGS; jj++)
    YGInti[jj] = udata[count++];

#ifdef SOLVE_TEMPERATURE
  TInti = udata[count++];
#endif
  */
  /**
  Set to zero gas-only species vaporization mass flow-rates. */
  /*
  for (int jj=0; jj<NLS; jj++)
    mEvapToti[LSI[jj]] = mEvapi[jj];

  for (int jj=0; jj<NGOS; jj++)
    mEvapToti[GOSI[jj]] = 0.;
  */
  /**
  Compute diffusive fluxes. */
  /*
  for (int jj=0; jj<NGS; jj++) {
    scalar YG = YGList[jj];
    double gtrgrad = embedgrad (point, YG, fL, fG, 2, YGInti[jj], &success);
    jG[jj] = -rho2*inDmix2[jj]*gtrgrad;
  }
  for (int jj=0; jj<NLS; jj++) {
    scalar YL = YLList[jj];
    double ltrgrad = embedgrad (point, YL, fL, fG, 1, YLInti[jj], &success);
    jL[jj] = rho1*inDmix1[jj]*ltrgrad;
  }

#ifdef SOLVE_TEMPERATURE
  gradTGn = embedgrad (point, TG, fL, fG, 2, TInti, &success);
  gradTLn = embedgrad (point, TL, fL, fG, 1, TInti, &success);
#endif
  */
  /**
  Compute sum of vaporization mass flow-rate. */
  /*
  mEvapSum = 0.;
  for (int jj=0; jj<NLS; jj++) {
    mEvapSum += mEvapi[jj];
  }

  double sum_jG = 0.;
  double sum_YGInt = 0.;
  for (int jj=0; jj<NLS; jj++) {
    sum_jG += jG[LSI[jj]];
    sum_YGInt += YGInti[LSI[jj]];
  }

  //mEvapSum = sum_jG/(1. - sum_YGInt);

  count = 0;
  */
  /**
  [NLS] vaporization rates (from mass balance on gas phase). */
  /*
  for (int jj=0; jj<NLS; jj++) {

    fdata[count++] = mEvapi[jj]
                   - mEvapSum*YGInti[LSI[jj]]
                   - jG[LSI[jj]]
                   ;
  }
  */
  /**
  [NLS] liquid interface mass fractions (from mass balance on liq phase). */
  /*
  for (int jj=0; jj<NLS; jj++) {

    fdata[count++] = mEvapi[jj]
                   - mEvapSum*YLInti[jj]
                   - jL[jj]
                   ;
  }
  */
  /**
  [NLS] gas interface mass fractions (from VLE). */
  /*
  for (int jj=0; jj<NLS; jj++) {

    fdata[count++] = (YLInti[jj]*inKeq[jj] - YGInti[LSI[jj]]);
  }
  */
  /**
  [NGOS] gas-only interface mass fractions (from mass balance on gas-only species). */
  /*
  for (int jj=0; jj<NGOS; jj++) {

    fdata[count++] = 0.
                   - mEvapSum*YGInti[GOSI[jj]]
                   - jG[GOSI[jj]]
                   ;
  }

#ifdef SOLVE_TEMPERATURE
  */
  /**
  Interface energy balance. */
  /*
  double vapheat = 0.;
  for (int jj=0; jj<NLS; jj++)
    vapheat -= dhev*mEvapi[jj];

  fdata[count++] = vapheat
                 + lambda1*gradTLn
                 + lambda2*gradTGn
                 ;

#endif

  return 0;
}
*/

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
      Solve the non-linear system of equations. */

      fsolve (Equations, arrUnk, point);

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

//void ijc_CoupledNls2 ()
//{
//
//  foreach() {
//
//    if (is_interfacial(point, f)) {
//
//      /**
//      Reset variables. */
//
//      Array * arrUnk = array_new();
//
//      // mInt field
//      for (int jj=0; jj<NLS; jj++) {
//        scalar mEvapi = mEvapList[LSI[jj]];
//        double vali = mEvapi[];
//        array_append (arrUnk, &vali, sizeof(double));
//      }
//
//      // YLInt field
//      for (int jj=0; jj<NLS; jj++) {
//        //scalar YLi = YLList[jj];
//        //double vali = YLi[]/max(f[], I_TOL);
//        scalar YLInti = YLIntList[jj];
//        double vali = YLInti[];
//        array_append (arrUnk, &vali, sizeof(double));
//      }
//
//      // YGInt field
//      for (int jj=0; jj<NGS; jj++) {
//        //scalar YGi = YGList[jj];
//        //double vali = YGi[]/max(1. - f[], I_TOL);
//        scalar YGInti = YGIntList[jj];
//        double vali = YGInti[];
//        array_append (arrUnk, &vali, sizeof(double));
//      }
//
//#ifdef SOLVE_TEMPERATURE
//      // Temperature field
//      {
//        double vali = TInt[];
//        array_append (arrUnk, &vali, sizeof(double));
//      }
//#endif
//
//      int neq = arrUnk->len/sizeof(double);
//#ifndef SOLVE_TEMPERATURE
//      assert (neq == (NLS+NLS+NGS));
//#else
//      assert (neq == (NLS+NLS+NGS+1));
//#endif
//
//      /**
//      Allocate variables for KINSol solution. */
//
//      UserDataNls data;
//      int retval = 0;
//      int flag = 0;
//      N_Vector u, s;
//      SUNMatrix J;
//      SUNLinearSolver LS;
//
//      u    = NULL;
//      s    = NULL;
//      J    = NULL;
//      LS   = NULL;
//      data = NULL;
//
//      /**
//      Assign UserData variables. */
//
//      data = (UserDataNls)malloc(sizeof *data);
//      data->neq = neq;
//      data->point = point;
//
//      /**
//      Vectors for solution and scaling. */
//
//      u  = N_VNew_Serial (neq);
//      s  = N_VNew_Serial (neq);
//
//      /**
//      Set first guess values. */
//
//      {
//        realtype * udata = N_VGetArrayPointer (u);
//        double * unk = (double *)arrUnk->p;
//
//        for (unsigned int jj=0; jj<neq; jj++) {
//          udata[jj] = unk[jj];
//        }
//      }
//
//      /**
//      Set scaling value to 1.0. */
//
//      N_VConst (1.0, s);
//
//      void * kmem;
//      kmem = KINCreate ();
//      KINSetUserData(kmem, data);
//      KINSetFuncNormTol(kmem, KIN_FTOL);
//      KINSetScaledStepTol(kmem, KIN_STOL);
//      KINInit (kmem, Equations, u);
//      J = SUNDenseMatrix (neq, neq);
//      LS = SUNLinSol_Dense (u, J);
//      KINSetLinearSolver (kmem, LS, J);
//      KINSetMaxSetupCalls(kmem, 1);
//      //KINSetNumMaxIters(kmem,1e+5);
//      //KINSetScaledStepTol();
//      KINSetPrintLevel(kmem, 0);
//
//      /**
//      Solve non-linear system of equations. */
//
//      retval = KINSol(kmem,           // KINSol memory block
//                      u,              // initial guess on input; solution vector
//                      KIN_NONE,     // global strategy choice
//                      //KIN_LINESEARCH, // global strategy choice
//                      s,              // scaling vector, for the variable cc
//                      s);            // scaling vector for function values fval
//
//      /**
//      Recover Nls solution. */
//
//      {
//        realtype * udata = N_VGetArrayPointer (u);
//        double * unk = (double *)arrUnk->p;
//        for (unsigned jj=0; jj<neq; jj++) {
//          unk[jj] = udata[jj];
//        }
//
//        int count = 0;
//
//        // mInt field
//        for (int jj=0; jj<NLS; jj++) {
//          scalar mEvapi = mEvapList[LSI[jj]];
//          mEvapi[] = unk[count++];
//        }
//
//        // YLInt field
//        for (int jj=0; jj<NLS; jj++) {
//          scalar YLInti = YLIntList[jj];
//          YLInti[] = unk[count++];
//        }
//
//        // YGInt field
//        for (int jj=0; jj<NGS; jj++) {
//          scalar YGInti = YGIntList[jj];
//          YGInti[] = unk[count++];
//        }
//
//#ifdef SOLVE_TEMPERATURE
//        // Temperature field
//        {
//          TInt[] = unk[count++];
//        }
//#endif
//      }
//
//      /**
//      Destroy KINSol variables. */
//
//      /*
//      N_VDestroy (u);
//      N_VDestroy (s);
//      KINFree(&kmem);
//      SUNLinSolFree(LS);
//      SUNMatDestroy(J);
//      free(data);
//      */
//      /**
//      Destroy array with first-guess/solution. */
//
//      //array_free (arrUnk);
//    }
//  }
//
//}

#ifdef SOLVE_TEMPERATURE

static int EqTemperature (N_Vector usol, N_Vector fsol, void *user_data)
{
  realtype * udata = NULL;
  realtype * fdata = NULL;

  udata = N_VGetArrayPointer(usol);
  fdata = N_VGetArrayPointer(fsol);

  UserData data;
  data = (UserData)user_data;
  Point point = data->point;

  double TInti = udata[0];
  bool success = false;

  double gtrgrad = ebmgrad (point, TG, fL, fG, fsL, fsG, true,  TInti, &success);
  double ltrgrad = ebmgrad (point, TL, fL, fG, fsL, fsG, false, TInti, &success);

  double vapheat = 0.;
  double mEvapi[NLS];
  for (int jj=0; jj<NLS; jj++) {
    scalar mEvapSi = mEvapList[LSI[jj]];
    mEvapi[jj] = mEvapSi[];
    vapheat -= mEvapi[jj]*dhev;
    //vapheat -= mEvapi[jj]*(dhev + (cp1 - cp2)*(Tboil[jj] - TInti));
  }

  fdata[0] = vapheat
           + lambda1*ltrgrad
           + lambda2*gtrgrad
           ;

  return 0;
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
      fsolve (EqTemperature, arrUnk, point);
      {
        //double * unk = (double *) array_shrink (arrUnk);
        double * unk = (double *)arrUnk->p;
        TInt[] = unk[0];
        //free (unk);
      }
      array_free (arrUnk);
    }
  }
  boundary({TInt});
}

#endif // SOLVE_TEMPERATURE

#endif // USE_SUNDIALS
#endif // PCM_MULTICOMPONENT

