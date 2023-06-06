
#include "curvature.h"
#include "intgrad.h"
#include "fsolve-gsl.h"
#include "thermodynamics.h"

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

#ifdef SOLVE_TEMPERATURE

void EqTemperature (const double * x, double * f, size_t size, void * params) {

  Point point = data->point;

  double TInti = x[0];
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

  f[0] = vapheat
       + lambda1*ltrgrad
       + lambda2*gtrgrad
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
#ifdef USE_GSL
      fsolve (EqTemperatureGsl, arrUnk, point);
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

