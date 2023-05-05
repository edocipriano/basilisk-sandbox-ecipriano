/**
# Non-Linear System of Equations Solver

This module defines a function *fsolve()* which, in analogy with
the MATLAB function, provides a high-level interface for the solution
of non-linear systems of equations.
*/

/**
## KINSol Interface

if the [SUNDIALS library](https://github.com/LLNL/sundials) is used,
the function *fsolve()* rely on the KINSol solver.
*/

#ifdef USE_SUNDIALS

#include <kinsol/kinsol.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sundials/sundials_types.h>

#define KIN_FTOL   1.e-6 // function tolerance 
#define KIN_STOL   1.e-6 // step tolerance

#pragma autolink -L$SUNDIALS_LIBRARY_PATH/build/src/kinsol -lsundials_kinsol

typedef struct {
  int neq;
  Point point;
} *UserData;

struct Fsolve {
  KINSysFn func;
  Array * arrUnk;
  Point point;
};

//void fsolve (struct Fsolve p)
void fsolve (KINSysFn func, Array * arrUnk, Point point)
{
  //KINSysFn func = p.func;
  //Array * arrUnk = p.arrUnk;
  //Point point = p.point;
  int neq = arrUnk->len/sizeof(double);

  UserData data;
  N_Vector u, s;
  SUNMatrix J;
  SUNLinearSolver LS;

  u    = NULL;
  s    = NULL;
  J    = NULL;
  LS   = NULL;
  data = NULL;

  data = (UserData)malloc(sizeof *data);
  data->neq = neq;
  data->point = point;

  u = N_VNew_Serial (neq);
  s = N_VNew_Serial (neq);

  {
    realtype * udata = N_VGetArrayPointer (u);
    double * unk = (double *)arrUnk->p;

    for (unsigned int jj=0; jj<neq; jj++)
      udata[jj] = unk[jj];
  }

  N_VConst (1.0, s);

  void * kmem;
  kmem = KINCreate();
  KINSetUserData (kmem, data);
  KINSetFuncNormTol (kmem, KIN_FTOL);
  KINSetScaledStepTol (kmem, KIN_STOL);
  KINInit (kmem, func, u);
  J = SUNDenseMatrix (neq, neq);
  LS = SUNLinSol_Dense (u, J);
  KINSetLinearSolver (kmem, LS, J);
  KINSetMaxSetupCalls(kmem, 1);
  KINSetPrintLevel (kmem, 0);

  //TODO following two lines give problems
  //{
  //  char name[80];
  //  sprintf (name, "KINSolErr");
  //  const FILE * fKinErr = fopen (name, "a");
  //  KINSetErrFile (kmem, fKinErr);
  //}

  /**
  Solve non-linear system of equations. */

  KINSol (kmem, u, KIN_NONE, s, s);
  //KINSol (kmem, u, KIN_LINESEARCH, s, s);

  /**
  Recover Nls solution. */

  {
    realtype * udata = N_VGetArrayPointer (u);
    double * unk = (double *)arrUnk->p;
    for (unsigned jj=0; jj<neq; jj++) {
      unk[jj] = udata[jj];
    }
  }

  /**
  Free memory. */

  N_VDestroy (u);
  N_VDestroy (s);
  KINFree (&kmem);
  SUNLinSolFree (LS);
  SUNMatDestroy (J);
  free (data);
}

#endif
