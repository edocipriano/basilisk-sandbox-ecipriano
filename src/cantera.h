#include "cantera/clib/ct.h"
#include "cantera/clib/ctreactor.h"
#include "cantera/clib/ctfunc.h"
#include "cantera/clib/ctmultiphase.h"
#include "cantera/clib/ctonedim.h"
#include "cantera/clib/ctrpath.h"
#include "cantera/clib/ctsurf.h"

#define CANTERA 1

#pragma autolink -L$CANTERA_LIBRARY_PATH -lcantera

int soln = -1;
int thermo = -1, kin = -1, tran = -1;

void kinetics (char * path, bool liquid = false) {
  soln = soln_newSolution (path, "gas", "default");
  thermo = soln_thermo (soln);
  kin = soln_kinetics (soln);
  tran = soln_transport (soln);
  size_t nr = kin_nReactions (soln);
}

void delete_kinetics (void) {
  ct_appdelete();
}

