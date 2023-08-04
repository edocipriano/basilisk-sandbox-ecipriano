/**
# Mass Conservation Balances

This module implements mass consevation tests for the phase
change models. The operations implemented here are too boring
and long to be written in every simulation. Therefore, an external
header that automatically computes the balances is useful. The
implementation is based on the nomenclature and on the data
structures used in [multicomponent.h](multicomponent.h). However, it
stores everything in a separate structure, making this module general
also for other phase change models.

The mass conservation tests can be defined in different ways. If a
liquid phase is placed in a closed domain, and no phase change
occurs, the mass balance reads:
$$
  \int_V\rho_l dV = \text{constant}
$$
which means that the total liquid mass (assumed to be a function of
the volume fraction only) must be constant.
If phase change is present, the variation of the liquid mass depends
on the total amount of evaporated liquid:
$$
  \int_V \rho_l dV - \int_0^{\tau}dt\int_V\dot{m}\delta_\Gamma dV
  = \text{constant}
$$
where $\dot{m} < 0$ if the liquid evaporates. If the domain is not
closed and the liquid (or the gas-phase) is allowed to leave the
domain, the mass balance must include also the contribution of the
fluxes across the boundaries:
$$
  \int_V \rho_l dV +
  \int_0^{\tau}dt\oint_{S_b} \rho_l \mathbf{u}\cdot\mathbf{n} dS_b -
  \int_0^{\tau}dt\int_V\dot{m}\delta_\Gamma dV =
  \text{constant}
$$
Exactly the same concept applies for the gas phase. If the same
balances are written for every chemical species, the equation must
account for the mass fraction of the chemical species and the
diffusive flux:
$$
  \int_V \rho_l\omega_{i,l} dV +
  \int_0^{\tau}dt\oint_{S_b} \rho_l\omega_{i,j}
  \mathbf{u}\cdot\mathbf{n} dS_b -
  \int_0^{\tau}dt\oint_{S_b}
  \rho_l\mathcal{D}_{i,l}\nabla\omega_{i,j}\cdot\mathbf{n}dS_b -
  \int_0^{\tau}dt\int_V\dot{m}\delta_\Gamma dV =
  \text{constant}
$$

*/

#define BALANCES

extern int maxlevel;
FILE * fb;

/**
We define a structure that gathers data that are used or
modified in this module. */

struct MassBalance {
  // Data
  scalar * YLList;
  scalar * YGList;
  scalar * mEvapList;
  double * liq_start;
  double * gas_start;
  char ** liq_species;
  char ** gas_species;
  bool boundaries;
  // Results
  double totmass1;
  double totmass2;
  double totmass1init;
  double totmass2init;
  double evapmass;
  double * liq_mass;
  double * gas_mass;
  double * liq_massinit;
  double * gas_massinit;
  double * liq_massbd;
  double * gas_massbd;
  double * tot_evap;
  double * liq_relerr;
  double * gas_relerr;
};

struct MassBalance mb = {};

/**
We initialize the varibales of the mass balances struct. */

event defaults (i = 0) {
  mb.YLList = NULL;
  mb.YGList = NULL;
  mb.mEvapList = NULL;
  mb.liq_start = NULL;
  mb.gas_start = NULL;
  mb.totmass1 = 0.;
  mb.totmass2 = 0.;
  mb.totmass1init = 0.;
  mb.totmass2init = 0.;
  mb.evapmass = 0.;
  mb.liq_mass = NULL;
  mb.gas_mass = NULL;
  mb.liq_massinit = NULL;
  mb.gas_massinit = NULL;
  mb.liq_massbd = NULL;
  mb.gas_massbd = NULL;
  mb.tot_evap = NULL;
  mb.liq_relerr = NULL;
  mb.gas_relerr = NULL;
  mb.boundaries = false;

  /**
  We allocate the fields that store the resuls of the mass
  balance calculations. */

  mb.liq_mass     = (double *)malloc(NLS*sizeof(double));
  mb.gas_mass     = (double *)malloc(NGS*sizeof(double));
  mb.tot_evap     = (double *)malloc(NGS*sizeof(double));
  mb.liq_relerr   = (double *)malloc(NLS*sizeof(double));
  mb.gas_relerr   = (double *)malloc(NGS*sizeof(double));
  mb.liq_massbd   = (double *)malloc(NLS*sizeof(double));
  mb.gas_massbd   = (double *)malloc(NGS*sizeof(double));
  mb.liq_massinit = (double *)malloc(NLS*sizeof(double));
  mb.gas_massinit = (double *)malloc(NGS*sizeof(double));
}

/**
## Initialization

In the initialization event, we write the header of the
output file, in order to facilitate the identification
of the correct column to be plotted. */

/**
### *write_balances()*: write the results in the correct order.
*/

static void write_balances (void) {
  fprintf (fb, "%g %g %g %g ",
      t, mb.totmass1, mb.totmass2, mb.evapmass);
  foreach_elem (mb.YLList, jj)
    fprintf (fb, "%g ", mb.liq_mass[jj]);
  foreach_elem (mb.YGList, jj)
    fprintf (fb, "%g ", mb.gas_mass[jj]);
  foreach_elem (mb.mEvapList, jj)
    fprintf (fb, "%g ", mb.tot_evap[jj]);
  foreach_elem (mb.YLList, jj)
    fprintf (fb, "%g ", mb.liq_relerr[jj]);
  foreach_elem (mb.YGList, jj)
    fprintf (fb, "%g ", mb.gas_relerr[jj]);
  fprintf (fb, "\n");
}

/**
Write the *balances* file header and compute the intial mass of each
phase and species. */

event init0 (i = 0) {
  char name[80];
  sprintf (name, "balances-%d", maxlevel);
  fb = fopen (name, "w");
  fprintf (fb, "time[s](1) totmass1(2) totmass2(3) evapmass(4) ");
  int counter = 5;
  foreach_elem (mb.YLList, jj)
    fprintf (fb, "mL_%s(%d) ", mb.liq_species[jj], counter++);
  foreach_elem (mb.YGList, jj)
    fprintf (fb, "mG_%s(%d) ", mb.gas_species[jj], counter++);
  foreach_elem (mb.mEvapList, jj)
    fprintf (fb, "mE_%s(%d) ", mb.gas_species[jj], counter++);
  foreach_elem (mb.YLList, jj)
    fprintf (fb, "liqerr_%s(%d) ", mb.liq_species[jj], counter++);
  foreach_elem (mb.YGList, jj)
    fprintf (fb, "gaserr_%s(%d) ", mb.gas_species[jj], counter++);
  fprintf (fb, "\n");

  foreach_elem (mb.YLList, jj) {
    mb.liq_mass[jj] = 0.;
    mb.liq_massbd[jj] = 0.;
    mb.liq_massinit[jj] = 0.;
    mb.liq_relerr[jj] = 0.;
  }

  foreach_elem (mb.YGList, jj) {
    mb.gas_mass[jj] = 0.;
    mb.gas_massbd[jj] = 0.;
    mb.gas_massinit[jj] = 0.;
    mb.tot_evap[jj] = 0.;
    mb.gas_relerr[jj] = 0.;
  }

  scalar fb1[], fb2[];
  foreach() {
    fb1[] = f[];
    fb2[] = 1. - f[];
  }
  mb.totmass1init = rho1*statsf(fb1).sum;
  mb.totmass2init = rho2*statsf(fb2).sum;

  foreach(nowarning) {
    foreach_elem (YLList, jj) {
      scalar YL = YLList[jj];
      mb.liq_massinit[jj] += rho1*YL[]*dv();
    }
    foreach_elem (YGList, jj) {
      scalar YG = YGList[jj];
      mb.gas_massinit[jj] += rho2*YG[]*dv();
    }
  }

@if _MPI
  MPI_Allreduce (MPI_IN_PLACE, mb.liq_massinit, NLS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, mb.gas_massinit, NGS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
@endif

  //write_balances();
}

/**
## Cleanup

We de-allocate the dynamic arrays from memory. */

event cleanup (t = end) {
  fclose (fb);
  free (mb.liq_mass);
  free (mb.gas_mass);
  free (mb.tot_evap);
  free (mb.liq_relerr);
  free (mb.gas_relerr);
  free (mb.liq_massbd);
  free (mb.gas_massbd);
  free (mb.liq_massinit);
  free (mb.gas_massinit);
}

/**
## Helper Functions

We define functions that facilitate loops over all the domain
boundaries, adjusting the stencil index according to the boundary ID.

### *face_gradient_bid()*: Return the face gradient at a given boundary.
*/

static double face_gradient_bid (Point point, scalar Y, int bid) {
  double grad = 0.;
  switch (bid) {
    case 0: grad = (Y[1,0]  - Y[])/Delta;  break;  // right
    case 1: grad = (Y[-1,0] - Y[])/Delta;  break;  // left
    case 2: grad = (Y[0,1]  - Y[])/Delta;  break;  // top
    case 3: grad = (Y[0,-1] - Y[])/Delta;  break;  // bottom
  }
  return grad;
}

/**
### *face_flux_bid()*: Return the face flux at a given boundary.
*/

static double face_flux_bid (Point point, scalar Y, int bid) {
  double flux = 0.;
  switch (bid) {
    case 0: flux = 0.5*(Y[1,0]  + Y[])*uf.x[]*Delta; break;  // right
    case 1: flux = 0.5*(Y[-1,0] + Y[])*uf.x[]*Delta; break;  // left
    case 2: flux = 0.5*(Y[0,1]  + Y[])*uf.y[]*Delta; break;  // top
    case 3: flux = 0.5*(Y[0,-1] + Y[])*uf.y[]*Delta; break;  // bottom
  }
  return flux;
}

/**
### *diffusion_boundary()*: Sum the diffusive fluxes across a given boundary.
*/

static void diffusion_boundary (Point point, int bid) {
  foreach_elem (YGList, jj) {
    scalar YG = YGList[jj];
    double gradYG = face_gradient_bid (point, YG, bid);
    mb.gas_massbd[jj] -= rho2*inDmix2[jj]*gradYG*Delta*dt;
  }
}

/**
### *advection_boundary()*: Sum the convective flux across a given boundary.
*/

static void advection_boundary (Point point, int bid) {
  foreach_elem (YGList, jj) {
    scalar YG = YGList[jj];
    double fluxYG = face_flux_bid (point, YG, bid);
    mb.gas_massbd[jj] += rho2*fluxYG*dt;
  }
}

/**
## Compute Balances

First, we write a function that is responsible for the
writing operations, based on the data in the structrue
*mb*.
*/

event balances (i++,last) {

  /**
  We compute the total mass of the gas and of the liquid phase. */

  scalar fb1[], fb2[];
  foreach() {
    fb1[] = f[];
    fb2[] = 1. - f[];
  }
  mb.totmass1 = (mb.totmass1init - rho1*statsf(fb1).sum);
  mb.totmass2 = (rho2*statsf(fb2).sum - mb.totmass2init);

  /**
  We compute the total evaporated mass. */

  scalar mEvapPerVol[];
  foreach() {
    vofrecon vr = vof_reconstruction (point, f);
    mEvapPerVol[] = mEvapTot[]*vr.dirac;
  }
  mb.evapmass += statsf(mEvapPerVol).sum*dt;

  /**
  The total liquid and gas phase masses for each species are set to
  zero before summing up over every cell and boundary. */

  foreach_elem (YGList, jj)
    mb.gas_mass[jj] = 0.;
  foreach_elem (YLList, jj)
    mb.liq_mass[jj] = 0.;

  /**
  We compute the total mass in liquid and gas phase for every
  chemical species of the domain. The total evaporated mass is also
  computed for each species. */

  foreach(nowarning) {
    vofrecon vr = vof_reconstruction (point, f);
    foreach_elem (mb.mEvapList, jj) {
      scalar mEvap = mEvapList[jj];
      scalar YG    = YGList[jj];

      mb.tot_evap[jj] += mEvap[]*vr.dirac*dv()*dt;
      mb.gas_mass[jj] += rho2*YG[]*dv();
    }
    foreach_elem (mb.YLList, jj) {
      scalar YL = mb.YLList[jj];
      mb.liq_mass[jj] += rho1*YL[]*dv();
    }
  }

  /**
  If *boundaries* are included (not by default), the total masses are
  corrected according to the convective and diffusive fluxes through
  the boundaries. */

  if (mb.boundaries) {
    for (int b = 0; b < nboundary; b++) {
      foreach_boundary(b) {
        diffusion_boundary (point, b);
        advection_boundary (point, b);
      }
    }
  }

  /**
  We sum the masses contributions from every processor. */

@if _MPI
  MPI_Allreduce (MPI_IN_PLACE, mb.liq_mass,     NLS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, mb.gas_mass,     NGS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, mb.tot_evap,     NGS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, mb.liq_massbd,   NLS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, mb.gas_massbd,   NGS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, mb.liq_massinit, NLS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, mb.gas_massinit, NGS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
@endif

  foreach_elem (YLList, jj)
    mb.liq_mass[jj] = (mb.liq_massinit[jj] - mb.liq_mass[jj]);
  foreach_elem (YGList, jj)
    mb.gas_mass[jj] = ((mb.gas_mass[jj] + mb.gas_massbd[jj]) - mb.gas_massinit[jj]);

  /**
  We write a new line in the *balances* file, with the updated mass
  values. */

  write_balances();
}

/**
## TODO

1. Add compatibility with metrics
2. Check the bid numbering system, what happens if new bid are
defined? And in 3D?
3. For some reasons it does not work with MPI, and it should be
fixed.
4. For some other reasons *foreach_boundary()* does not work correctly when
a grid/multigrid.h is used.
*/
