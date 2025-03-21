#include "variable-properties.h"
#include "fracface.h"
#include "diffusion.h"

#define UNDEFINED -1
#ifndef F_ERR
# define F_ERR 1e-10
#endif

typedef struct {
  // Public attributes
  char * name;
  scalar T;
  scalar P;
  scalar STimp;
  scalar STexp;
  scalar * YList;
  scalar * XList;
  vector * JList;
  scalar * SYimpList;
  scalar * SYexpList;
  size_t n;
  bool isothermal;
  bool isomassfrac;
  bool inverse;
  scalar * tracers;
  // Material properties
  scalar rho;
  scalar mu;
  scalar MW;
  scalar lambda;
  scalar cp;
  scalar dhev;
  scalar * DList;
  double * MWs;
} Phase;

macro foreach_scalar_in (Phase * phase) {
  scalar T = phase->T; NOT_UNUSED (T);
  scalar P = phase->P; NOT_UNUSED (P);
  scalar rho = phase->rho; NOT_UNUSED (rho);
  scalar mu = phase->mu; NOT_UNUSED (mu);
  scalar MW = phase->MW; NOT_UNUSED (MW);
  scalar lambda = phase->lambda; NOT_UNUSED (lambda);
  scalar cp = phase->cp; NOT_UNUSED (cp);
  scalar dhev = phase->dhev; NOT_UNUSED (dhev);
  scalar STexp = phase->STexp; NOT_UNUSED (STexp);
  scalar STimp = phase->STimp; NOT_UNUSED (STimp);
  {...}
}

macro foreach_species_in (Phase * phase) {
  for (size_t i = 0; i < phase->n; i++) {
    scalar Y = phase->YList[i]; NOT_UNUSED (Y);
    scalar X = phase->XList[i]; NOT_UNUSED (X);
    vector J = phase->JList[i]; NOT_UNUSED (J);
    scalar D = phase->DList[i]; NOT_UNUSED (D);
    scalar SYexp = phase->SYexpList[i]; NOT_UNUSED (SYexp);
    scalar SYimp = phase->SYimpList[i]; NOT_UNUSED (SYimp);
    {...}
  }
}

#define new_field_scalar(Y, phase)                                  \
  {                                                                 \
    scalar Y = new scalar;                                          \
    phase->Y = Y;                                                   \
    Y.inverse = phase->inverse;                                     \
    char name[80];                                                  \
    sprintf (name, "%s%s", #Y, phase->name);                        \
    free (Y.name);                                                  \
    Y.name = strdup (name);                                         \
  }

#define new_field_type(type, Y, phase)                              \
  new_field_##type(Y, phase);

#define new_field_scalar_name(Y, phase)                             \
  scalar Y = new scalar;                                            \
  Y.inverse = phase->inverse;                                       \
  char name[80];                                                    \
  sprintf (name, "%s%s%zu", #Y, phase->name, i);                    \
  free (Y.name);                                                    \
  Y.name = strdup (name);

#define new_list_scalar_name(Y, phase, list)                        \
  for (size_t i = 0; i < phase->n; i++) {                           \
    new_field_scalar_name (Y, phase);                               \
    list = list_add (list, Y);                                      \
  }

#define new_field_vector_name(Y, phase)                             \
  vector Y = new vector;                                            \
  char name[80];                                                    \
  foreach_dimension() {                                             \
    Y.x.inverse = phase->inverse;                                   \
    sprintf (name, "%s%s%zu%s", #Y, phase->name, i, ext.x);         \
    free (Y.x.name);                                                \
    Y.x.name = strdup (name);                                       \
  }                                                                 \

#define new_list_vector_name(Y, phase, list)                        \
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};          \
  for (size_t i = 0; i < phase->n; i++) {                           \
    new_field_vector_name(Y, phase);                                \
    list = vectors_add (list, Y);                                   \
  }

#define new_list_type_name(type, Y, phase, list)                    \
  new_list_##type##_name(Y, phase, list)

Phase * new_phase_empty (char * name = "", bool inverse = false) {
  Phase * phase = (Phase *)malloc (sizeof (Phase));
  phase->name = strdup (name);
  phase->n = 0;
  phase->YList = NULL;
  phase->XList = NULL;
  phase->JList = NULL;
  phase->MWs = NULL;
  phase->isothermal = true;
  phase->isomassfrac = true;
  phase->inverse = inverse;
  phase->tracers = NULL;

  return phase;
}

Phase * new_phase_minimal (char * name = "", size_t ns = 0,
    bool inverse = false)
{
  Phase * phase = new_phase_empty (name, inverse);
  phase->n = ns;

  // Create minimal set of scalar fields
  new_list_type_name (scalar, Y, phase, phase->YList);
  new_list_type_name (scalar, X, phase, phase->XList);
  new_field_type (scalar, T, phase);

  foreach() {
    foreach_scalar_in (phase) {
      T[] = 0.;
      foreach_species_in (phase) {
        Y[] = 0.;
        X[] = 0.;
      }
    }
  }

  return phase;
}

Phase * new_phase (char * name = "", size_t ns = 0, bool inverse = false) {
  Phase * phase = new_phase_empty (name, inverse);
  phase->isothermal = false;
  phase->isomassfrac = false;

  // Create scalar fields
  new_field_type (scalar, T, phase);
  new_list_type_name (scalar, Y, phase, phase->YList);
  new_list_type_name (scalar, X, phase, phase->XList);
  new_list_type_name (vector, J, phase, phase->JList);

  // Create material properties
  new_field_type (scalar, P, phase);
  new_field_type (scalar, rho, phase);
  new_field_type (scalar, mu, phase);
  new_field_type (scalar, MW, phase);
  new_list_type_name (scalar, D, phase, phase->DList);
  new_field_type (scalar, lambda, phase);
  new_field_type (scalar, cp, phase);
  new_field_type (scalar, dhev, phase);

  // Create source terms
  new_field_type (scalar, STimp, phase);
  new_field_type (scalar, STexp, phase);
  new_list_type_name (scalar, SYimp, phase, phase->SYimpList);
  new_list_type_name (scalar, SYexp, phase, phase->SYexpList);

  foreach()
    foreach_scalar_in (phase) {
      P[] = 0.;
      rho[] = 0.;
      mu[] = 0.;
      MW[] = 1.;
      lambda[] = 0.;
      cp[] = 0.;
      dhev[] = 0.;
      STimp[] = 0.;
      STexp[] = 0.;
      foreach_species_in (phase) {
        Y[] = 0.;
        X[] = 0.;
        D[] = 0.;
        SYimp[] = 0.;
        SYexp[] = 0.;
        foreach_dimension()
          J.x[] = 0.;
      }
    }

  // Create species molecular weights and set it to 1
  phase->MWs = (double *)malloc (phase->n*sizeof (double));
  foreach_species_in (phase)
    phase->MWs[i] = 1.;

  return phase;
}

void delete_phase (Phase * phase) {
  free (phase->name), phase->name = NULL;
  delete (phase->YList), free (phase->YList), phase->YList = NULL;
  delete (phase->XList), free (phase->XList), phase->XList = NULL;
  delete ((scalar *)phase->JList), free (phase->JList), phase->JList = NULL;
  delete ({phase->P});
  free (phase->MWs), phase->MWs = NULL;
  // Free the temperature using the attribute _freeme
  free (phase->tracers), phase->tracers = NULL;
  free (phase), phase = NULL;
}

void phase_set_tracers (Phase * phase) {
  if (!phase->isomassfrac)
    phase->tracers = list_concat (phase->tracers, phase->YList);
  if (!phase->isothermal)
    phase->tracers = list_add (phase->tracers, phase->T);
}

void phase_reset_sources (Phase * phase) {
  foreach() {
    foreach_scalar_in (phase) {
      STexp[] = 0.;
      STimp[] = 0.;
      foreach_species_in (phase) {
        SYexp[] = 0.;
        SYimp[] = 0.;
      }
    }
  }
}

void phase_diffusion (Phase * phase, (const) scalar f = unity,
    bool varcoeff = false)
{
  scalar ff[];
  foreach()
    ff[] = phase->inverse ? 1. - f[] : f[];

  face vector fs[];
  face_fraction (ff, fs); // fixme: can't use f in this function

  foreach_scalar_in (phase) {
    if (!phase->isothermal) {
      scalar thetaT[];
      foreach()
        thetaT[] = varcoeff ? cm[]*max (ff[]*rho[]*cp[], F_ERR)
                            : cm[]*max (ff[], F_ERR);

      face vector lambdaf[];
      foreach_face() {
        lambdaf.x[] = varcoeff ? face_value (lambda, 0)
          : face_value (lambda, 0)/face_value (rho, 0)/face_value (cp, 0);
        lambdaf.x[] *= fs.x[]*fm.x[];
      }

      if (!varcoeff)
        foreach() {
          STexp[] = (rho[]*cp[] > 0.) ? STexp[]/(rho[]*cp[]) : STexp[];
          STimp[] = (rho[]*cp[] > 0.) ? STimp[]/(rho[]*cp[]) : STimp[];
        }

      diffusion (T, dt, D=lambdaf, r=STexp, beta=STimp, theta=thetaT);
    }

    if (!phase->isomassfrac) {
      scalar thetaY[];
      foreach()
        thetaY[] = varcoeff ? cm[]*max (ff[]*rho[], F_ERR)
                            : cm[]*max (ff[], F_ERR);

      foreach_species_in (phase) {
        face vector Df[];
        foreach_face() {
          Df.x[] = varcoeff ? face_value (rho, 0)*face_value (D, 0)
                            : face_value (D, 0);
          Df.x[] *= fs.x[]*fm.x[];
        }

        if (!varcoeff)
          foreach() {
            SYexp[] = (rho[] > 0.) ? SYexp[]/rho[] : SYexp[];
            SYimp[] = (rho[] > 0.) ? SYimp[]/rho[] : SYimp[];
          }

        diffusion (Y, dt, D=Df, r=SYexp, beta=SYimp, theta=thetaY);
      }

    }
  }

  phase_reset_sources (phase);  // fixme: maybe I still need those terms
}

void phase_set_thermo_state (Phase * phase, const ThermoState * ts) {
  foreach() {
    foreach_scalar_in (phase) {
      T[] = ts->T;
      P[] = ts->P;
      foreach_species_in (phase)
        Y[] = ts->x[i];
    }
  }
}

#define provided(x) (x != UNDEFINED)

#define provided_list(x) (x != NULL)

#define error_unprovided(message)                                             \
  {                                                                           \
    fprintf (ferr, "src/phase.h:%d: error: %s not provided\n",                \
        LINENO, #message);                                                    \
    fflush(ferr);                                                             \
    abort();                                                                  \
  }

#define check_provided(x)                                                     \
  if (!provided (x))                                                          \
    error_unprovided (x)

#define check_provided_list(x)                                                \
  if (!provided_list (x))                                                     \
    error_unprovided (x)

void phase_set_properties (Phase * phase,
    double rho = UNDEFINED, double mu = UNDEFINED,
    double lambda = UNDEFINED, double cp = UNDEFINED,
    double dhev = UNDEFINED, double * dhevs = NULL,
    double * D = NULL, double * MWs = NULL,
    double * cps = NULL)
{
  if (provided_list (MWs))
    foreach_species_in (phase)
      phase->MWs[i] = MWs[i];

  foreach() {
    if (provided (rho)) {
      scalar phase_rho = phase->rho;
      phase_rho[] = rho;
    }
    if (provided (mu)) {
      scalar phase_mu = phase->mu;
      phase_mu[] = mu;
    }
    if (provided (dhev)) {
      scalar phase_dhev = phase->dhev;
      phase_dhev[] = dhev;
    }

    if (!phase->isomassfrac) {
      for (size_t i = 0; i < phase->n; i++) {
        if (provided_list (D)) {
          scalar phase_D = phase->DList[i];
          phase_D[] = D[i];
        }
      }
    }

    if (!phase->isothermal) {
      if (provided (lambda)) {
        scalar phase_lambda = phase->lambda;
        phase_lambda[] = lambda;
      }
      if (provided (cp)) {
        scalar phase_cp = phase->cp;
        phase_cp[] = cp;
      }
    }

  }
}

void phase_tracers_to_scalars (Phase * phase, scalar f, double tol = 1e-10) {
  foreach() {
    double ff = phase->inverse ? 1. - f[] : f[];
    foreach_scalar_in (phase) {
      if (!phase->isothermal)
        T[] = (ff > tol) ? T[]/ff : 0.;
      foreach_species_in (phase)
        if (!phase->isomassfrac)
          Y[] = (ff > tol) ? Y[]/ff : 0.;
    }
  }
}

void phase_scalars_to_tracers (Phase * phase, scalar f) {
  foreach() {
    double ff = phase->inverse ? 1. - f[] : f[];
    foreach_scalar_in (phase) {
      if (!phase->isothermal)
        T[] *= phase->isothermal ? 1. : ff;
      foreach_species_in (phase)
        if (!phase->isomassfrac)
          Y[] *= ff;
    }
  }
}

void phase_update_mw_moles (Phase * phase, (const) scalar f = unity,
    double tol = 1e-10, bool extend = true)
{
  double * ymass = (double *)malloc (phase->n*sizeof (double));
  double * xmolar = (double *)malloc (phase->n*sizeof (double));
  foreach() {
    double ff = phase->inverse ? 1. - f[] : f[];
    if (ff > tol) {
      foreach_scalar_in (phase) {
        foreach_species_in (phase)
          ymass[i] = Y[];
        correctfrac (ymass, phase->n);
        mass2molefrac (xmolar, ymass, phase->MWs, phase->n);
        MW[] = mass2mw (ymass, phase->MWs, phase->n);
        foreach_species_in (phase)
          X[] = xmolar[i];
      }
    }
  }
  free (ymass);
  free (xmolar);
}

void phase_update_properties (Phase * phase, const ThermoProps * tp,
    (const) scalar f = unity, double tol = 1e-10)
{
  double * xmolar = (double *)malloc (phase->n*sizeof (double));
  foreach() {
    double ff = phase->inverse ? 1. - f[] : f[];
    if (ff > tol) {
      ThermoState ts;
      foreach_scalar_in (phase) {
        foreach_species_in (phase)
          xmolar[i] = X[];
        ts.T = T[];
        ts.P = P[];
        ts.x = (phase->n == 1) ? (double[]){1.} : xmolar;

        if (tp->rhov) rho[] = tp->rhov (&ts);
        if (tp->muv) mu[] = tp->muv (&ts);
        if (tp->lambdav) lambda[] = tp->lambdav (&ts);
        if (tp->cpv) cp[] = tp->cpv (&ts);

        foreach_species_in (phase)
          if (tp->diff) D[] = tp->diff (&ts, i);
      }
    }
  }
  free (xmolar);
}

void phase_extend_properties (Phase * phase,
    (const) scalar f = unity, double tol = 1e-10)
{}

