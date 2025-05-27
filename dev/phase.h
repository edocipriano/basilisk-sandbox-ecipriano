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
  char ** species;
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
  scalar * cpList;
  scalar * dhevList;
  double * MWs;
  ThermoState * ts0;
  // Flow divergence
  scalar divu;
  scalar betaT;
  scalar * betaYList;
  scalar DTDt;
  scalar * DYDtList;
} Phase;

macro foreach_scalar_in (Phase * phase) {
  {
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
    scalar divu = phase->divu; NOT_UNUSED (divu);
    scalar betaT = phase->betaT; NOT_UNUSED (betaT);
    scalar DTDt = phase->DTDt; NOT_UNUSED (DTDt);
    {...}
  }
}

macro foreach_species_in (Phase * phase) {
  for (size_t i = 0; i < phase->n; i++) {
    scalar Y = phase->YList[i]; NOT_UNUSED (Y);
    scalar X = phase->XList[i]; NOT_UNUSED (X);
    vector J = phase->JList[i]; NOT_UNUSED (J);
    scalar D = phase->DList[i]; NOT_UNUSED (D);
    scalar cps = phase->cpList[i]; NOT_UNUSED (cps);
    scalar dhevs = phase->dhevList[i]; NOT_UNUSED (dhevs);
    scalar SYexp = phase->SYexpList[i]; NOT_UNUSED (SYexp);
    scalar SYimp = phase->SYimpList[i]; NOT_UNUSED (SYimp);
    scalar betaY = phase->betaYList[i]; NOT_UNUSED (betaY);
    scalar DYDt = phase->DYDtList[i]; NOT_UNUSED (DYDt);
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
  sprintf (name, "%s%s%zu", #Y, phase->name, i+1);                  \
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
    sprintf (name, "%s%s%zu%s", #Y, phase->name, i+1, ext.x);       \
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

void phase_species_names (Phase * phase, char ** names = NULL) {
  if (names) {
    foreach_species_in (phase)
      phase->species[i] = strdup (names[i]);
  }
  else {
    phase->species = malloc (phase->n*sizeof (char *));
    foreach_species_in (phase) {
      char name[80];
      sprintf (name, "%zu", i+1);
      phase->species[i] = strdup (name);
    }
  }
}

Phase * new_phase_empty (char * name = "", bool inverse = false) {
  Phase * phase = (Phase *)malloc (sizeof (Phase));
  phase->name = strdup (name);
  phase->n = 0;
  phase->YList = NULL;
  phase->XList = NULL;
  phase->JList = NULL;
  phase->MWs = NULL;
  phase->ts0 = NULL;
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

  // Default set of names
  phase_species_names (phase);

  // Create minimal set of scalar fields
  new_list_type_name (scalar, Y, phase, phase->YList);
  new_list_type_name (scalar, X, phase, phase->XList);
  new_field_type (scalar, T, phase);

  // Create minimal set of material properties
  new_field_type (scalar, MW, phase);

  foreach() {
    foreach_scalar_in (phase) {
      T[] = 0.;
      MW[] = 0.;
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
  phase->n = ns;
  phase->isothermal = false;
  phase->isomassfrac = false;

  // Default set of names
  phase_species_names (phase);

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
  new_list_type_name (scalar, cp, phase, phase->cpList);
  new_list_type_name (scalar, dhev, phase, phase->dhevList);
  new_field_type (scalar, lambda, phase);
  new_field_type (scalar, cp, phase);
  new_field_type (scalar, dhev, phase);
  new_field_type (scalar, divu, phase);
  new_field_type (scalar, betaT, phase);
  new_field_type (scalar, DTDt, phase);

  // Create source terms
  new_field_type (scalar, STimp, phase);
  new_field_type (scalar, STexp, phase);
  new_list_type_name (scalar, SYimp, phase, phase->SYimpList);
  new_list_type_name (scalar, SYexp, phase, phase->SYexpList);
  new_list_type_name (scalar, betaY, phase, phase->betaYList);
  new_list_type_name (scalar, DYDt, phase, phase->DYDtList);

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
      divu[] = 0.;
      betaT[] = 0.;
      DTDt[] = 0.;
      foreach_species_in (phase) {
        Y[] = 0.;
        X[] = 0.;
        D[] = 0.;
        cps[] = 0.;
        dhevs[] = 0.;
        SYimp[] = 0.;
        SYexp[] = 0.;
        betaY[] = 0.;
        DYDt[] = 0.;
        foreach_dimension()
          J.x[] = 0.;
      }
    }

  // Create species molecular weights and set it to 1
  // fixme: check MWs are correctly provided for each species
  phase->MWs = (double *)malloc (phase->n*sizeof (double));
  foreach_species_in (phase)
    phase->MWs[i] = 1.;

  // Create the initial thermo state
  phase->ts0 = new_thermo_state (phase->n);

  return phase;
}

void delete_phase (Phase * phase) {
  free (phase->name), phase->name = NULL;
  foreach_species_in (phase)
    free (phase->species[i]), phase->species[i] = NULL;
  free (phase->species), phase->species = NULL;
  delete (phase->YList), free (phase->YList), phase->YList = NULL;
  delete (phase->XList), free (phase->XList), phase->XList = NULL;
  delete ((scalar *)phase->JList), free (phase->JList), phase->JList = NULL;
  delete ({phase->P});
  delete ({phase->divu, phase->betaT, phase->DTDt});
  delete (phase->betaYList), free (phase->betaYList), phase->betaYList = NULL;
  delete (phase->DYDtList), free (phase->DYDtList), phase->DYDtList = NULL;
  delete (phase->DList), free (phase->DList), phase->DList = NULL;
  delete (phase->cpList), free (phase->cpList), phase->cpList = NULL;
  delete (phase->dhevList), free (phase->dhevList), phase->dhevList = NULL;
  free (phase->MWs), phase->MWs = NULL;
  free_thermo_state (phase->ts0);
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
      divu[] = 0.;
      DTDt[] = 0.;
      foreach_species_in (phase) {
        DYDt[] = 0.;
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
          : face_value (lambda, 0) /
            (face_value (rho, 0)*face_value (cp, 0) + 1e-16);
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
      foreach_species_in (phase) {
        scalar thetaY[];
        foreach()
          thetaY[] = varcoeff ? cm[]*max (ff[]*rho[], F_ERR)
                              : cm[]*max (ff[], F_ERR);

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
  copy_thermo_state (phase->ts0, ts, phase->n);
  foreach() {
    foreach_scalar_in (phase) {
      T[] = ts->T;
      P[] = ts->P;
      foreach_species_in (phase)
        if (ts->x) Y[] = ts->x[i];
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
        if (provided_list (dhevs)) {
          scalar phase_dhev = phase->dhevList[i];
          phase_dhev[] = dhevs[i];
        }
        if (provided_list (cps)) {
          scalar phase_cp = phase->cpList[i];
          phase_cp[] = cps[i];
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

size_t phase_species_index (Phase * phase, char * species) {
  foreach_species_in (phase) {
    if (strcmp (Y.name, species) == 0)
      return i;
  }
  fprintf (ferr, "src/phase.h:%d: error: species %s not found\n",
      LINENO, species), fflush (ferr);
  exit(1);
}

// Usage: phase_set_composition_from_string (phase, "NC7H16 0.2 N2 0.8");
void phase_set_composition_from_string (Phase * phase, char * s) {
  foreach() {
    foreach_species_in (phase) {
      phase->ts0->x[i] = 0.;
      Y[] = 0.;
    }
  }

  char * input = strdup (s);
  char * token = strtok (input, " ");
  unsigned int count = 0;
  while (token != NULL) {
    count++;
    char * species;
    double val = 0.;
    if (count%2 != 0) {
      species = token;
      val = 0;
    }
    else
      val = atof (token);

    if (val) {
      size_t index = phase_species_index (phase, species);
      scalar Y = phase->YList[index];
      phase->ts0->x[index] = val;
      foreach()
        Y[] = val;
    }
    token = strtok (NULL, " ");
  }
  free (input);
}

bool phase_is_uniform (Phase * phase) {
  bool uniform = false;
  foreach_scalar_in (phase) {
    uniform |= (statsf (T).stddev == 0);
    foreach_species_in (phase)
      uniform |= (statsf (Y).stddev == 0);
  }
  return uniform;
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
    double tol = 1e-10, bool extend = false)
{
  double * ymass = (double *)malloc (phase->n*sizeof (double));
  double * xmole = (double *)malloc (phase->n*sizeof (double));
  foreach_scalar_in (phase) {
    foreach() {
      //double ff = phase->inverse ? 1. - f[] : f[];
      //if (ff > tol) {
        foreach_species_in (phase)
          ymass[i] = (phase->n == 1.) ? 1. : Y[];
        correctfrac (ymass, phase->n);
        mass2molefrac (xmole, ymass, phase->MWs, phase->n);
        MW[] = mass2mw (ymass, phase->MWs, phase->n);
        foreach_species_in (phase)
          X[] = xmole[i];
      //}
    }

    if (extend) {
      foreach() {
        double ff = phase->inverse ? 1. - f[] : f[];
        if (ff <= tol) {
          int counter = 0;
          double ext_MW = 0.;
          foreach_neighbor(1) {
            double ffnei = phase->inverse ? 1. - f[] : f[];
            if (ffnei > tol) {
              counter++;
              ext_MW += MW[];
            }
          }
          MW[] = (counter != 0.) ? ext_MW/counter : 0.;
        }
      }
    }
  }
  free (ymass);
  free (xmole);
}

void phase_normalize_fractions (Phase * phase) {
  double * ymass = (double *)malloc (phase->n*sizeof (double));
  foreach() {
    foreach_species_in (phase)
      ymass[i] = Y[];
    correctfrac (ymass, phase->n);
    foreach_species_in (phase)
      Y[] = ymass[i];
  }
  free (ymass);
}

void phase_update_properties (Phase * phase, const ThermoProps * tp,
    (const) scalar f = unity, double tol = 1e-10)
{
  double * xmole = (double *)malloc (phase->n*sizeof (double));
  double * arrdiff = (double *)malloc (phase->n*sizeof (double));
  double * arrbetaY = (double *)malloc (phase->n*sizeof (double));
  double * arrdhev = (double *)malloc (phase->n*sizeof (double));
  double * arrcps  = (double *)malloc (phase->n*sizeof (double));
  foreach() {
    double ff = phase->inverse ? 1. - f[] : f[];
    if (ff > tol) {
      ThermoState ts;
      foreach_scalar_in (phase) {
        foreach_species_in (phase)
          xmole[i] = X[];
        ts.T = T[];
        ts.P = P[];
        ts.x = (phase->n == 1) ? (double[]){1.} : xmole;

        if (tp->rhov) rho[] = tp->rhov (&ts);
        if (tp->muv) mu[] = tp->muv (&ts);
        if (tp->lambdav) lambda[] = tp->lambdav (&ts);
        if (tp->cpv) cp[] = tp->cpv (&ts);
        if (tp->betaT) betaT[] = tp->betaT (tp, &ts);

        if (tp->diff) tp->diff (&ts, arrdiff);
        if (tp->betaY) tp->betaY (tp, &ts, arrbetaY);
        if (tp->dhev) tp->dhev (&ts, arrdhev);
        if (tp->cps)  tp->cps  (&ts, arrcps);

        foreach_species_in (phase) {
          if (tp->diff) D[] = arrdiff[i];
          if (tp->betaY) betaY[] = arrbetaY[i];
          if (tp->dhev) dhevs[] = arrdhev[i];
          if (tp->cps) cps[] = arrcps[i];
        }
      }
    }
  }
  free (xmole);
  free (arrdiff);
  free (arrbetaY);
  free (arrdhev);
  free (arrcps);
}

#define increment_property(ext_s, s) \
  ext_s += (s.i > 0) ? s[] : 0;

#define assign_property(ext_s, s, counter) \
  if (s.i > 0) s[] = (counter > 0) ? ext_s/counter : 0.;

void phase_extend_properties (Phase * phase,
    (const) scalar f = unity, double tol = 1e-10)
{
  double * ext_diff = (double *)malloc (phase->n*sizeof (double));
  double * ext_cps = (double *)malloc (phase->n*sizeof (double));
  double * ext_dhevs = (double *)malloc (phase->n*sizeof (double));
  double * ext_betaY = (double *)malloc (phase->n*sizeof (double));

  foreach_scalar_in (phase) {
    foreach() {
      double ff = phase->inverse ? 1. - f[] : f[];
      if (ff <= tol) {
        double ext_rho = 0.;
        double ext_mu = 0.;
        double ext_MW = 0.;
        double ext_lambda = 0.;
        double ext_cp = 0.;
        double ext_dhev = 0.;
        double ext_betaT = 0.;

        foreach_species_in (phase) {
          ext_diff[i] = 0.;
          ext_cps[i] = 0.;
          ext_dhevs[i] = 0.;
          ext_betaY[i] = 0.;
        }

        int counter = 0;
        foreach_neighbor(1) {
          double ffnei = phase->inverse ? 1. - f[] : f[];
          if (ffnei > tol) {
            counter++;
            increment_property (ext_rho, rho);
            increment_property (ext_mu, mu);
            increment_property (ext_MW, MW);
            increment_property (ext_lambda, lambda);
            increment_property (ext_cp, cp);
            increment_property (ext_dhev, dhev);
            increment_property (ext_betaT, betaT);

            foreach_species_in (phase) {
              increment_property (ext_diff[i], D);
              increment_property (ext_cps[i], cps);
              increment_property (ext_dhevs[i], dhevs);
              increment_property (ext_betaY[i], betaY);
            }
          }
        }
        assign_property (ext_rho, rho, counter);
        assign_property (ext_mu, mu, counter);
        assign_property (ext_MW, MW, counter);
        assign_property (ext_lambda, lambda, counter);
        assign_property (ext_cp, cp, counter);
        assign_property (ext_dhev, dhev, counter);
        assign_property (ext_betaT, betaT, counter);

        foreach_species_in (phase) {
          assign_property (ext_diff[i], D, counter);
          assign_property (ext_cps[i], cps, counter);
          assign_property (ext_dhevs[i], dhevs, counter);
          assign_property (ext_betaY[i], betaY, counter);
        }
      }
    }
  }
  free (ext_diff);
  free (ext_cps);
  free (ext_dhevs);
  free (ext_betaY);
}

void phase_update_divergence (Phase * phase,
    (const) scalar f = unity,
    bool fick_corrected = true, bool molar_diffusion = true)
{
  scalar ff[];
  foreach()
    ff[] = phase->inverse ? 1. - f[] : f[];

  face vector fs[];
  face_fraction (ff, fs); // fixme: can't use f in this function

  foreach_scalar_in (phase) {

    /**
    We calculate the Lagrangian derivative of the temperature field. */

    face vector lambdagrad[];
    foreach_face()
      lambdagrad.x[] = face_value (lambda, 0)*face_gradient_x (T, 0) *
        fm.x[]*fs.x[];

    foreach() {
      foreach_dimension()
        DTDt[] += (lambdagrad.x[1] - lambdagrad.x[])/Delta;
      DTDt[] += STexp[];
    }

    /**
    We calculate the Lagrangian derivative of the chemical species mass fractions.
    */

    foreach_species_in (phase) {
      face vector rhoDmixY[];
      foreach_face() {
        double rhof = face_value (rho, 0);
        double Df = face_value (D, 0);
        rhoDmixY.x[] = rhof*Df*face_gradient_x (Y, 0)*fm.x[]*fs.x[];
      }

      foreach() {
        foreach_dimension()
          DYDt[] += (rhoDmixY.x[1] - rhoDmixY.x[])/Delta;
        DYDt[] += (SYexp[] + SYimp[]*Y[]); // fixme: I was using the YGInt
      }
    }

    /**
    We add diffusion correction contribution to the chemical species mass
    fraction derivatives. */

    face vector phictot[];
    foreach_face() {
      phictot.x[] = 0.;
      if (fick_corrected) {
        foreach_species_in (phase) {
          double rhof = face_value (rho, 0);
          double Df = face_value (D, 0);
          if (molar_diffusion) {
            double MWf = face_value (MW, 0);
            phictot.x[] += (MWf > 0.) ?
              rhof*Df*phase->MWs[i]/MWf*face_gradient_x (X, 0)*fm.x[]*fs.x[] :
              0.;
          }
          else {
            phictot.x[] += rhof*Df*face_gradient_x (Y, 0)*fm.x[]*fs.x[];
          }
        }
      }
    }

    foreach_species_in (phase) {
      face vector phic[];
      foreach_face() {
        phic.x[] = phictot.x[];
        if (molar_diffusion) {
          double rhof = face_value (rho, 0);
          double Df = face_value (D, 0);
          double MWf = face_value (MW, 0);

          phic.x[] -= (MWf > 0.) ?
            rhof*Df/MWf*face_gradient_x (MW, 0)*fm.x[]*fs.x[] : 0.;
        }
        phic.x[] *= face_value (Y, 0);
      }

      foreach()
        foreach_dimension()
          DYDt[] -= (phic.x[1] - phic.x[])/Delta;
    }

    /**
    We calculate the divergence of the phase, just in the region occupied by the
    phase, defined by the volume fraction `f`. */

    foreach() {
      divu[] = 0.;

      // Add temperature contribution
      divu[] += (rho[]*cp[] > 0.) ? betaT[]/(rho[]*cp[])*DTDt[] : 0.;

      // Add chemical species contribution
      double divuspecies = 0.;
      foreach_species_in (phase)
        divuspecies += (rho[] > 0.) ? betaY[]/rho[]*DYDt[] : 0.;
      divu[] += divuspecies;

      // Volume-averaged divergence
      divu[] *= ff[];

      // Adjust sign for internal convention
      divu[] *= -1;
    }

  }
}

void phase_add_heat_species_diffusion (Phase * phase, const scalar f = unity,
    bool molar_diffusion = false, double tol = 1e-10)
{
  if (!phase->isothermal) {
    foreach_scalar_in (phase) {
      foreach() {
        double mde = 0.;
        coord gT = {0.,0.,0.};
        coord gY = {0.,0.,0.};
        coord gYsum = {0.,0.,0.};

        foreach_dimension()
          gT.x = (T[1] - T[-1])/(2.*Delta);

        foreach_dimension() {
          foreach_species_in (phase) {
            if (molar_diffusion)
              gYsum.x -= (MW[] > 0.) ?
                rho[]*D[]*phase->MWs[i]/MW[]*(X[1] - X[-1])/(2.*Delta) : 0.;
            else
              gYsum.x -= rho[]*D[]*(Y[1] - Y[-1])/(2.*Delta);
          }

          foreach_species_in (phase) {
            if (molar_diffusion)
              gY.x = (MW[] > 0.) ?
                -rho[]*D[]*phase->MWs[i]/MW[]*(X[1] - X[-1])/(2.*Delta) : 0.;
            else
              gY.x = -rho[]*D[]*(Y[1] - Y[-1])/(2.*Delta);
            mde += cps[]*(gY.x - Y[]*gYsum.x)*gT.x;
          }
        }
        double ff = phase->inverse ? 1. - f[] : f[];
        STexp[] -= mde*cm[]*(ff == 1);  // fixme: use interfacial gradients if
                                        // you want to apply it to interfacial
                                        // cells
      }
    }
  }
}

void phase_diffusion_velocity (Phase * phase, const scalar f = unity,
    bool fick_corrected = true, bool molar_diffusion = true)
{
  if (!phase->isomassfrac) {
    scalar ff[];
    foreach()
      ff[] = phase->inverse ? 1. - f[] : f[];

    face vector fs[];
    face_fraction (ff, fs); // fixme: can't use f in this function

    foreach_scalar_in (phase) {
      face vector phictot[];
      foreach_face() {
        phictot.x[] = 0.;
        if (fick_corrected) {
          foreach_species_in (phase) {
            double rhof = face_value (rho, 0);
            double Df = face_value (D, 0);
            if (molar_diffusion) {
              double MWf = face_value (MW, 0);
              phictot.x[] += (MWf > 0.) ?
                rhof*Df*phase->MWs[i]/MWf*face_gradient_x (X, 0) : 0.;
            }
            else
              phictot.x[] += rhof*Df*face_gradient_x (Y, 0);
          }
        }
      }

      foreach_species_in (phase) {
        face vector phic[];
        foreach_face() {
          phic.x[] = phictot.x[];
          if (molar_diffusion) {
            double rhof = face_value (rho, 0);
            double Df = face_value (D, 0);
            double MWf = face_value (MW, 0);
            phic.x[] -= (MWf > 0.) ? rhof*Df/MWf*face_gradient_x (MW, 0) : 0.;
          }
          phic.x[] *= fm.x[]*fs.x[];
        }
        double (* gradient_backup)(double, double, double) = Y.gradient;
        Y.gradient = NULL;
        face vector flux[];
        tracer_fluxes (Y, phic, flux, dt, zeroc);
        Y.gradient = gradient_backup;

        foreach()
          foreach_dimension()
            Y[] += (rho[] > 0.) ? dt/rho[]*(flux.x[] - flux.x[1])/(Delta*cm[]) : 0.;
      }
    }
  }
}

