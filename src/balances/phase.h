PhaseMassBalance * balance = NULL;
FILE * fbalances = NULL;

#define FORMAT_S "%-15s"
#define FORMAT_G "%-15g"

event init (i = 0) {
  if (nboundary != 4)
    fprintf (ferr,
        "src/balances/two-phase.h:%d: warning: balances considers just 4 boundaries\n",
        LINENO), fflush (ferr);

  balance = new_phase_mass_balance (phase);

  char name[80];
  sprintf (name, "balances-%d", grid->maxdepth);
  fbalances = fopen (name, "w");

  fprintf (fbalances, FORMAT_S, "time(1)");
  fprintf (fbalances, FORMAT_S, "totmass(2)");

  int counter = 3;
  char species[80];
  foreach_species_in (phase) {
    sprintf (species, "m_%s(%d)", phase->species[i], counter++);
    fprintf (fbalances, FORMAT_S, species);
  }
  sprintf (species, "energy(%d)", counter++);
  fprintf (fbalances, FORMAT_S, species);
  fprintf (fbalances, "\n");
}

event cleanup (t = end) {
  delete_phase_mass_balances (balance), balance = NULL;
  fclose (fbalances), fbalances = NULL;
}

static void write_balances (void) {
  fprintf (fbalances, FORMAT_G FORMAT_G, t,
      (balance->mtot + balance->mftot) - balance->mtot0);
  foreach_species_in (phase)
    fprintf (fbalances, FORMAT_G,
        (balance->m[i] + balance->mf[i]) - balance->m0[i]);
  fprintf (fbalances, FORMAT_G, (balance->e + balance->ef) - balance->e0);
  fprintf (fbalances, "\n");
}

// note: convective fluxes before the projection and the diffusion step,
// contents and diffusive fluxes at the end of the time step
event tracer_advection (i++) {
  phase_balance_transport (balance, phase, NULL, dt, uf);
}

event balances (i++) {
  phase_balance_update (balance, phase, dt);

  if (pid() == 0)
    write_balances();
}

