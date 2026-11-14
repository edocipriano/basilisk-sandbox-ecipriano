PhaseBalance * liq_balance = NULL, * gas_balance = NULL;
FILE * fbalances = NULL;

#define FORMAT_S "%-15s"
#define FORMAT_G "%-15g"

event init (i = 0) {
  if (nboundary != 4)
    fprintf (ferr,
        "src/balances/two-phase.h:%d: warning: balances considers just 4 boundaries\n",
        LINENO), fflush (ferr);
  // enthalpy: zero for the liquid at TIntVal, dhev for the gas at TIntVal
  liq_balance = new_phase_balance (liq, TIntVal, 0.);
  gas_balance = new_phase_balance (gas, TIntVal, dhev);

  char name[80];
  sprintf (name, "balances-%d", grid->maxdepth);
  fbalances = fopen (name, "w");

  fprintf (fbalances, "time[s](1) totmass1(2) totmass2(3) evapmass(4) ");
  int counter = 5;
  foreach_species_in (liq)
    fprintf (fbalances, "mL_%s(%d) ", liq->species[i], counter++);
  foreach_species_in (gas)
    fprintf (fbalances, "mG_%s(%d) ", gas->species[i], counter++);
  foreach_species_in (gas)
    fprintf (fbalances, "mE_%s(%d) ", gas->species[i], counter++);
  fprintf (fbalances, "totenergy1(%d) totenergy2(%d) qint1(%d) qint2(%d) "
      "evapenergy(%d) ", counter, counter + 1, counter + 2, counter + 3,
      counter + 4);
  fprintf (fbalances, "\n");
}

event cleanup (t = end) {
  delete_phase_balance (liq_balance), liq_balance = NULL;
  delete_phase_balance (gas_balance), gas_balance = NULL;
  fclose (fbalances), fbalances = NULL;
}

static void write_balances (void) {
  fprintf (fbalances, FORMAT_G FORMAT_G FORMAT_G FORMAT_G, t,
      liq_balance->mtot0 - (liq_balance->mtot + liq_balance->mftot),
      (gas_balance->mtot + gas_balance->mftot) - gas_balance->mtot0,
      gas_balance->mevaptot);
  foreach_species_in (liq)
    fprintf (fbalances, FORMAT_G,
        liq_balance->m0[i] - (liq_balance->m[i] + liq_balance->mf[i]));
  foreach_species_in (gas)
    fprintf (fbalances, FORMAT_G,
        (gas_balance->m[i] + gas_balance->mf[i]) - gas_balance->m0[i]);
  foreach_species_in (gas)
    fprintf (fbalances, FORMAT_G, gas_balance->mevap[i]);
  fprintf (fbalances, FORMAT_G FORMAT_G FORMAT_G FORMAT_G FORMAT_G,
      liq_balance->e0 - (liq_balance->e + liq_balance->ef),
      (gas_balance->e + gas_balance->ef) - gas_balance->e0,
      liq_balance->eint, gas_balance->eint,
      gas_balance->mevaptot*gas_balance->h0);
  fprintf (fbalances, "\n");
}

// note: interfacial terms and convective fluxes before the interface moves
// (events with the same name run in reverse order of definition: this one runs
// before the VOF advection), contents and diffusive fluxes at the end of the
// time step
event vof (i++) {
  face vector ufl = (nv > 1) ? uflist[1] : uflist[0];
  face vector ufg = uflist[0];
  phase_balance_transport (liq_balance, liq, NULL, dt, ufl, f);
  phase_balance_transport (gas_balance, gas, mEvapList, dt, ufg, f);
}

event balances (i++) {
  phase_balance_update (liq_balance, liq, dt, f);
  phase_balance_update (gas_balance, gas, dt, f);

  if (pid() == 0)
    write_balances();
}

