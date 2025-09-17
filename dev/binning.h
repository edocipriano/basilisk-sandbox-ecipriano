/**
# Binning/Cell-Agglomeration Technique

We wish to accelerate the direct integration of the reactive step by grouping
together cells with similar thermodynamic states in *bins*.  Finally, the bins
are integrated as if they were a single computational cell, and the results are
then re-mapped into the cell values in a conservative manner.  We follow the
algorithm proposed by [Goldin et al. 2009](#goldin2009cell).

## Usage

This algorithm can be written in a few lines of code using the functions and the
smart iterators defined in this module. The logic is:

* *binning*: divide the domain in a number of bins are return the table
* *integration*: integrate the ODE system each bin
* *remap*: map the bin values back to the cells

This procedure corresponds to the following code:

~~~bash
BinTable * table = binning (target, eps);

foreach_bin (table) {
  foreach_bin_target (bin) {
    target0 = target;
    target += RHS*dt;
  }
}

binning_remap (table, targets);
binning_cleanup (table);
~~~
*/

/**
## Bin

A bin is an aggregate of cells. In practice, it is a structure containing a
global ID, the number of target fields, the number of cells within the bin, and
the arrays `phi` and `phi0` which are the mass-averaged target scalar value
within the bin (with dimension equal to `ntargets`).  The cells are stored using
the coordinate of the centroids, stored in the `cells` list of coordinates. */

typedef struct {
  size_t id;
  size_t ntargets;
  size_t ncells;
  double * phi, * phi0;
  coord * cells;
} Bin;

Bin * new_bin (size_t id, size_t ntargets) {
  Bin * bin = malloc (sizeof (Bin));
  bin->id = id;
  bin->ntargets = ntargets;
  bin->ncells = 0;
  bin->phi = malloc (ntargets*sizeof (double));
  bin->phi0 = malloc (ntargets*sizeof (double));
  bin->cells = NULL;
  return bin;
}

void free_bin (Bin * bin) {
  free (bin->phi);
  free (bin->phi0);
  free (bin->cells);
  free (bin), bin = NULL;
}

void bin_append_cell (Point point, Bin * bin) {
  bin->ncells++;
  bin->cells = realloc (bin->cells, bin->ncells*sizeof (coord));

  coord o = {x,y,z};
  foreach_dimension()
    bin->cells[bin->ncells-1].x = o.x;
}

/**
## Bin Table

A bin table collects all the bins in the domain, and it facilitates loops and
operations over the bins. */

typedef struct {
  size_t ntargets;
  size_t nbins;
  Bin ** bins;
} BinTable;

BinTable * new_bintable (size_t nbins, size_t ntargets) {
  BinTable * table = malloc (sizeof (BinTable));
  table->ntargets = ntargets;
  table->nbins = nbins;
  table->bins = malloc (nbins*sizeof (Bin *));
  for (size_t i = 0; i < nbins; i++)
    table->bins[i] = new_bin (i, ntargets);
  return table;
}

void free_bintable (BinTable * table) {
  for (size_t i = 0; i < table->nbins; i++)
    free (table->bins[i]);
  free (table), table = NULL;
}

void bintable_update (BinTable * table, size_t nbins) {
  table->nbins = nbins;
  table->bins = malloc (nbins*sizeof (Bin *));
  for (size_t i = 0; i < nbins; i++)
    table->bins[i] = new_bin (i, table->ntargets);
}

/**
## Iterators

We define useful iterators to simplify loops over the bins and their cells.  In
particular `foreach_bin` loops over each bin in the bin table;
`foreach_bin_target` iterates over the targets of a bin; `foreach_bin_cell`
iterates over the cells of a bin. */

macro foreach_bin (BinTable * table) {
  for (size_t i = 0; i < table->nbins; i++) {
    Bin * bin = table->bins[i]; NOT_UNUSED (bin);
    if (bin->ncells > 0) // active bin
      {...}
  }
}

macro foreach_bin_target (Bin * bin) {
  for (size_t j = 0; j < bin->ntargets; j++) {
    double target = bin->phi[j], target0 = bin->phi0[j];
    {...}
    bin->phi[j] = target, bin->phi0[j] = target0;
  }
}

macro foreach_bin_cell (Bin * bin) {
  for (size_t _k = 0; _k < bin->ncells; _k++) {
    coord o = bin->cells[_k];
    foreach_point (o.x,o.y,o.z) {
      {...}
    }
  }
}

/**
## Program interface

The following functions can be called by the user to perform the binning,
remap the bins, and clean the memory. */

void binning_normalize (scalar * targets) {
  int len = list_len (targets);
  double * tmax = malloc (len*sizeof (double));
  double * tmin = malloc (len*sizeof (double));

  for (int i = 0; i < len; i++) {
    scalar t = targets[i];
    tmax[i] = statsf (t).max;
    tmin[i] = statsf (t).min;
  }

  foreach() {
    for (int i = 0; i < len; i++) {
      scalar t = targets[i];
      t[] = (t[] - tmin[i])/(tmax[i] - tmin[i]);
    }
  }

  free (tmax);
  free (tmin);
}

// fixme: avoid empty bins
BinTable * binning_build_table (scalar * targets, double eps) {
  double L = 1./eps;
  size_t len = (size_t)list_len (targets);

  // All the possible combinations
  size_t nbins = 1;
  for (size_t i = 0; i < len; i++)
    nbins *= (L + 1);

  BinTable * table = new_bintable (nbins, len);

  foreach() {
    size_t bin_j = 0;
    for (size_t i = 0; i < len; i++) {
      scalar t = targets[i];
      size_t bin_i = (size_t)(t[]*L);
      if (bin_i > L) bin_i = L;
      bin_j += bin_i * pow (L + 1, i);
    }

    Bin * bin = table->bins[bin_j];
    bin_append_cell (point, bin);
  }
  return table;
}

// fixme: mass averaged
void binning_average (BinTable * table, scalar * targets) {
  foreach_bin (table) {
    foreach_bin_target (bin) {
      scalar t = targets[j];
      double num = 0., den = 0.;
      foreach_bin_cell (bin) {
        num += t[]*dv();
        den += dv();
      }
      //assert (den > 0.);  // this is an empty bin
      target = (den > 0.) ? num / den : 0;
      target0 = target;
    }
  }
}

void binning_remap (BinTable * table, scalar * targets) {
  foreach_bin (table) {
    foreach_bin_cell (bin) {
      foreach_bin_target (bin) {
        scalar t = targets[j];
        t[] += (target - target0);
      }
    }
  }
}

void binning_cleanup (BinTable * table) {
  free_bintable (table);
}

void binning_ids (const BinTable * table, scalar ids) {
  foreach_bin (table) {
    foreach_bin_cell (bin) {
      ids[] = bin->id;
    }
  }
}

typedef struct {
  // total number of (active/non-empty) bins
  size_t nactive;
  // average number of cells per bin
  size_t navg;
  // maximum number of cells in a bin
  size_t nmax;
  // minimum number of cells in a bin
  size_t nmin;
} bstats;

bstats binning_stats (const BinTable * table) {
  bstats s = {0};
  s.nmax = 0, s.nmin = 1e9;
  foreach_bin (table) {
    s.nactive++;
    s.nmax = max (s.nmax, bin->ncells);
    s.nmin = min (s.nmin, bin->ncells);
  }
  s.navg = (s.nactive > 0) ? grid->tn / s.nactive : 0;
  return s;
}

BinTable * binning (scalar * targets, double eps) {
  scalar * tnorm = list_clone (targets);
  foreach()
    for (size_t i = 0; i < list_len (targets); i++) {
      scalar tn = tnorm[i], t = targets[i];
      tn[] = t[];
    }
  binning_normalize (tnorm);
  BinTable * table = binning_build_table (tnorm, eps);
  binning_average (table, tnorm);
  free (tnorm);
  return table;
}

/**
## References

~~~bib
@article{goldin2009cell,
  title={A cell agglomeration algorithm for accelerating detailed chemistry in CFD},
  author={Goldin, Graham M and Ren, Zhuyin and Zahirovic, Selma},
  journal={Combustion Theory and Modelling},
  volume={13},
  number={4},
  pages={721--739},
  year={2009},
  publisher={Taylor \& Francis}
}
~~~
*/
