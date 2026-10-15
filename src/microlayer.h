/**
# Microlayer post-processing

We implement the post-processing of microlayer features in nucleate boiling
simulations. We assume that the simulation is axial-symmetric, and that the
bubble is placed on the lower-left corner of the domain, with a solid boundary
on the left. */

typedef struct {
  // bubble volume
  double volume;
  // bubble equivalent diameter
  double deq;
  // velocity of the contact line
  double ucl;
  // radius of the contact line
  double rcl;
  // radius of the bubble (nose)
  double rb;
  // radius of the microlayer
  double rml;
  // microlayer thickness @ theta degrees
  double delta;
  // average microlayer thickness up to theta degrees
  double avgdelta;
  // total vaporization rate microlayer
  double fluxml;
  // total vaporization rate bulk
  double fluxbulk;
} mlstats;

static double interface_position (Point point, scalar f, bool boundary = false) {
  if (f.height.y.i) {
    vector h = f.height;
    if (boundary && h.y[] != nodata && h.y[-1] != nodata)
      return y + Delta*0.5*(height (h.y[]) + height (h.y[-1]));
    else if (h.y[] != nodata)
      return y + Delta*height (h.y[]);
    else {
      coord m = interface_normal (point, f), p;
      double alpha = plane_alpha (f[], m);
      plane_area_center (m, alpha, &p);
      return y + Delta*p.y;
    }
  }
  return nodata;
}

static inline
int is_contact (Point point, scalar f) {
  extern double theta0;
  if (f[] == 0. || f[] == 1.)
    return false;
  else {
    for (int i = -1; i <= 1; i += 2)
      if (theta0 < 90. && f[0,i] == 0. && is_boundary (neighbor(-1)))
        return true;
      else if (theta0 >= 90. && f[0,i] == 1. && is_boundary (neighbor(-1)))
        return true;
  }
  return false;
}

static inline
int is_adjacent (Point point, scalar f) {
  if (f[] == 0. || f[] == 1.)
    return false;
  else {
    for (int i = -1; i <= 1; i += 2)
      if (is_contact (neighborp(0,i), f))
        return true;
  }
  return false;
}

scalar deg[], region[];
enum regions { IS_MICROLAYER = 1, IS_BULK };

trace
mlstats microlayer_stats (scalar f, double theta) {
  mlstats ml = {0};

  // bubble equivalent diameter
  double volume = 0.;
  foreach (reduction(+:volume))
    volume += (1. - f[])*dv();
  ml.volume = volume;
  ml.deq = 2.*cbrt (3./2.*volume);

  // radius of the contact line
  double rcl = 0.;
  foreach_boundary (left, reduction(max:rcl))
    if (is_contact (point, f))
      rcl = max (rcl, interface_position (point, f, true));
  ml.rcl = rcl;

  // radius of the bubble nose
  double rb = 0.;
  foreach (reduction(max:rb))
    if (f[] > F_ERR && f[] < 1. - F_ERR)
      rb = max (rb, interface_position (point, f, false));
  ml.rb = rb;

  // Orientation of the interface in degrees: this is a strange definition of
  // quadrants, where the microlayer might have positive or negative values,
  // while every point with negative n.x is increased by 180 deg. This makes it
  // easier to distinguish different regions
  //scalar deg[];
  foreach() {
    deg[] = nodata;
    if (f[] > F_ERR && f[] < 1. - F_ERR) {
      coord m = interface_normal (point, f);
      normalize (&m);
      foreach_dimension()
        m.x *= -1;
      deg[] = 90. - acos (m.y)*180./pi;
      if (m.x > 0.)
        deg[] += 180.;
    }
  }

  // radius of the microlayer
  double rml = 0., delta = 0.;
  foreach (reduction(max:rml) reduction(max:delta))
    if (deg[] != nodata && deg[] > 0. && deg[] <= theta) {
      rml = max (rml, interface_position (point, f, false));
      vector h = f.height;
      delta = max (delta, x + Delta*height (h.x[]));
    }
  ml.rml = rml;
  ml.delta = delta;

  // velocity of the contact line
  double half = 0.5*L0/(1 << grid->maxdepth);
  ml.ucl = interpolate (u.y, half, rcl);

  // Divide the interface in two regions (microlayer and bulk)
  //scalar region[];
  foreach() {
    region[] = nodata;
    if (deg[] != nodata)
      region[] = IS_BULK;
    if (y <= rml && deg[] < 180.)
      region[] = IS_MICROLAYER;
  }

  // Compute mass flux in the different regions
  double fluxml = 0., fluxbulk = 0.;
  foreach (reduction(+:fluxml) reduction(+:fluxbulk)) {
    if (region[] != nodata) {
      coord m = interface_normal (point, f), p;
      double alpha = plane_alpha (f[], m);
      double area = plane_area_center (m, alpha, &p);
#if AXI
      double dirac = area*(y + p.y*Delta)/(Delta*y)*cm[];
#else
      double dirac = area/Delta*cm[];
#endif
      if (region[] == IS_MICROLAYER)
        fluxml += mEvapTot[]*dirac*dv();

      if (region[] == IS_BULK)
        fluxbulk += mEvapTot[]*dirac*dv();
    }
  }
  ml.fluxml = fluxml;
  ml.fluxbulk = fluxbulk;

  // Average microlayer thickness
  int count = 0;
  double avgdelta = 0.;
  foreach (reduction(+:count) reduction(+:avgdelta))
    if (region[] == IS_MICROLAYER) {
      coord m = interface_normal (point, f), p;
      double alpha = plane_alpha (f[], m);
      plane_area_center (m, alpha, &p);
      avgdelta += x + p.x*Delta;
      count++;
    }
  ml.avgdelta = (count > 0) ? avgdelta/count : 0.;

  return ml;
}

void microlayer_print (FILE * fp, mlstats ml, bool header = false) {
  const int W = 20;
  if (header) {
    fprintf(fp,
        "%-*s %-*s %-*s %-*s %-*s %-*s %-*s %-*s %-*s %-*s %-*s\n",
        W, "time(1) [ms]", W, "volume(2) [m3]", W, "Deq(3) [mm]", W, "Ca(4) [-]",
        W, "Rcl(5) [mm]", W, "Rb(6) [mm]", W, "Rml(7) [mm]", W, "delta(8) [µm]",
        W, "fluxml(9) [kg/s]", W, "fluxbulk(10) [kg/s]", W, "avgdelta(11) [µm]");
  }

  fprintf(fp,
      "%-*g %-*g %-*g %-*g %-*g %-*g %-*g %-*g %-*g %-*g %-*g\n",
      W, t*1e3,
      W, ml.volume,
      W, ml.deq*1e3,
      W, ml.ucl*mu1/f.sigma,
      W, ml.rcl*1e3,
      W, ml.rb*1e3,
      W, ml.rml*1e3,
      W, ml.delta*1e6,
      W, ml.fluxml,
      W, ml.fluxbulk,
      W, ml.avgdelta*1e6);
}

