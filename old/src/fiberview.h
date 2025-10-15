/**
# *colorbar()*: draws a colorbar.

* *map*: the colormap.
* *size*: the size.
* *pos*: the relative screen position (default is lower-left corner).
* *label*: a label.
* *lscale*: the label scale factor.
* *min*, *max*: the range.
* *horizontal*: true for anb horizontal colorbar.
* *border*: adds a border.
* *mid*: adds a mid-value label.
* *lc*: the line color.
* *lw*: the line width.
* *fsize*: another size.
* *format*: how to format numbers.
* *levels*: the number of subdivisions.

Note that this cannot be used with [jview](jview/README) (yet) because
it mixes surface and line rendering. */

trace
bool draw_solid (scalar Z, double min, double max,
    Colormap map = jet, float lc[3] = {0}, float lw = 1., int levels = 50)
{
  bview * view = draw();

  double x0 = -0.5*L0, y0 = 0., th = Y0;
  double dh = L0/(double)levels, xi = x0, yi = y0;

  double cmap[NCMAP][3];
  (* map) (cmap);
  glBegin (GL_QUADS);
  for (int i = 0; i < levels; i++) {
    double Zh = 0.;
    foreach_point (xi, Y0)
      Zh = Z[];
    Color c = colormap_color (cmap, Zh, min, max);
    glColor3f (c.r/255., c.g/255., c.b/255.);
    glVertex2d (xi, yi);
    glVertex2d (xi + dh, yi);
    glVertex2d (xi + dh, yi + th);
    glVertex2d (xi, yi + th);
    xi = x0 + (i+1)*dh;
    view->ni++;
  }
  glEnd();

  draw_lines (view, lc, lw) {
    glBegin (GL_LINE_LOOP);
    glvertex2d (view, x0, y0);
    glvertex2d (view, x0 + L0, y0);
    glvertex2d (view, x0 + L0, y0 + th);
    glvertex2d (view, x0, y0 + th);
    glEnd();
    view->ni++;
  }

  return true;
}

