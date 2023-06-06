/**
# Pin a Liquid Droplet at a Specific Point

This module defines the boundary conditions required to
pin the droplet a specific point of the domain.
*/

#include "contact.h"

vector h[];

/**
Known the position of the interface $a$ we compute the corresponding
heigth function at the neighbouring cell.*/

foreach_dimension()
static double line_x (Point point, scalar h, double a)
{
  return (h[] == nodata ? nodata : -h[] + 2.*(a-x)/Delta);
}

#define contact_line(theta) contact_line_ (point, neighbor, _s, theta)

double contact_line_ (Point point, Point neighbor, scalar h, double a)
{
  if (neighbor.i != point.i)
    return line_x (point, h, a);
  if (neighbor.j != point.j)
    return line_y (point, h, a);
  assert (false); // not reached
  return 0.;
}

double ap = 0.;
