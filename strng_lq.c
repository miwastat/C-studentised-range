/*
 *  double strng_lq(double p, int k, int nu,
 *                  double xeps, double peps, int *itr)
 *    returns lower quantile of
 *    the Studentised range distribution.
 *
 *  Arguments:
 *    p:    lower probability
 *    k:    number of treatments
 *    nu:   error degrees of freedom (nu<=0 means nu=infinity)
 *    xeps: precision for quantile x
 *    peps: precision for probability p
 *    *itr: number of calls of strng_lp()
 *
 *  Required functions:
 *    double strng_lp()
 *
 *  Include files:
 *    <math.h>
 *
 *  Note
 *    1) Solves the root of quadratic interpolation.
 *       References
 *         Muller, D. E. (1956).
 *         "A method for solving algebraic equations
 *         using an automatic computer",
 *         Mathematical Tables and Other Aids to Computation,
 *         Vol. 10, 208-215.
 *
 *  Stored in:
 *    strng_lq.c
 *
 *  History
 *    c. 1994:    First written in Fortran.
 *    2018-11-11: Created for the new version.
 *
 *  License
 *    GPLv3 (Free and No Warranty)
 *    https://www.gnu.org/licenses/
 *
 *  Coded by Tetsuhisa Miwa.
 */


#include  <math.h>
#define   YERR  1.0e-12 // accuracy of Studentised range probabilities

extern double strng_lp(double q, int k, int nu);


double strng_lq(double p, int k, int nu,
                double xeps, double peps, int *itr)
{
  double  x1, x2, x3, y1, y2, y3;
  double  a, b, x, y;
  int     i;

  (*itr) = 0;
  if(p <= 0.0)
    return (0.0);
  if(p >= 1.0)
    return (1.0e+99);

  // x1 < x2 (x3 <= x1 or x2 <= x3)
  // y1 < p <= y2
  x1 = 0.0;
  y1 = 0.0;
  x2 = 4.0;
  y2 = strng_lp(x2, k, nu);
  (*itr)++;
  while(y2 < p) {
    x1 = x2;
    y1 = y2;
    x2 *= 2.0;
    y2 = strng_lp(x2, k, nu);
    (*itr)++;
  }
  x3 = x2;  // (x3, y3) is used for quadratic interpolation.
  y3 = y2;

  for(i=1; i < 201; i++) {
    // bisection for odd i, or small fabs(y2-y1)
    if(i%2 == 1 || fabs(y2-y1) < YERR)
      x = 0.5*(x1+x2);

    // quadratic interpolation for even i
    else {
      if(fabs(x1-x3) < xeps || fabs(x2-x3) < xeps)
        a = 0.0;
      else
        a = ((y3-y1)/(x3-x1) - (y2-y1)/(x2-x1)) / (x3-x2);
      b = (y2-y1)/(x2-x1) - a*(x2-x1);
      if(a > 0.0)
        x = x1 + (-b + sqrt(b*b + 4.0*a*(p-y1)))/(2.0*a);
      else
        x = x1 + 2.0*(p-y1)/(b + sqrt(b*b + 4.0*a*(p-y1)));
      if(x < x1 || x > x2)
        x = 0.5*(x1+x2);
    }

    y = strng_lp(x, k, nu);
    (*itr)++;
    if(fabs(x2-x1) < xeps && fabs(y-p) < peps)
      break;

    if(y >= p) {
      x3 = x2;
      y3 = y2;
      x2 = x;
      y2 = y;
    }
    else {
      x3 = x1;
      y3 = y1;
      x1 = x;
      y1 = y;
    }
  }
  return(x);
}
