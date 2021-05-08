/*
 *  Test program for strng_lq().
 *    Command format: ./strng_lq_tst k nu alpha xeps
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern double strng_lq(double p, int k, int nu,
                       double xeps, double peps, int *itr);

int main(int argc, char **argv)
{
  int k, nu, itr;
  double x, x0, x1, alpha, xeps, peps;

  if(argc < 5) {
    printf("Command format: srng_lq_tst k nu alpha xeps\n");
    exit (1);
  }
  k = atoi(argv[1]);
  nu = atoi(argv[2]);
  alpha = atof(argv[3]);
  xeps = atof(argv[4]);
  peps = alpha*xeps;

  x = strng_lq(1.0-alpha, k, nu, xeps, peps, &itr);
  printf("itr = %4d, quantile = %24.16g\n", itr, x);

  // Interpolation between nu=240 and infinity.
  if(nu > 240)
    {
      x0 = strng_lq(1.0-alpha, k, 0, xeps, peps, &itr);
      x1 = strng_lq(1.0-alpha, k, 240, xeps, peps, &itr);
      x = x0 + (240.0/nu) * (x1-x0);
      printf("Interpolation in 1/nu\n"
             "itr = %4d, quantile = %24.16g\n", itr, x);
    }
  exit (0);
}
