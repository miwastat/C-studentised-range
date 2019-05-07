/*
 *  This program tabulates the upper quantiles
 *    of the Studentised range distribution.
 *
 *  command format: strng_tbl k_end alpha [index]
 *
 *  Arguments
 *    k_end    : k = 2, ..., k_end.
 *               If k_end > 100,
 *               k = 2, ..., 20, 50, 100, 200, 500, 1000.
 *    alpha    : upper probability
 *    [index]  : If index==2, nu runs from 1 to 40.
 *
 *  Required functions:
 *    double strng_lq()
 *    static void line(int i)
 *
 *  Include files:
 *    <stdio.h>
 *    <stdlib.h>
 *    <math.h>
 *
 *  Note
 *    The tables can be stored in a file by piping such as
 *      ./strng_tbl 20 0.05 2 > strng05.txt
 *
 *  Stored in:
 *    strng_tbl.c
 *
 *  History
 *    2002-10-04: Last modified.
 *    2018-11-10: Created for the new version.
 *    2019-04-26: k_end > 100
 *
 *  Coded by Tetsuhisa Miwa.
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define EPS (1.0e-8)

double strng_lq(double p, int k, int nu,
                double xeps, double peps, int *itr);

static void line(int i)
{
  for( ; i > 0; i--)
    printf("-");
  printf("\n");
}

int main(int argc, char **argv)
{
  double  alpha, xeps, peps, q;
  int     kupper[5]={50, 100, 200, 500, 1000}, k[99], ke, j;
  int     index=1, nu[106], i, itr, itrmax=0;
  FILE    *fout;

  if(argc < 3) {
    printf("command format: strng_tbl"
           " k_end alpha [index=1 or 2]\n");
    exit(1);
  }

  ke = atoi(argv[1]) - 2; // end value of k
  if(ke <= 98) {
    for(j=0; j <= ke; j++)
      k[j] = j+2;
  } else {
    ke = 23;
    for(j=0; j <= 18; j++)
      k[j] = j+2;
    for( ; j <= ke; j++)
      k[j] = kupper[j-19];
  }

  alpha = atof(argv[2]);
  xeps = EPS;
  peps = alpha*EPS;

  if(argc >= 4) {
    index = atoi(argv[3]);
    if(index != 1)  // index value should be 1 or 2
      index = 2;
  }
  for(i=0; i < 20*index; i++)
    nu[i] = i+1;
  for(i=0; i < 5; i++)
    nu[i+20*index] = 120*index/(5-i);
  nu[5+20*index] = 0;

  printf("The Studentised range upper quantiles "
          "q(k, nu; %4.2lf)\n", alpha);
  line(7*ke + 12);
  printf(" nu  k->%3i", k[0]);
  for(j=1; j <= ke; j++)
    printf("%7i", k[j]);
  printf("\n");
  line(7*ke + 12);

  for(i=0; i < 6+20*index; i++){
    if(nu[i] == 0)
      printf("Inf  ");
    else
      printf("%3i  ", nu[i]);

    for(j=0; j <= ke; j++){
      q = strng_lq(1.0-alpha, k[j], nu[i], xeps, peps, &itr);
      if(q < 100.0)
        printf("%7.3lf", q);
      else
        printf("%7.2lf", q);
      if(itr > itrmax)
        itrmax = itr;
    }
    printf("\n");

    if((i+1)%10==0)
      line(7*ke+12);
    if((i+1)==20 && index==2){
      printf(" nu  k->%3i", k[0]);
      for(j=1; j <= ke; j++)
        printf("%7i", k[j]);
      printf("\n");
      line(7*ke+12);
    }
  }
  line(7*ke+12);

  printf("max.iterations = %5i\n", itrmax);
}
