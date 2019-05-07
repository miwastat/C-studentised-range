/*
 *  double strng_lp(double q, int k, int nu)
 *    returns lower probability
 *    of the Studentised range distribution.
 *
 *  Arguments
 *    q:  Studentised range value
 *    k:  number of treatments
 *    nu: error degrees of freedom (nu<=0 means nu=infinity)
 *
 *  Required functions
 *    extern double rng_lp()
 *    static double chi2u()
 *    static double chi2l()
 *    static double coef()
 *    static double rlower()
 *    static double f()
 *
 *  Include files
 *    <math.h>
 *
 *  Note
 *    1) Uses the 40-node Legendre-Gauss quadrature.
 *    2) The accuracy is of order e-11 or more (I hope).
 *    3) This accuracy is not guaranteed for k > 1000.
 *    4) Integrates twice if ru/q < su (ru: upper limit of range).
 *
 *  Stored in
 *   strng_lp.c
 *
 *  History
 *    c. 1994:    First written in Fortran.
 *    2018-11-02: Created for the new version.
 *
 *  License
 *    GPLv3 (Free and No Warranty)
 *    https://www.gnu.org/licenses/
 *
 *  Coded by Tetsuhisa Miwa.
 */


#include <math.h>
#define LOGSQRTPI 0.572364942924700087071713675676529356  // log(sqrt(pi))

extern double rng_lp(double r, int k);

/* Upper limit for chi^2(nu) with approx upper prob=0.5e-13.
 */
static double chi2u(int nu)
{
  double  first[5]={56.73, 61.26, 65.01, 68.38, 71.50};
  double  w, z, dnu=2.0/9.0/nu;

  if(nu <= 5)
    return(first[nu-1]);
  if(nu <= 20)
    w = 7.391 - 3.050/nu + 5.208/(nu*nu);
  else
    w = 7.441 - 5.209/nu + 29.27/(nu*nu);

  // Wilson-Hilferty approximation.
  z = nu*pow(w*sqrt(dnu) + (1.0 - dnu), 3.0);
  return(z);
}

/* Lower limit for chi^2(nu) with approx lower prob=0.5e-13.
 */
static double chi2l(int nu)
{
  double  first[5]={3.926e-27, 1.0e-13, 3.281e-09, 6.324e-07, 1.546e-05};
  double  w, z, dnu=2.0/9.0/nu;

  if(nu <= 5)
    return(first[nu-1]);
  if(nu <= 20) {
    // Log approximation.
    w = -8.645 - 70.72/nu + 77.47/(nu*nu);
    z = nu*exp(w/sqrt(0.5*nu)-1.0/nu);
  }
  else {
    // Wilson-Hilferty approximation.
    w = -7.451 + 10.07/nu + 82.83/(nu*nu);
    z = nu*pow(w*sqrt(dnu) + (1.0 - dnu), 3.0);
  }    
  return(z);
}

/* Coefficient of chi distribution (Note: not chi^2 distribution).
 */
static double coef(int nu)
{
  int n;
  double g = (nu%2 == 1) ? LOGSQRTPI : 0.0;

  for(n=nu-2; n > 0; n -= 2)
    g += log(0.5*n);

  return (2.0 * exp(0.5*nu*(log(0.5*nu)-1.0) - g));
}

/* Lower limit of range with approx lower prob=0.5e-13.
 */
static double rlower(int k)
{
  double w=log(k);

  if(k <= 4)
    return(0.0);
  else if(k <= 40)
    return(exp(1.008 + 7.927/w - 46.97/(w*w) + 21.93/(w*w*w)));
  else
    return(exp(2.7404 - 10.144/w + 16.373/(w*w) - 52.728/(w*w*w)));
}

/* Integrand function
 */
static double f(double s, int k, int nu, double q, int isw)
{
  double y = exp((nu-1.0)*log(s) + 0.5*nu*(1.0-s*s));
  if(isw == 0)
    return (y);
  else
    return (y*rng_lp(s*q, k));
}


double strng_lp(double q, int k, int nu)
{
  // 40 nodes and weights for Gauss-Legendre quadrature.
  const double nd[20]={
    0.998237709710559200349622702420586492,
    0.990726238699457006453054352221372155,
    0.977259949983774262663370283712903807,
    0.957916819213791655804540999452759285,
    0.932812808278676533360852166845205716,
    0.902098806968874296728253330868493104,
    0.865959503212259503820781808354619964,
    0.824612230833311663196320230666098774,
    0.778305651426519387694971545506494848,
    0.727318255189927103280996451754930549,
    0.671956684614179548379354514961494110,
    0.612553889667980237952612450230694877,
    0.549467125095128202075931305529517970,
    0.483075801686178712908566574244823005,
    0.413779204371605001524879745803713683,
    0.341994090825758473007492481179194310,
    0.268152185007253681141184344808596183,
    0.192697580701371099715516852065149895,
    0.116084070675255208483451284408024114,
    0.0387724175060508219331934440246232947
  };
  const double wt[20]={
    0.00452127709853319125847173287818533273,
    0.0104982845311528136147421710672796524,
    0.0164210583819078887128634848823639273,
    0.0222458491941669572615043241842085732,
    0.0279370069800234010984891575077210773,
    0.0334601952825478473926781830864108490,
    0.0387821679744720176399720312904461623,
    0.0438709081856732719916746860417154958,
    0.0486958076350722320614341604481463881,
    0.0532278469839368243549964797722605046,
    0.0574397690993915513666177309104259856,
    0.0613062424929289391665379964083985959,
    0.0648040134566010380745545295667527300,
    0.0679120458152339038256901082319239860,
    0.0706116473912867796954836308552868324,
    0.0728865823958040590605106834425178359,
    0.0747231690579682642001893362613246732,
    0.0761103619006262423715580759224948230,
    0.0770398181642479655883075342838102485,
    0.0775059479784248112637239629583263270
  };
  double  sl, su, cnst, rlq, ruq, sll, x;
  double  p=0.0, p1, cntr, wdth;
  int     isw=0, i;

  if(q <= 0.0)
    return(0.0);
  // nu = infinity
  if(nu <= 0)
    return(rng_lp(q, k));

  // Upper and lower integral limits
  sl = sqrt(chi2l(nu)/nu);
  su = sqrt(chi2u(nu)/nu);
  cnst = coef(nu);

  // Lower limit of range R.
  rlq = rlower(k)/q;
  if(rlq >= su)
    return(0.0);
  if(rlq > sl)
    sl = rlq;

  // Upper limit of range R.
  ruq = (0.4494*pow(log(0.5*k), 0.87) + 10.652)/q;
  if(ruq <= sl)
    return(1.0);

  // If ru/q < su, then integrate twice:
  //   1) \int_{sl}^{ru/q}
  //   2) \int_{ru/q}^{su}, where rng_p(s*q)=1.0
  // First integrate the latter (with isw=0).
  if(ruq < su) {
    sll = sl;
    sl = ruq;
  }
  else
    isw = 1;

  for( ; isw < 2; isw++) {
    p1 = 0.0;
    cntr = 0.5*(sl+su);
    wdth = 0.5*(su-sl);
    for(i=0; i < 20; i++) {
      x = wdth*nd[i];
      p1 += wt[i] * (f(cntr-x, k, nu, q, isw) + f(cntr+x, k, nu, q, isw));
    }
    p += wdth*p1;

    if(isw == 0) {
      su = ruq;
      sl = sll;
    }
  }

  return (cnst*p);
}
