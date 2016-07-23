/*

    coalseam -- a program for simulating coalescence in pedigrees
    Copyright 2016 Peter Wilton

    This file is part of coalseam.

    coalseam is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    coalseam is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with coalseam.  If not, see <http://www.gnu.org/licenses/>.

*/

#include <stdint.h>
#include <math.h>
#include <time.h>
#include "random.h"

/* runifd() returns an int between lower and upper, inclusive. */
inline int32_t runifd(int32_t lower, int32_t upper)
{
	/* This is just as efficient as a modulo version, plus it's 
	 * 100% accurate. */ 
	return((int32_t)((runif()*(float)(upper-lower+1)))+lower);
}

void setseed(unsigned long i1,unsigned long i2)
{
	z=i1; w=i2;
}


inline int32_t rbern(double prob)
{
	if(runif()<prob)
		return 1;
	else
		return 0;
}

double localgammln(double xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

double rpois(double xm) 
{ 
	static double sq,alxm,g,oldm=(-1.0); 
	double em,t,y; 
 
	if (xm < 12.0) { 
		if (xm != oldm) { 
			oldm=xm; 
			g=exp(-xm); 
		} 
		em = -1; 
		t=1.0; 
		do { 
			++em; 
			t *= runif(); 
		} while (t > g); 
	} else { 
		if (xm != oldm) { 
			oldm=xm; 
			sq=sqrt(2.0*xm); 
			alxm=log(xm); 
			g=xm*alxm-localgammln(xm+1.0); 
		} 
		do { 
			do { 
				y=tan(PI*runif()); 
				em=sq*y+xm; 
			} while (em < 0.0); 
			em=floor(em); 
			t=0.9*(1.0+y*y)*exp(em*alxm-localgammln(em+1.0)-g); 
		} while (runif() > t); 
	} 
	return em; 
} 

/* The following is an implementation of the Ziggurat method */
/* Setting up the ziggurat tables */
static long hz;
static unsigned long iz, kn[128], ke[256];
static float wn[128],fn[128], we[256],fe[256];

/* nfix() generates variates from the residue when rejection in RNOR occurs. */

float nfix(void)
{
	const float r = 3.442620f;     /* The start of the right tail */
	static float x, y;
	 for(;;)
	  {  x=hz*wn[iz];      /* iz==0, handles the base strip */
		 if(iz==0)
		   { do{ x=-log(runif())*0.2904764; y=-log(runif());}	/* .2904764 is 1/r */
			while(y+y<x*x);
			return (hz>0)? r+x : -r-x;
		   }
							 /* iz>0, handle the wedges of other strips */
		  if( fn[iz]+runif()*(fn[iz-1]-fn[iz]) < exp(-.5*x*x) ) return x;

		 /* initiate, try to exit for(;;) for loop*/
		  hz=(znew+wnew);
		  iz=hz&127;
		  if(fabs(hz)<kn[iz]) return (hz*wn[iz]);
	  }
}

/* efix() generates variates from the residue when rejection in REXP occurs. */

float efix(void)
{ float x;
 for(;;)
  {  if(iz==0) return (7.69711-log(runif()));          /* iz==0 */
     x=z*we[iz]; if( fe[iz]+runif()*(fe[iz-1]-fe[iz]) < exp(-x) ) return (x);

      /* initiate, try to exit for(;;) loop */
   z=(znew+wnew);
   iz=(z&255);
   if(z<ke[iz]) return (z*we[iz]);
  }
}
/*--------This procedure sets the seed and creates the tables------*/

void zigset(unsigned long wseed)
{  
	const double m1 = 2147483648.0, m2 = 4294967296.;
	double dn=3.442619855899,tn=dn,vn=9.91256303526217e-3, q;
	double de=7.697117470131487, te=de, ve=3.949659822581572e-3;
	int i;
	w^=wseed;

/* Set up tables for rnorm() */
   
   q=vn/exp(-.5*dn*dn);
   kn[0]=(dn/q)*m1;
   kn[1]=0;

   wn[0]=q/m1;
   wn[127]=dn/m1;

   fn[0]=1.;
   fn[127]=exp(-.5*dn*dn);

    for(i=126;i>=1;i--)
    {dn=sqrt(-2.*log(vn/dn+exp(-.5*dn*dn)));
     kn[i+1]=(dn/tn)*m1;
     tn=dn;
     fn[i]=exp(-.5*dn*dn);
     wn[i]=dn/m1;
    }
	
/* Set up tables for rexp() */
    q = ve/exp(-de);
    ke[0]=(de/q)*m2;
    ke[1]=0;

    we[0]=q/m2;
    we[255]=de/m2;

    fe[0]=1.;
    fe[255]=exp(-de);

	for(i=254;i>=1;i--)
    {
	de=-log(ve/de+exp(-de));
	ke[i+1]= (de/te)*m2;
	te=de;
	fe[i]=exp(-de);
	we[i]=de/m2;
    }
}

/* Just call randseed in main() to seed the generator. */

void randseed(void)
{
	int i;
	z=getpid();
	w=time(0);
	/* By default, ziggurat tables are set up. */
	zigset((unsigned long)z*w);
	// Warm up the RNG...
	for(i = 0; i < 500; i++)
		randint32();
	
}

/* The following three functions are taken directly from R. See random.c in
 * R-[version]/src/main/. */


/* This particular function uses the function revsort, which sorts an array of
 * doubles into reverse (decreasing order), so that on average, the most
 * probably outcome is first to be looked up. */

/*
void unequal_sample_replace(int n, double *p, int *perm, int nans, int *ans)
{
    double rU;
    int i, j;
    int nm1 = n - 1;

    // record element identities 
    for (i = 0; i < n; i++)
	perm[i] = i + 1;

    // sort the probabilities into descending order 
	// This is R's function, too. See sort.h (PW) 
    revsort(p, perm, n);

   	//  compute cumulative probabilities 
    for (i = 1 ; i < n; i++)
	p[i] += p[i - 1];

    // compute the sample 
    for (i = 0; i < nans; i++)
	{
		rU = runif();
		for (j = 0; j < nm1; j++) 
		{
			if (rU <= p[j])
			break;
		}
		ans[i] = perm[j];
    }
}
*/


/* Samples {0, 1, ..., n-1} without replacement.
 *
 * k   number of elements to sample
 * n   number of elements in pool
 * y   vector (length k) of 
 * x   bucket vector (length n) for holding possibilities 
 */
void equal_sample_noreplace(int32_t k, int32_t n, int32_t * y, int32_t * x)
{
    int32_t i, j;
    for (i = 0; i < n; i++)
		x[i] = i;
    for (i = 0; i < k; i++)
	{
		j = (int32_t)((double)n * runif());
		y[i] = x[j];
		x[j] = x[--n];
    }
}

/* Another version of the sampling {0, 1, ..., n-1} without replacement, but
 * more useful for big n and small k, (e.g., sampling 2 elements from {0, 1,
 * ..., 100000}), where it's improbable that two elements would be sampled more
 * than once.
 *
 * k          number of elements to sample
 * n          number of elements in pool
 * toreturn   output vector (length k)
 */
void equal_sample_noreplace_smallkbign(int k, int n, int32_t * toreturn)
{
	int32_t i,j;
	for(i = 0; i < k; i++)
	{
		toreturn[i] = runifd(0,n-1);
		for(j = 0; j < i; j++)
		{
			if(toreturn[i] == toreturn[j])
			{
				toreturn[i] = runifd(0,n-1);
				j = -1;		// Sets the clock back to zero at the next time through the loop.
			}
		}
	}
	return;
}

/* Samples Uniformly from {0, ..., n-1} with replacement
 *
 * k   number of elements to sample
 * n   number of elements in pool
 * y   output vector (length n) 
 */
void equal_sample_replace(int k, int n, int * y)
{
    int i;
    for (i = 0; i < k; i++)
		y[i] = (int)((double)n * runif());
}

void revsort(double * a, int * ib, int n)
{
/* Sort a[] into descending order by "heapsort";
 * sort ib[] alongside;
 * if initially, ib[] = 1...n, it will contain the permutation finally
 */

    int l, j, ir, i;
    double ra;
    int ii;

    if (n <= 1) return;

    a--; ib--;

    l = (n >> 1) + 1;
    ir = n;

    for (;;) {
	if (l > 1) {
	    l = l - 1;
	    ra = a[l];
	    ii = ib[l];
	}
	else {
	    ra = a[ir];
	    ii = ib[ir];
	    a[ir] = a[1];
	    ib[ir] = ib[1];
	    if (--ir == 1) {
		a[1] = ra;
		ib[1] = ii;
		return;
	    }
	}
	i = l;
	j = l << 1;
	while (j <= ir) {
	    if (j < ir && a[j] > a[j + 1]) ++j;
	    if (ra > a[j]) {
		a[i] = a[j];
		ib[i] = ib[j];
		j += (i = j);
	    }
	    else
		j = ir + 1;
	}
	a[i] = ra;
	ib[i] = ii;
    }
}

void ProbSampleReplace(int n, double *p, int *perm, int nans, int *ans)
{
    double rU;
    int i, j;
    int nm1 = n - 1;

    /* record element identities */
    for (i = 0; i < n; i++)
		perm[i] = i;

    /* sort the probabilities into descending order */
    revsort(p, perm, n);

    /* compute cumulative probabilities */
    for (i = 1 ; i < n; i++)
		p[i] += p[i - 1];

    /* compute the sample */
    for (i = 0; i < nans; i++)
	{
		rU = runif();
		for (j = 0; j < nm1; j++) 
		{
			if (rU <= p[j])
			break;
		}
		ans[i] = perm[j];
    }
}

/* Binomial distribution taken from gsl. */


#define SMALL_MEAN 14           /* If n*p < SMALL_MEAN then use BINV
                                   algorithm. The ranlib
                                   implementation used cutoff=30; but
                                   on my computer 14 works better */

#define BINV_CUTOFF 110         /* In BINV, do not permit ix too large */

#define FAR_FROM_MEAN 20        /* If ix-n*p is larger than this, then
                                   use the "squeeze" algorithm.
                                   Ranlib used 20, and this seems to
                                   be the best choice on my machine as
                                   well */

//#define LNFACT(x) gsl_sf_lnfact(x)
inline double Stirling (double y1)
{
  double y2 = y1 * y1;
  double s =
    (13860.0 -
     (462.0 - (132.0 - (99.0 - 140.0 / y2) / y2) / y2) / y2) / y1 / 166320.0;
  return s;
}

double gsl_pow_uint(double x, unsigned int n)
{
  double value = 1.0;

  /* repeated squaring method 
   * returns 0.0^0 = 1.0, so continuous in x
   */
  do {
     if(n & 1) value *= x;  /* for n odd */
     n >>= 1;
     x *= x;
  } while (n);

  return value;
}

unsigned int gsl_ran_binomial (unsigned int n, double p)
{
  int ix;                       /* return value */
  int flipped = 0;
  double q, s, np;

  if (n == 0)
    return 0;

  if (p > 0.5)
    {
      p = 1.0 - p;              /* work with small p */
      flipped = 1;
    }

  q = 1 - p;
  s = p / q;
  np = n * p;

  /* Inverse cdf logic for small mean (BINV in K+S) */

  if (np < SMALL_MEAN)
    {
      double f0 = gsl_pow_uint (q, n);   /* f(x), starting with x=0 */

      while (1)
        {
          /* This while(1) loop will almost certainly only loop once; but
           * if u=1 to within a few epsilons of machine precision, then it
           * is possible for roundoff to prevent the main loop over ix to
           * achieve its proper value.  following the ranlib implementation,
           * we introduce a check for that situation, and when it occurs,
           * we just try again.
           */

          double f = f0;
          double u = runif();

          for (ix = 0; ix <= BINV_CUTOFF; ++ix)
            {
              if (u < f)
                goto Finish;
              u -= f;
              /* Use recursion f(x+1) = f(x)*[(n-x)/(x+1)]*[p/(1-p)] */
              f *= s * (n - ix) / (ix + 1);
            }

          /* It should be the case that the 'goto Finish' was encountered
           * before this point was ever reached.  But if we have reached
           * this point, then roundoff has prevented u from decreasing
           * all the way to zero.  This can happen only if the initial u
           * was very nearly equal to 1, which is a rare situation.  In
           * that rare situation, we just try again.
           *
           * Note, following the ranlib implementation, we loop ix only to
           * a hardcoded value of SMALL_MEAN_LARGE_N=110; we could have
           * looped to n, and 99.99...% of the time it won't matter.  This
           * choice, I think is a little more robust against the rare
           * roundoff error.  If n>LARGE_N, then it is technically
           * possible for ix>LARGE_N, but it is astronomically rare, and
           * if ix is that large, it is more likely due to roundoff than
           * probability, so better to nip it at LARGE_N than to take a
           * chance that roundoff will somehow conspire to produce an even
           * larger (and more improbable) ix.  If n<LARGE_N, then once
           * ix=n, f=0, and the loop will continue until ix=LARGE_N.
           */
        }
    }
  else
    {
      /* For n >= SMALL_MEAN, we invoke the BTPE algorithm */

      int k;

      double ffm = np + p;      /* ffm = n*p+p             */
      int m = (int) ffm;        /* m = int floor[n*p+p]    */
      double fm = m;            /* fm = double m;          */
      double xm = fm + 0.5;     /* xm = half integer mean (tip of triangle)  */
      double npq = np * q;      /* npq = n*p*q            */

      /* Compute cumulative area of tri, para, exp tails */

      /* p1: radius of triangle region; since height=1, also: area of region */
      /* p2: p1 + area of parallelogram region */
      /* p3: p2 + area of left tail */
      /* p4: p3 + area of right tail */
      /* pi/p4: probability of i'th area (i=1,2,3,4) */

      /* Note: magic numbers 2.195, 4.6, 0.134, 20.5, 15.3 */
      /* These magic numbers are not adjustable...at least not easily! */

      double p1 = floor (2.195 * sqrt (npq) - 4.6 * q) + 0.5;

      /* xl, xr: left and right edges of triangle */
      double xl = xm - p1;
      double xr = xm + p1;

      /* Parameter of exponential tails */
      /* Left tail:  t(x) = c*exp(-lambda_l*[xl - (x+0.5)]) */
      /* Right tail: t(x) = c*exp(-lambda_r*[(x+0.5) - xr]) */

      double c = 0.134 + 20.5 / (15.3 + fm);
      double p2 = p1 * (1.0 + c + c);

      double al = (ffm - xl) / (ffm - xl * p);
      double lambda_l = al * (1.0 + 0.5 * al);
      double ar = (xr - ffm) / (xr * q);
      double lambda_r = ar * (1.0 + 0.5 * ar);
      double p3 = p2 + c / lambda_l;
      double p4 = p3 + c / lambda_r;

      double var, accept;
      double u, v;              /* random variates */

    TryAgain:

      /* generate random variates, u specifies which region: Tri, Par, Tail */
      u = runif () * p4;
      v = runif ();

      if (u <= p1)
        {
          /* Triangular region */
          ix = (int) (xm - p1 * v + u);
          goto Finish;
        }
      else if (u <= p2)
        {
          /* Parallelogram region */
          double x = xl + (u - p1) / c;
          v = v * c + 1.0 - fabs (x - xm) / p1;
          if (v > 1.0 || v <= 0.0)
            goto TryAgain;
          ix = (int) x;
        }
      else if (u <= p3)
        {
          /* Left tail */
          ix = (int) (xl + log (v) / lambda_l);
          if (ix < 0)
            goto TryAgain;
          v *= ((u - p2) * lambda_l);
        }
      else
        {
          /* Right tail */
          ix = (int) (xr - log (v) / lambda_r);
          if (ix > (double) n)
            goto TryAgain;
          v *= ((u - p3) * lambda_r);
        }

      /* At this point, the goal is to test whether v <= f(x)/f(m) 
       *
       *  v <= f(x)/f(m) = (m!(n-m)! / (x!(n-x)!)) * (p/q)^{x-m}
       *
       */

      /* Here is a direct test using logarithms.  It is a little
       * slower than the various "squeezing" computations below, but
       * if things are working, it should give exactly the same answer
       * (given the same random number seed).  */

      /* More efficient determination of whether v < f(x)/f(M) */

      k = abs (ix - m);

      if (k <= FAR_FROM_MEAN)
        {
          /* 
           * If ix near m (ie, |ix-m|<FAR_FROM_MEAN), then do
           * explicit evaluation using recursion relation for f(x)
           */
          double g = (n + 1) * s;
          double f = 1.0;

          var = v;

          if (m < ix)
            {
              int i;
              for (i = m + 1; i <= ix; i++)
                {
                  f *= (g / i - s);
                }
            }
          else if (m > ix)
            {
              int i;
              for (i = ix + 1; i <= m; i++)
                {
                  f /= (g / i - s);
                }
            }

          accept = f;
        }
      else
        {
          /* If ix is far from the mean m: k=ABS(ix-m) large */

          var = log (v);

          if (k < npq / 2 - 1)
            {
              /* "Squeeze" using upper and lower bounds on
               * log(f(x)) The squeeze condition was derived
               * under the condition k < npq/2-1 */
              double amaxp =
                k / npq * ((k * (k / 3.0 + 0.625) + (1.0 / 6.0)) / npq + 0.5);
              double ynorm = -(k * k / (2.0 * npq));
              if (var < ynorm - amaxp)
                goto Finish;
              if (var > ynorm + amaxp)
                goto TryAgain;
            }

          /* Now, again: do the test log(v) vs. log f(x)/f(M) */

//#if USE_EXACT
          /* This is equivalent to the above, but is a little (~20%) slower */
          /* There are five log's vs three above, maybe that's it? */

          //accept = LNFACT (m) + LNFACT (n - m)
            //- LNFACT (ix) - LNFACT (n - ix) + (ix - m) * log (p / q);

//#else /* USE STIRLING */
          /* The "#define Stirling" above corresponds to the first five
           * terms in asymptoic formula for
           * log Gamma (y) - (y-0.5)log(y) + y - 0.5 log(2*pi);
           * See Abramowitz and Stegun, eq 6.1.40
           */

          /* Note below: two Stirling's are added, and two are
           * subtracted.  In both K+S, and in the ranlib
           * implementation, all four are added.  I (jt) believe that
           * is a mistake -- this has been confirmed by personal
           * correspondence w/ Dr. Kachitvichyanukul.  Note, however,
           * the corrections are so small, that I couldn't find an
           * example where it made a difference that could be
           * observed, let alone tested.  In fact, define'ing Stirling
           * to be zero gave identical results!!  In practice, alv is
           * O(1), ranging 0 to -10 or so, while the Stirling
           * correction is typically O(10^{-5}) ...setting the
           * correction to zero gives about a 2% performance boost;
           * might as well keep it just to be pendantic.  */

          {
            double x1 = ix + 1.0;
            double w1 = n - ix + 1.0;
            double f1 = fm + 1.0;
            double z1 = n + 1.0 - fm;

            accept = xm * log (f1 / x1) + (n - m + 0.5) * log (z1 / w1)
              + (ix - m) * log (w1 * p / (x1 * q))
              + Stirling (f1) + Stirling (z1) - Stirling (x1) - Stirling (w1);
          }
        }


      if (var <= accept)
        {
          goto Finish;
        }
      else
        {
          goto TryAgain;
        }
    }

Finish:

  return (flipped) ? (n - ix) : (unsigned int)ix;
}

int32_t gsl_ran_binomial_tpe(int32_t n, double p)
{
  return gsl_ran_binomial(n, p);
}

int32_t is_in_int(int32_t testNum, int32_t * array, int32_t arrayLen)
{
	static int i;
	for(i = 0; i < arrayLen; i++)
	{
		if(testNum == array[i])
			return 1;
	}
	return 0;
}

/* rgeom() adapted from gsl_ran_geometric in GSL library
 * https://fossies.org/dox/gsl-1.16/randist_2geometric_8c_source.html
 *
 * the returned value is the number of trials until success, where p is the
 * probability of success.
 *
 *  pmf: 
 *    p(k) =  p(1-p)^(k-1)
 *
 *
 * (note that the pmf of the R geometric distn is defined differently, with 
 *    p(k) = p(1-p)^k,
 *  the number of failures before success)
 * 
 */

int32_t rgeom(const double p)
{
    double u = runif();

    int32_t k;

    if (p == 1)
        k = 1;

    else
        k = (int32_t)(log (u) / log (1 - p) + 1);

    return k;
}

