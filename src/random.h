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

#ifndef RANDOM_HEADER
#define RANDOM_HEADER

/* 
 Peter Wilton's random number generator functions
*/

#include <stdint.h>

/* Marsaglia's carry-and-multiply random number generator.
 * According to Marsaglia, has passed his own tests of randomness and is about 3/2
 * as fast as the Marsenne-Twister. Period of 2^62. 
 * (This isn't as long a period as I thought...) */

#define PI 3.141592653589793 

#define MAX_INT 4294967296   // (2^32)

#define znew  ((z=36969*(z&65535)+(z>>16))<<16)
#define wnew  ((w=18000*(w&65535)+(w>>16))&65535)

/* randint32() returns a random integer over 2^32 */
#define randint32()  (znew+wnew)

/* runif() returns a Unif r.v. (double) over [0,1) */
#define runif()  (znew+wnew)*2.328306e-10

/* rnorm() returns a standard Normal variate using the Ziggurat method */
#define rnorm() (hz=(znew+wnew), iz=hz&127, (fabs(hz) < kn[iz]) ? hz*wn[iz] : nfix())

/* rexp() returns an Exponential variate with mean 1 using the Ziggurat method */
#define rexp() (z=(znew+wnew), iz=z&255, (z < ke[iz]) ? z*we[iz] : efix())

#define RBERN(p) ((runif() < p) ? 1 : 0)

uint32_t z,w;
inline int32_t runifd(int32_t lower, int32_t upper);
void setseed(unsigned long i1,unsigned long i2);
inline int32_t rbern(double prob);
double localgammln(double xx);
double rpois(double xm) ;
float nfix();
float efix();
void zigset(unsigned long wseed);
void randseed(void);
void equal_sample_noreplace(int32_t k, int32_t n, int32_t * y, int32_t * x);
void equal_sample_noreplace_smallkbign(int k, int n, int32_t * toreturn);
void equal_sample_replace(int k, int n, int *y);
void revsort(double * a, int * ib, int n);
void ProbSampleReplace(int n, double *p, int *perm, int nans, int *ans);
inline double Stirling (double y1);
double gsl_pow_uint(double x, unsigned int n);
unsigned int gsl_ran_binomial (unsigned int n, double p);
int32_t gsl_ran_binomial_tpe(int32_t n, double p);
int32_t is_in_int(int32_t testNum, int32_t * array, int32_t arrayLen);

#endif
