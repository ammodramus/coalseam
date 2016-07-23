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

#ifndef DEFINITIONS_HEADER
#define DEFINITIONS_HEADER

#include <stdio.h>

#define PERROR(msg,...) do {fprintf(stderr,"\n\nProgram error:\n%s\n\n",msg); exit(-1);} while(0)
#define MAX(a,b) (a > b ? a : b)
#define MIN(a,b) (a > b ? b : a)
#define GETMLRGCCoalTime(MLRGC,l,i1,i2) (MLRGC).coaltimes[l][i1][i2]
#define CHECKPOINTER(ptr) if(ptr == NULL) do {fprintf(stderr,"\nFailed memory allocation: %s\n",#ptr); exit(-1);} while(0)
#define NOTININTERVAL(x, a, b) (x < a || x > b)
#define REPORTI(x) printf(#x " = %i\n", x)
#define REPORTF(x) printf(#x " = %f\n", x)

typedef int32_t twoints[2]; /* makes allocating the 3D pedigree (numgens by N [which may vary with time] by 2) easier */
typedef int8_t twoeights[2]; /* makes allocating the 3D pedigree (numgens by N [which may vary with time] by 2) easier */

#endif
