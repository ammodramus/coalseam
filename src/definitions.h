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
