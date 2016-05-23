#ifndef INTVECTOR_HEADER
#define INTVECTOR_HEADER

#include <stdint.h>

typedef struct intvector
{
	int32_t * vector;
	int32_t curIdx;
	int32_t curSize;
} IntVector;

void IntVector_init(IntVector * iv, int32_t initSize);
void IntVector_resize(IntVector * iv, int32_t newSize);
void IntVector_add(IntVector * iv, int32_t toAdd);
void IntVector_add_multiple(IntVector * iv, int32_t * toAdd, int32_t numToAdd);
void IntVector_reset(IntVector * iv, int32_t newSize);
void IntVector_print(IntVector * iv, FILE * fout);
int32_t IntVector_max(IntVector * iv);
int32_t IntVector_sum(IntVector * iv);
double IntVector_mean(IntVector * iv);

#endif
