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

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include "definitions.h"
#include "intvector.h"

void IntVector_init(IntVector * iv, int32_t initSize)
{
	iv->vector = (int32_t *)malloc(sizeof(int32_t) * (size_t)initSize);
	CHECKPOINTER(iv->vector);
	iv->curIdx = 0;
	iv->curSize = initSize;
	return;
}

void IntVector_resize(IntVector * iv, int32_t newSize)
{
	iv->vector = realloc(iv->vector, sizeof(int32_t) * (size_t)newSize);
	CHECKPOINTER(iv->vector);
	iv->curSize = newSize;
	return;
}

void IntVector_add(IntVector * iv, int32_t toAdd)
{
	if(iv->curIdx+1 > iv->curSize)
		IntVector_resize(iv, iv->curSize+100);
	iv->vector[iv->curIdx] = toAdd;
	(iv->curIdx)++;
	return;
}

// Note: this is untested.
void IntVector_add_multiple(IntVector * iv, int32_t * toAdd, int32_t numToAdd)
{
	if(iv->curIdx+numToAdd > iv->curSize)
		IntVector_resize(iv, iv->curSize+2*numToAdd);
	memcpy(&(iv->vector[iv->curIdx]), toAdd, sizeof(int32_t)*(size_t)numToAdd);
	iv->curIdx += numToAdd;
	return;
}

void IntVector_reset(IntVector * iv, int32_t newSize)
{
	int32_t i;
	free(iv->vector);
	iv->vector = (int32_t *)malloc(sizeof(int32_t) * (size_t)newSize);
	CHECKPOINTER(iv->vector);
	for(i = 0; i < newSize; i++)
		iv->vector[i] = 0;
	iv->curIdx = 0;
	iv->curSize = newSize;
	return;
}

void IntVector_print(IntVector * iv, FILE * fout)
{
	int32_t i;
	for(i = 0; i < iv->curIdx; i++)
		fprintf(fout, "%i ", iv->vector[i]);
	fprintf(fout, "\n");
	return;
}

int32_t IntVector_max(IntVector * iv)
{
	int32_t i, max = -1;
	for(i = 0; i < iv->curIdx; i++)
		if(iv->vector[i] > max)
			max = iv->vector[i];
	return(max);
}

int32_t IntVector_sum(IntVector * iv)
{
	int32_t sum = 0, i;
	for(i = 0; i < iv->curIdx; i++)
		sum += iv->vector[i];
	return sum;
}

double IntVector_mean(IntVector * iv)
{
	int32_t sum = IntVector_sum(iv);
	int32_t len = iv->curIdx;
	double mean = (double)sum/(double)iv->curIdx;
	return mean;
}



