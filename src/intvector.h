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
