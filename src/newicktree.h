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

#ifndef NEWICKTREE_HEADER
#define NEWICKTREE_HEADER

#include <stdint.h>

typedef struct newicktree
{
	char * line;
	int32_t curIdx;
	int32_t curSize;
} NewickTree;

void NewickTree_init(NewickTree * nt);
void NewickTree_add(NewickTree * nt, char * toadd, int32_t toaddLength);
void NewickTree_reset(NewickTree * nt);
void NewickTree_free(NewickTree * nt);

#endif
