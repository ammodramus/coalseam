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

#ifndef MUTLIST_HEADER
#define MUTLIST_HEADER

#include <stdio.h>
#include "node.h"

typedef struct mutlist
{
	int32_t ** list;
	int32_t * curIdxs;
	int32_t * curSizes;
	Node ** nodeIDs;
	int32_t sampleSize;
	int32_t S;		// Number of segregating sites / mutations.
} MutList;

typedef struct sfs
{
	int32_t * counts;
	int32_t sampleSize;
} SFS;

void MutList_init(MutList * ml, int32_t sampleSize);
void MutList_reset(MutList * ml);
void MutList_add(MutList * ml, Node * mutNode, int32_t mutIdx);
void MutList_print_sequences(MutList * ml, FILE * fout);
void MutList_print_genetree_format(MutList * ml, int32_t * populations, FILE * fout);
void MutList_print_ms_format(MutList * ml, FILE * fout);
void MutList_print_list(MutList * ml, FILE * fout);
double MutList_calculate_pi(MutList * ml);
double MutList_calculate_Tajimas_D(MutList * ml);
void MutList_free(MutList * ml);
void SFS_init(SFS * sfs, int32_t sampleSize);
void SFS_free(SFS * sfs);
void SFS_count(SFS * sfs, MutList * ml);
void SFS_print_X(SFS * sfs, int32_t N, double mu, FILE * fout);

#endif
