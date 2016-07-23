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
#include <math.h>

#include "definitions.h"
#include "mutlist.h"

void MutList_init(MutList * ml, int32_t sampleSize)
{
	int32_t i,j;
	ml->list = (int32_t **)malloc(sizeof(int32_t *) * (size_t)sampleSize);
	CHECKPOINTER(ml->list);
	ml->curIdxs = (int32_t *)malloc(sizeof(int32_t) * (size_t)sampleSize);
	CHECKPOINTER(ml->curIdxs);
	ml->curSizes = (int32_t *)malloc(sizeof(int32_t) * (size_t)sampleSize);
	CHECKPOINTER(ml->curSizes);
	for(i = 0; i < sampleSize; i++)
	{
		ml->list[i] = (int32_t *)malloc(sizeof(int32_t) * (size_t)100);
		CHECKPOINTER(ml->list[i]);
		for(j = 0; j < 100; j++)
			ml->list[i][j] = -1;		// Mutation ID is never -1.
		ml->curIdxs[i] = 0;
		ml->curSizes[i] = 100;
	}
	ml->nodeIDs = (Node **)malloc(sizeof(Node *) * (size_t)sampleSize);
	CHECKPOINTER(ml->nodeIDs);
	ml->sampleSize = sampleSize;
	ml->S = 0;
	return;
}

void MutList_free(MutList * ml)
{
	int32_t i;
	for(i = 0; i < ml->sampleSize; i++)
		free(ml->list[i]);
	free(ml->curIdxs);
	free(ml->curSizes);
	free(ml->nodeIDs);
	free(ml->list);
	return;
}


// When you don't want the extra malloc/free cycles from _init, but multiple reps are performed.
void MutList_reset(MutList * ml)
{
	int32_t i,j;
	for(i = 0; i < ml->sampleSize; i++)
	{
		ml->list[i] = (int32_t *)realloc(ml->list[i],sizeof(int32_t) * (size_t)100);
		CHECKPOINTER(ml->list[i]);
		for(j = 0; j < 100; j++)
			ml->list[i][j] = -1;
		ml->curIdxs[i] = 0;
		ml->curSizes[i] = 100;
	}
	ml->S = 0;
	return;
}

void MutList_add(MutList * ml, Node * mutNode, int32_t mutIdx)
{
	int32_t i, j, nodeIdx = -1;
	size_t newSize;
	nodeIdx = mutNode->idx;
	if(nodeIdx == -1)
		PERROR("nodeIdx not found in MutList_add()");
	if(ml->curIdxs[nodeIdx] + 1 > ml->curSizes[nodeIdx])
	{
		newSize = (size_t)(ml->curSizes[nodeIdx] + 100);
		ml->list[nodeIdx] = (int32_t *)realloc(ml->list[nodeIdx], sizeof(int32_t) * newSize);
		CHECKPOINTER(ml->list[nodeIdx]);
		for(j = ml->curSizes[nodeIdx]; j < newSize; j++)
			ml->list[nodeIdx][j] = -1;		// Initialize to -1.
		ml->curSizes[nodeIdx] += 100;
	}
	ml->list[nodeIdx][(ml->curIdxs[nodeIdx])++] = mutIdx;
	return;
}

void MutList_print_sequences(MutList * ml, FILE * fout)
{
	int32_t i, j, curPos, curMutIdx, maxMutIdx = 0;
	if(ml->S > 0)		// only print if there are segregating sites
	{
		for(i = 0; i < ml->sampleSize; i++)
		{
			curPos = 0;
			for(j = 0; j < ml->S; j++)
			{
				if(j == ml->list[i][curPos])
				{
					fprintf(fout,"1");
					curPos++;
				}
				else
					fprintf(fout,"0");
			}
			fprintf(fout,"\n");
		}
	}
	fprintf(fout,"\n");
	return;
}

// must be converted using seq2tr from the genetree package.
// populations is a ml->sampleSize-length array of ints with population IDs.
void MutList_print_genetree_format(MutList * ml, int32_t * populations, FILE * fout)
{
	int32_t i, j, curPos, curMutIdx, maxMutIdx = 0;
	if(ml->S > 0)		// only print if there are segregating sites
	{
		for(i = 0; i < ml->sampleSize; i++)
		{
			fprintf(fout, "%i : ", populations[i]);
			curPos = 0;
			for(j = 0; j < ml->S; j++)
			{
				if(j == ml->list[i][curPos])
				{
					fprintf(fout,"1 ");
					curPos++;
				}
				else
					fprintf(fout,"0 ");
			}
			fprintf(fout,"\n");
		}
	}
	return;
}

void MutList_print_ms_format(MutList * ml, FILE * fout)
{
	int32_t i, j, curPos, curMutIdx, maxMutIdx = 0;
	if(ml->S > 0)		// only print if there are segregating sites
	{
		for(i = 0; i < ml->sampleSize; i++)
		{
			curPos = 0;
			for(j = 0; j < ml->S; j++)
			{
				if(j == ml->list[i][curPos])
				{
					fprintf(fout,"1");
					curPos++;
				}
				else
					fprintf(fout,"0");
			}
			fprintf(fout,"\n");
		}
	}
	return;
}

void MutList_print_list(MutList * ml, FILE * fout)
{
	int32_t i, j, curPos, curMutIdx, maxMutIdx = 0;
	for(i = 0; i < ml->sampleSize; i++)
	{
		for(j = 0; j < ml->curIdxs[i]; j++)
			printf("%i ",ml->list[i][j]);
		fprintf(fout,"\n");
	}
	fprintf(fout,"\n");
	return;
}

double MutList_calculate_pi(MutList * ml)
{
	int32_t i, j, k, sampleSize;
	int32_t differences = 0;
	double pi;

	sampleSize = ml->sampleSize;
	// Easy: make a table of counts of mutations. For each indiv in the pair, add a count fo
	// 		 each of their mutations. Then count the number of 1's.
	if(ml->S < 1)
		PERROR("ml->S < 1 in MutList_calculate_pi");
	int32_t * varCounts = (int32_t *)malloc(sizeof(int32_t) * (size_t) ml->S);
	CHECKPOINTER(varCounts);
	for(i = 0; i < sampleSize-1; i++)
	{
		for(j = i+1; j < sampleSize; j++)
		{
			for(k = 0; k < ml->S; k++)
				varCounts[k] = 0;
			for(k = 0; k < ml->curIdxs[i]; k++)
				varCounts[ ml->list[i][k] ]++;
			for(k = 0; k < ml->curIdxs[j]; k++)
				varCounts[ ml->list[j][k] ]++;
			for(k = 0; k < ml->S; k++)
				if(varCounts[k] == 1)
					differences++;
		}
	}
	int32_t sampChoose2 = sampleSize * (sampleSize - 1) / 2;
	pi = (double)differences/sampChoose2;
	free(varCounts);
	return(pi);
}

double MutList_calculate_Tajimas_D(MutList * ml)
{
	int32_t i, n = ml->sampleSize;
	int32_t S = ml->S;
	if(S == 0)
		return 0.0;

	if(ml->S == -1 || ml->S == 0)
		PERROR("Bad ml->S in MutList_calculate_Tajimas_D");
	double pi = MutList_calculate_pi(ml);
	double g1 = 0.0;
	double g2 = 0.0;
	for(i = 1; i < n; i++)
	{
		g1 += 1.0/(double)i;
		g2 += 1.0/(double)(i*i);
	}
	double b1 = (double)(n+1)/(double)(3*(n-1));
	double b2 = (double)(2*(n*n + n + 3))/(double)(9*n*(n-1));
	double c1 = b1 - 1.0/g1;
	double c2 = b2 - (double)(n+2)/(double)(n*g1) + g2/(g1*g1);
	double vhat = (c1*S)/g1 + (c2*S*(S-1))/(g1*g1 + g2);
	double denom = sqrt(vhat);
	double D = (pi - (double)S/g1)/denom;
	return D;
}

void SFS_init(SFS * sfs, int32_t sampleSize)
{
	int32_t i;
	sfs->counts = (int32_t *)malloc(sizeof(int32_t) * sampleSize);
	for(i = 0; i < sampleSize; i++)
		sfs->counts[i] = 0;
	sfs->sampleSize = sampleSize;
	return;
}

void SFS_free(SFS * sfs)
{
	free(sfs->counts);
	sfs->sampleSize = 0;
	return;
}

void SFS_count(SFS * sfs, MutList * ml)
{
	int32_t i,j,k;
	int32_t maxMutationIdx = 0;
	int32_t currentCount;
	for(i = 0; i < ml->sampleSize; i++)
		for(j = 0; j < ml->curSizes[i]; j++)
			if(ml->list[i][j] > maxMutationIdx)
				maxMutationIdx = ml->list[i][j];
	int32_t numMutations = maxMutationIdx + 1;		// For clarity.
	for(i = 0; i < sfs->sampleSize; i++)
		sfs->counts[i] = 0;
	for(i = 0; i < numMutations; i++)		// i is the current mutation.
	{
		currentCount = 0;
		for(j = 0; j <= ml->sampleSize; j++)
		{
			for(k = 0; k < ml->curIdxs[j]; k++)
			{
				if(ml->list[j][k] == i)
					currentCount++;
			}
		}
		sfs->counts[currentCount-1]++;
	}
	return;
}

void SFS_print_X(SFS * sfs, int32_t N, double mu, FILE * fout)
{
	int32_t i;
	int32_t currentCount;
	double zeroedNum;
	double theta = (double)N * 4.0 * mu;
	for(i = 1; i < sfs->sampleSize; i++)
	{
		zeroedNum = (sfs->counts[i-1] - theta/(double)i)/(theta/(double)i);
		fprintf(fout,"%i\t%f\n",i, zeroedNum);
	}
	return;
}
