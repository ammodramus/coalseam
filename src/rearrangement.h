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

#ifndef REARRANGEMENT_HEADER
#define REARRANGEMENT_HEADER

#include "coalpedigree.h"

typedef struct rearrgroup_
{
	int32_t groupSize;
	int32_t * group;
	int32_t deme;
    int32_t countFactor;
} RearrGroup;

typedef struct rearrangement_
{
	int32_t numGroups;
	int32_t maxSize;	// n.b. = # chromosomes = maximum number of groups = maximum individual group size
	int32_t numDemes;
	RearrGroup * groups;
	double prob;
	int32_t count;
    int32_t countFactor;
} Rearrangement;

typedef struct rearrcontainer_
{
    Rearrangement ** rearrs;
    int32_t curNumRearrs;
    int32_t maxNumRearrs;
    int32_t maxSize; // for each of the individual Rearrangement objects
    int32_t numDemes; // number of demes for every Rearrangement object
    int32_t totalCount;
} RearrContainer;

void Rearrangement_simulate_reconfiguration_probs(CoalPedigree * ped, int32_t * indivs, int32_t numIndivs, int32_t numGens, int32_t numReps);
void RearrContainer_print(RearrContainer * rearrs, FILE * fout);
void Rearrangement_simulate_reconfiguration_probs_diploid(CoalPedigree * ped, int32_t * indivs, int32_t numIndivs, int32_t numGens, int32_t numReps, FILE * output);
void Rearrangement_simulate_reconfiguration_probs_diploid_twodemes(CoalPedigree * ped, int32_t * indivs, int32_t numIndivs, int32_t numGens, int32_t numReps, FILE * output);
void Rearrangement_calculate_reconfiguration_probs_diploid(CoalPedigree * ped, int32_t * indivs, int32_t numIndivs, int32_t numGens, FILE * output);

#endif
