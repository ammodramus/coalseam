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
