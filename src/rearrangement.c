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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "coalpedigree.h"
#include "genegenealogy.h"
#include "rearrangement.h"

int compare_int32ts_(const void * p1, const void * p2)
{
	return (*(int32_t *)p1 - *(int32_t *)p2);
}

// can compare just the first item because each index is *unique*... there can be no ties.
int compare_rearrgroups_(const void * p1, const void * p2)
{
	return ((*(RearrGroup *)p1).group[0] - ((*(RearrGroup *)p2).group[0]));
}

/* ///////////////////////////////////////////////////////////////////////
 * Simultation-based calculation of rearrangement probabilities
 * /////////////////////////////////////////////////////////////////////// */

// note that numIndivs is the number of *chromosomes*, not diploid individuals
void Rearrangement_init(Rearrangement * rearr, int32_t numIndivs, int32_t numDemes)
{
	int32_t i;
	rearr->maxSize = numIndivs;
	rearr->numDemes = numDemes;
	rearr->groups = (RearrGroup *)malloc(sizeof(RearrGroup) * rearr->maxSize);
	CHECKPOINTER(rearr->groups);
	for(i = 0; i < rearr->maxSize; i++)
	{
		rearr->groups[i].group = (int32_t *)malloc(sizeof(int32_t) * rearr->maxSize);
		CHECKPOINTER(rearr->groups[i].group);
		rearr->groups[i].groupSize = 0;
		rearr->groups[i].deme = -1;
        rearr->groups[i].countFactor = 1.0;
	}
	rearr->numGroups = 0;
	rearr->count = 0;
	rearr->prob = 0.0;
    rearr->countFactor = 1;
	return;
}

void Rearrangement_reset(Rearrangement * rearr)
{
	int32_t i, j;
	for(i = 0; i < rearr->maxSize; i++)
	{
		for(j = 0; j < rearr->groups[i].groupSize; j++)
			rearr->groups[i].group[j] = -1;
		rearr->groups[i].groupSize = 0;
		rearr->groups[i].deme = -1;
	}
	rearr->numGroups = 0;
	rearr->count = 0;
	rearr->prob = 0.0;
	return;
}

void Rearrangement_free(Rearrangement * rearr)
{
	int32_t i;
	for(i = 0; i < rearr->maxSize; i++)
    {
		free(rearr->groups[i].group);
    }
	free(rearr->groups);
	return;
}

int32_t Rearrangement_check_equal(Rearrangement * rearr1, Rearrangement * rearr2)
{
	int32_t i, j;
	assert(rearr1->numDemes == rearr2->numDemes);
	if(rearr1->numGroups != rearr2->numGroups)
		return 0;
	for(i = 0; i < rearr1->numGroups; i++)
	{
		if(rearr1->groups[i].groupSize != rearr2->groups[i].groupSize)
			return 0;
	}
	for(i = 0; i < rearr1->numGroups; i++)
	{
		for(j = 0; j < rearr1->groups[i].groupSize; j++)
		{
			if(rearr1->groups[i].group[j] != rearr2->groups[i].group[j])
				return 0;
		}
	}
	return 1;
}

int32_t Rearrangement_check_equal_twodemes(Rearrangement * rearr1, Rearrangement * rearr2)
{
	int32_t i, j;

	assert(rearr1->numDemes == rearr2->numDemes);
    assert(rearr1->numDemes == 2);
    assert(rearr1->numGroups > 0 && rearr2->numGroups > 0);

	if(rearr1->numGroups != rearr2->numGroups)
		return 0;
	for(i = 0; i < rearr1->numGroups; i++)
	{
		if(rearr1->groups[i].groupSize != rearr2->groups[i].groupSize)
			return 0;
		if(rearr1->groups[i].deme != rearr2->groups[i].deme)
			return 0;
	}
	for(i = 0; i < rearr1->numGroups; i++)
	{
		for(j = 0; j < rearr1->groups[i].groupSize; j++)
		{
			if(rearr1->groups[i].group[j] != rearr2->groups[i].group[j])
				return 0;
		}
	}
	return 1;
}

void Rearrangement_print(Rearrangement * rearr, FILE * output)
{
	int32_t i, j;
	fprintf(output, "|");
	for(i = 0; i < rearr->numGroups; i++)
	{
		for(j = 0; j < rearr->groups[i].groupSize-1; j++)
			fprintf(output, "%i ", rearr->groups[i].group[j]);
        fprintf(output, "%i", rearr->groups[i].group[rearr->groups[i].groupSize-1]);
		fprintf(output, "|");
	}
	//fprintf(output, "\n");
	fprintf(output, "\t");
	return;
}

void Rearrangement_get_descendent_nodes_recursive_(Node * curNode, Rearrangement * rearr, int32_t idx, int32_t * pcurArrayIdx)
{
	if(curNode->desc1 == NULL)
	{
		//printf("NULL\n");
		rearr->groups[idx].group[*pcurArrayIdx] = curNode->idx;
		rearr->groups[idx].groupSize++;
		//REPORTI(curNode->idx);
		//REPORTI(rearr->groups[idx].groupSize);
		(*pcurArrayIdx)++;
		return;
	}
	//printf("NOT NULL\n");
	Rearrangement_get_descendent_nodes_recursive_(curNode->desc1, rearr, idx, pcurArrayIdx);
	Rearrangement_get_descendent_nodes_recursive_(curNode->desc2, rearr, idx, pcurArrayIdx);
	return;
}

void Rearrangement_get_descendent_nodes_(GeneGenealogy * tree, Rearrangement * rearr, int32_t idx)
{
	// use tree->nodeListContainer->nodes[i] to fill out group i of rearr.
	// call a recursive function that adds an index to rearr if the node is a
	// base node (without descendents). keep track of the index (count) of the
	// number of nodes that has been added.
	int32_t curArrayIdx = 0;
	Rearrangement_get_descendent_nodes_recursive_(tree->nodeListContainer->nodeList[idx], rearr, idx, &curArrayIdx);
	return;
}


void Rearrangement_get_rearr_from_genealogy(GeneGenealogy * tree, Rearrangement * rearr) 
{
	int32_t i;
	NodeList * list = tree->nodeListContainer;
	rearr->numGroups = list->numNodes;
	//REPORTI(rearr->numGroups);
	for(i = 0; i < rearr->numGroups; i++)
		Rearrangement_get_descendent_nodes_(tree, rearr, i);
	//Rearrangement_print(rearr, stdout);
	return;
}

void Rearrangement_get_descendent_nodes_diploid_recursive_(Node * curNode, Rearrangement * rearr, int32_t idx, int32_t * pCurArrayIdx)
{
	if(curNode->desc1 == NULL)
	{
		//printf("NULL\n");
		//rearr->groups[idx][*pCurArrayIdx] = curNode->idx / 2;
		rearr->groups[idx].group[*pCurArrayIdx] = curNode->idx;
		rearr->groups[idx].groupSize++;
		//REPORTI(curNode->idx);
		//REPORTI(rearr->groups[idx].groupSize);
		(*pCurArrayIdx)++;
		return;
	}
	//printf("NOT NULL\n");
	Rearrangement_get_descendent_nodes_diploid_recursive_(curNode->desc1, rearr, idx, pCurArrayIdx);
	Rearrangement_get_descendent_nodes_diploid_recursive_(curNode->desc2, rearr, idx, pCurArrayIdx);
	return;
}

void Rearrangement_get_descendent_nodes_diploid_(GeneGenealogy * tree, Rearrangement * rearr, int32_t idx)
{
	// use tree->nodeListContainer->nodes[i] to fill out group i of rearr.
	// call a recursive function that adds an index to rearr if the node is a
	// base node (without descendents). keep track of the index (count) of the
	// number of nodes that has been added.
	int32_t curArrayIdx = 0;
	Rearrangement_get_descendent_nodes_diploid_recursive_(tree->nodeListContainer->nodeList[idx], rearr, idx, &curArrayIdx);
	return;
}

void Rearrangement_get_rearr_from_genealogy_diploid(GeneGenealogy * tree, Rearrangement * rearr) 
{
	int32_t i;
	NodeList * list = tree->nodeListContainer;
	rearr->numGroups = list->numNodes;
	//REPORTI(rearr->numGroups);
	for(i = 0; i < rearr->numGroups; i++)
		Rearrangement_get_descendent_nodes_diploid_(tree, rearr, i);
	//Rearrangement_print(rearr, stdout);
	return;
}

void Rearrangement_sort_rearr(Rearrangement * rearr)
{
	int32_t i;
	for(i = 0; i < rearr->numGroups; i++)
		qsort((void *)rearr->groups[i].group, rearr->groups[i].groupSize, sizeof(int32_t), compare_int32ts_);
    qsort((void *)rearr->groups, rearr->numGroups, sizeof(RearrGroup), compare_rearrgroups_);
	return;
}

void Rearrangement_tally_rearrangement(Rearrangement ** rearrs, int32_t * numRearrs)
{
	int32_t i, equal, found = 0;
	for(i = 0; i < *numRearrs; i++)
	{
		equal = Rearrangement_check_equal(rearrs[i], rearrs[*numRearrs]);
		if(equal)
		{
			(rearrs[i]->count)++;
			// also have to reset the current rearr
			Rearrangement_reset(rearrs[*numRearrs]);
			found = 1;
			break;
		}
	}
	if(!found)
	{
		// allocate another Rearrangement
		rearrs[*numRearrs]->count = 1;
		(*numRearrs)++;
		rearrs[(*numRearrs)] = (Rearrangement *)malloc(sizeof(Rearrangement));
		CHECKPOINTER(rearrs[(*numRearrs)]);
		Rearrangement_init(rearrs[(*numRearrs)], rearrs[0]->maxSize, 1);		// 1 deme
	}
	return;
}

void Rearrangement_simulate_reconfiguration_probs(CoalPedigree * ped, int32_t * indivs, int32_t numIndivs, int32_t numGens, int32_t numReps)
{
	// simulate gene genealogies for some number of generations, stopping when
	// you may have more than one segregating ancestral lineage. the number of
	// lineages at that point is the number of 'groups' in the reconfiguration.
	int32_t i, j, k, done, indivDone, match;
	twoints * lastGen = (twoints *)malloc(sizeof(twoints) * (size_t)ped->popnSize);
	CHECKPOINTER(lastGen);
	// note that if you consider *all* pedigrees, 10000 may not be enough.
	// How many rearrangements can a single pedigree produce? This is not an
	// easy question. For most pedigrees here, 10000 is probably more than
	// sufficient.
	Rearrangement ** rearrs = (Rearrangement **)malloc(sizeof(Rearrangement *) * 10000);
	CHECKPOINTER(rearrs);
	rearrs[0] = (Rearrangement *)malloc(sizeof(Rearrangement));
	CHECKPOINTER(rearrs[0]);
	Rearrangement_init(rearrs[0], numIndivs, 1);		// 1 deme
	//
	// always have at least one Rearrangement: that's the one being filled in
	// according to the simulated genealogy after simulating a genealogy, take
	// rearrrs[numRearrs] and fill it in according to the genealogy then check
	// it against rearrs[0] through rearrs[numRearrs-1]; if it matches, add a
	// count to the matching rearrangement and use the same rearrs[numRearrs]
	// next time. If no match is found, allocate another rearr in rearrs and
	// increment numRearrs, then use that new rearr as the focal one next time
	// around.
	//
	int32_t numRearrs = 0;
	for(i = 0; i < ped->popnSize; i++)
	{
		lastGen[i][0] = ped->relationships[numGens][i][0];
		lastGen[i][1] = ped->relationships[numGens][i][1];
		ped->relationships[numGens][i][0] = 2*i;
		ped->relationships[numGens][i][1] = 2*i + 1;
	}
	CHECKPOINTER(rearrs);
	GeneGenealogy tree;
	for(i = 0; i < numReps; i++)
	{
		// numGens instead of numGens+1 because sim_gene_genealogy coalesces
		// with prob rbern(0.5) if they end up in the same individual. check
		// this.
		//GeneGenealogy_sim_gene_genealogy_head(ped, &tree, indivs, numIndivs, numGens+1);
		GeneGenealogy_sim_gene_genealogy_head(ped, &tree, indivs, numIndivs, numGens);
		// get the rearrangement from the genealogy
		Rearrangement_get_rearr_from_genealogy(&tree, rearrs[numRearrs]);
		Rearrangement_sort_rearr(rearrs[numRearrs]);
		// tally the rearrangement
		Rearrangement_tally_rearrangement(rearrs, &numRearrs);
		GeneGenealogy_free(&tree);
	}
	for(i = 0; i < ped->popnSize; i++)
	{
		ped->relationships[numGens][i][0] = lastGen[i][0];
		ped->relationships[numGens][i][1] = lastGen[i][1];
	}
	for(i = 0; i < numRearrs; i++)
	{
		rearrs[i]->prob = (double)rearrs[i]->count / (double)numReps;
		printf("prob[%i] = %f\n", i, rearrs[i]->prob);
	}
	for(i = 0; i <= numRearrs; i++)
	{
		Rearrangement_free(rearrs[i]);
		free(rearrs[i]);
	}
	free(rearrs);
	free(lastGen);
	return;
}

void RearrContainer_print(RearrContainer * rearrs, FILE * fout)
{
    int32_t i, j, k;
    fprintf(fout, "-------------------\nrearrContainer: %p\n", rearrs);
    fprintf(fout, "\tnumRearrs: %i\n", rearrs->curNumRearrs);
    fprintf(fout, "\ttotalCount: %i\n\n", rearrs->totalCount);
    for(i = 0; i < rearrs->curNumRearrs; i++)
    {
        fprintf(fout, "rearr ");
        Rearrangement_print(rearrs->rearrs[i], fout);
        fprintf(fout, "countFactor: %i\n", rearrs->rearrs[i]->countFactor);
        for(j = 0; j < rearrs->rearrs[i]->numGroups; j++)
        {
            for(k = 0; k < rearrs->rearrs[i]->groups[j].groupSize; k++)
                fprintf(fout, "%i ", rearrs->rearrs[i]->groups[j].group[k]);
            fprintf(fout, "-- countFactor: %i\n", rearrs->rearrs[i]->groups[j].countFactor);
        }
    }
}

void Rearrangement_simulate_reconfiguration_probs_diploid(CoalPedigree * ped, int32_t * indivs, int32_t numIndivs, int32_t numGens, int32_t numReps, FILE * output)
{
	// simulate gene genealogies for some number of generations, stopping when
	// you may have more than one segregating ancestral lineage. the number of
	// lineages at that point is the number of 'groups' in the reconfiguration.
	int32_t i, j, k, done, indivDone, match;
	twoints * lastGen;
	// note that if you consider *all* pedigrees, 10000 may not be enough.
	// How many rearrangements can a single pedigree produce? This is not an
	// easy question. For most pedigrees here, 10000 is probably more than
	// sufficient.
	Rearrangement ** rearrs = (Rearrangement **)malloc(sizeof(Rearrangement *) * 10000);
	CHECKPOINTER(rearrs);
	rearrs[0] = (Rearrangement *)malloc(sizeof(Rearrangement));
	CHECKPOINTER(rearrs[0]);
	Rearrangement_init(rearrs[0], 2*numIndivs, 1);		// 1 deme
	//
	// always have at least one Rearrangement: that's the one being filled in
	// according to the simulated genealogy after simulating a genealogy, take
	// rearrrs[numRearrs] and fill it in according to the genealogy then check
	// it against rearrs[0] through rearrs[numRearrs-1]; if it matches, add a
	// count to the matching rearrangement and use the same rearrs[numRearrs]
	// next time. If no match is found, allocate another rearr in rearrs and
	// increment numRearrs, then use that new rearr as the focal one next time
	// around.
	//
	int32_t numRearrs = 0;
	if(!ped->variablePopnSize)
	{
		lastGen = (twoints *)malloc(sizeof(twoints) * (size_t)ped->popnSize);
		CHECKPOINTER(lastGen);
		for(i = 0; i < ped->popnSize; i++)
		{
			lastGen[i][0] = ped->relationships[numGens-1][i][0];
			lastGen[i][1] = ped->relationships[numGens-1][i][1];
			ped->relationships[numGens-1][i][0] = 2*i;
			ped->relationships[numGens-1][i][1] = 2*i + 1;
		}
	}
	else // ped->variablePopnSize
	{
		lastGen = (twoints *)malloc(sizeof(twoints) * (size_t)ped->popnSizes[numGens-1]);
		CHECKPOINTER(lastGen);
		for(i = 0; i < ped->popnSizes[numGens-1]; i++)
		{
			lastGen[i][0] = ped->relationships[numGens-1][i][0];
			lastGen[i][1] = ped->relationships[numGens-1][i][1];
			ped->relationships[numGens-1][i][0] = 2*i;
			ped->relationships[numGens-1][i][1] = 2*i + 1;
		}
	}

	GeneGenealogy tree;
	for(i = 0; i < numReps; i++)
	{
		// numGens instead of numGens+1 because sim_gene_genealogy coalesces
		// with prob rbern(0.5) if they end up in the same individual. check
		// this.
		//GeneGenealogy_sim_gene_genealogy_head(ped, &tree, indivs, numIndivs, numGens+1);
		GeneGenealogy_sim_gene_genealogy_diploid_head(ped, &tree, indivs, numIndivs, numGens);
		// get the rearrangement from the genealogy
		Rearrangement_get_rearr_from_genealogy_diploid(&tree, rearrs[numRearrs]);
		Rearrangement_sort_rearr(rearrs[numRearrs]);
		// tally the rearrangement
		Rearrangement_tally_rearrangement(rearrs, &numRearrs);
		GeneGenealogy_free(&tree);
	}
	// putting the relationships back
	if(!ped->variablePopnSize)
	{
		for(i = 0; i < ped->popnSize; i++)
		{
			ped->relationships[numGens-1][i][0] = lastGen[i][0];
			ped->relationships[numGens-1][i][1] = lastGen[i][1];
		}
	}
	else // ped->variablePopnSize
	{
		for(i = 0; i < ped->popnSizes[numGens-1]; i++)
		{
			ped->relationships[numGens-1][i][0] = lastGen[i][0];
			ped->relationships[numGens-1][i][1] = lastGen[i][1];
		}
	}

	for(i = 0; i < numRearrs; i++)
	{
		rearrs[i]->prob = (double)rearrs[i]->count / (double)numReps;
		Rearrangement_print(rearrs[i], output);
		fprintf(output, "%f\n", rearrs[i]->prob);
	}
	for(i = 0; i <= numRearrs; i++)
	{
		Rearrangement_free(rearrs[i]);
		free(rearrs[i]);
	}
	free(rearrs);
	free(lastGen);
	return;
}

/* ----------------------------------------------
 * -----------------Two demes--------------------
 * ---------------------------------------------- */

void Rearrangement_print_twodemes(Rearrangement * rearr, FILE * output)
{
	int32_t i, j;
	fprintf(output, "|");
	for(i = 0; i < rearr->numGroups; i++)
	{
		for(j = 0; j < rearr->groups[i].groupSize-1; j++)
			fprintf(output, "%i ", rearr->groups[i].group[j]);
        fprintf(output, "%i", rearr->groups[i].group[rearr->groups[i].groupSize-1]);
		fprintf(output, ":%i|", rearr->groups[i].deme);
	}
	//fprintf(output, "\n");	// still going to print the probability
	fprintf(output, "\t");
	return;
}

void Rearrangement_get_rearr_from_genealogy_diploid_twodemes(GeneGenealogy * tree, Rearrangement * rearr, int32_t demeSize) 
{
	int32_t i;
	NodeList * list = tree->nodeListContainer;
	rearr->numGroups = list->numNodes;
	int32_t deme;
	for(i = 0; i < rearr->numGroups; i++)
	{
		// n.b.! it's 2*demeSize because all relationships of indiv i are made 
		// to be 2*i during Rearrangement calculations!
		deme = (list->nodeList[i]->indiv >= 2*demeSize);
		// !!!!
		rearr->groups[i].deme = deme;
		Rearrangement_get_descendent_nodes_diploid_(tree, rearr, i);
	}
	return;
}

void Rearrangement_tally_rearrangement_twodemes(Rearrangement ** rearrs, int32_t * pnumRearrs)
{
	int32_t i, equal, found = 0;
	for(i = 0; i < *pnumRearrs; i++)
	{
		equal = Rearrangement_check_equal_twodemes(rearrs[i], rearrs[*pnumRearrs]);
		if(equal)
		{
			(rearrs[i]->count)++;
			// also have to reset the current rearr
			Rearrangement_reset(rearrs[*pnumRearrs]);
			found = 1;
			break;
		}
	}
	if(!found)
	{
		// allocate another Rearrangement2D
		rearrs[*pnumRearrs]->count = 1;
		(*pnumRearrs)++;
		rearrs[(*pnumRearrs)] = (Rearrangement *)malloc(sizeof(Rearrangement));
		CHECKPOINTER(rearrs[(*pnumRearrs)]);
		Rearrangement_init(rearrs[(*pnumRearrs)], rearrs[0]->maxSize, 2);		// 2 demes
	}
	return;
}

void Rearrangement_simulate_reconfiguration_probs_diploid_twodemes(CoalPedigree * ped, int32_t * indivs, int32_t numIndivs, int32_t numGens, int32_t numReps, FILE * output)
{
    assert(ped->numDemes == 2);
	// simulate gene genealogies for some number of generations, stopping when
	// you may have more than one segregating ancestral lineage. the number of
	// lineages at that point is the number of 'groups' in the reconfiguration.
	int32_t i, j, k, done, indivDone, match;
	twoints * lastGen = (twoints *)malloc(sizeof(twoints) * (size_t)ped->popnSize);
	CHECKPOINTER(lastGen);
	// note that if you consider *all* pedigrees, 10000 may not be enough.
	// How many rearrangements can a single pedigree produce? This is not an
	// easy question. For most pedigrees here, 10000 is probably more than
	// sufficient.
	Rearrangement ** rearrs = (Rearrangement **)malloc(sizeof(Rearrangement *) * 10000);
	CHECKPOINTER(rearrs);
	rearrs[0] = (Rearrangement *)malloc(sizeof(Rearrangement));
	CHECKPOINTER(rearrs[0]);
	Rearrangement_init(rearrs[0], 2*numIndivs, 2);		// 2 demes
	//
	// always have at least one Rearrangement: that's the one being filled in
	// according to the simulated genealogy after simulating a genealogy, take
	// rearrrs[numRearrs] and fill it in according to the genealogy then check
	// it against rearrs[0] through rearrs[numRearrs-1]; if it matches, add a
	// count to the matching rearrangement and use the same rearrs[numRearrs]
	// next time. If no match is found, allocate another rearr in rearrs and
	// increment numRearrs, then use that new rearr as the focal one next time
	// around.
	//
	int32_t numRearrs = 0;
	for(i = 0; i < ped->popnSize; i++)
	{
		lastGen[i][0] = ped->relationships[numGens-1][i][0];
		lastGen[i][1] = ped->relationships[numGens-1][i][1];
		ped->relationships[numGens-1][i][0] = 2*i;
		ped->relationships[numGens-1][i][1] = 2*i + 1;
	}
	GeneGenealogy tree;
	for(i = 0; i < numReps; i++)
	{
		// numGens instead of numGens+1 because sim_gene_genealogy coalesces
		// with prob rbern(0.5) if they end up in the same individual. check
		// this.
		//GeneGenealogy_sim_gene_genealogy_head(ped, &tree, indivs, numIndivs, numGens+1);
		GeneGenealogy_sim_gene_genealogy_diploid_head(ped, &tree, indivs, numIndivs, numGens);
		// get the rearrangement from the genealogy
		Rearrangement_get_rearr_from_genealogy_diploid_twodemes(&tree, rearrs[numRearrs], ped->popnSize/2);	// n.b. numDemes==2
		//Rearrangement2D_print(rearrs[numRearrs], stdout);
		//printf("\n");
		Rearrangement_sort_rearr(rearrs[numRearrs]);
		//Rearrangement2D_print(rearrs[numRearrs], stdout);
		//printf("\n\n");
		// tally the rearrangement
		Rearrangement_tally_rearrangement_twodemes(rearrs, &numRearrs);
		GeneGenealogy_free(&tree);
	}
	for(i = 0; i < ped->popnSize; i++)
	{
		ped->relationships[numGens-1][i][0] = lastGen[i][0];
		ped->relationships[numGens-1][i][1] = lastGen[i][1];
	}
	for(i = 0; i < numRearrs; i++)
	{
		rearrs[i]->prob = (double)rearrs[i]->count / (double)numReps;
		Rearrangement_print_twodemes(rearrs[i], output);
		fprintf(output, "%f\n", rearrs[i]->prob);
	}
	for(i = 0; i <= numRearrs; i++)
	{
		Rearrangement_free(rearrs[i]);
		free(rearrs[i]);
	}
	free(rearrs);
	free(lastGen);
	return;
}

int32_t Rearrangement_increment_segregation_(int32_t * seg, int32_t len)
{
    int32_t i, j;
    for(i = len-1; i >= 0; i--)
    {
        if(seg[i] == 0)
        {
            seg[i] = 1;
            for(j = i+1; j < len; j++)
                seg[j] = 0;
            return 0;
        }
    }
    return 1;
}


//////////////////////////////////
///////// RearrContainer /////////
//////////////////////////////////

void RearrContainer_init(RearrContainer * rearrs, int32_t maxNumRearrs, int32_t maxGroupSize, int32_t numDemes)
{
    rearrs->curNumRearrs = 0;
    rearrs->maxNumRearrs = maxNumRearrs;
    rearrs->rearrs = (Rearrangement **)malloc(sizeof(Rearrangement *) * maxNumRearrs);
    CHECKPOINTER(rearrs->rearrs);
    rearrs->maxSize = maxGroupSize;
    rearrs->numDemes = numDemes;
    rearrs->totalCount = 0;
    return;
}

void RearrContainer_free(RearrContainer * rearrs)
{
    int32_t i;
    for(i = 0; i < rearrs->curNumRearrs; i++)
    {
        Rearrangement_free(rearrs->rearrs[i]);
        free(rearrs->rearrs[i]);
    }
    free(rearrs->rearrs);
    return;
}

void Rearrangement_copy_rearrangement(Rearrangement * src, Rearrangement * dest)
{
    int32_t i, j;
    dest->numGroups = src->numGroups;
    dest->maxSize = src->maxSize;
    dest->count = src->count;
    for(i = 0; i < src->numGroups; i++)
    {
        dest->groups[i].groupSize = src->groups[i].groupSize;
        for(j = 0; j < src->groups[i].groupSize; j++)
            dest->groups[i].group[j] = src->groups[i].group[j];
        dest->groups[i].deme = src->groups[i].deme;
        dest->groups[i].countFactor = src->groups[i].countFactor;
    }
    return;
}

int32_t Rearrangement_get_count_factor(Rearrangement * rearr)
{
    int32_t countFactor = 1, j;
    for(j = 0; j < rearr->numGroups; j++)
        countFactor *= rearr->groups[j].countFactor;
    return countFactor;
}

void RearrContainer_add_rearr(Rearrangement * rearr, RearrContainer * rearrContainer)
{
    rearrContainer->rearrs[rearrContainer->curNumRearrs] = (Rearrangement *)malloc(sizeof(Rearrangement));
    CHECKPOINTER(rearrContainer->rearrs[rearrContainer->curNumRearrs]);

    Rearrangement * addedRearr = rearrContainer->rearrs[rearrContainer->curNumRearrs];

    Rearrangement_init(addedRearr, rearr->maxSize, rearrContainer->numDemes);
    Rearrangement_copy_rearrangement(rearr, addedRearr);

    addedRearr->numDemes = rearr->numDemes;

    addedRearr->countFactor = rearr->countFactor;

    (rearrContainer->curNumRearrs)++;
    return;
}

void RearrContainer_tally_rearrangement(Rearrangement * rearr, RearrContainer * rearrContainer)
{
    assert(rearr->numDemes == rearrContainer->numDemes);
    assert(rearr->numDemes == 1 || rearr->numDemes == 2);
    rearr->countFactor = Rearrangement_get_count_factor(rearr);
	int32_t i, j, equal = 0, found = 0;
    Rearrangement_sort_rearr(rearr);
	for(i = 0; i < rearrContainer->curNumRearrs; i++)
	{
        if(rearr->numDemes == 1)
            equal = Rearrangement_check_equal(rearr, rearrContainer->rearrs[i]);
        else
        {
            equal = Rearrangement_check_equal_twodemes(rearr, rearrContainer->rearrs[i]);
        }
		if(equal)
		{
            // calculate count factor
            rearr->countFactor = Rearrangement_get_count_factor(rearr);
			rearrContainer->rearrs[i]->countFactor += rearr->countFactor;
			found = 1;
			break;
		}
	}
	if(!found)
        RearrContainer_add_rearr(rearr, rearrContainer);
    rearrContainer->totalCount += rearr->countFactor;
	return;
}

void RearrGroup_copy_group(RearrGroup * src, RearrGroup * dest)
{
    dest->groupSize = src->groupSize;
    dest->deme = src->deme;
    dest->group = src->group;
    dest->countFactor = src->countFactor;
}

void RearrGroup_coalesce_groups(RearrGroup * src, RearrGroup * dest, int32_t curGen, int32_t targetGen)
{
    int32_t i;
    dest->deme = src->deme;
    dest->group = realloc(dest->group, sizeof(int32_t) * (dest->groupSize + src->groupSize));
    CHECKPOINTER(dest->group);
    for(i = 0; i < src->groupSize; i++)
        dest->group[i + dest->groupSize] = src->group[i];
    dest->groupSize += src->groupSize;
    dest->countFactor = dest->countFactor * src->countFactor * (int32_t)pow(2.0, (double)(targetGen-1-curGen));
    // free the src rearrgroup
    free(src->group);
    return;
}

void Rearrangement_copy_lineages(RearrGroup * lineages, RearrGroup * lineagesCons, int32_t numLineages)
{
    // copy entire array
    memcpy(lineagesCons, lineages, sizeof(RearrGroup) * numLineages);
    // now copy individual groups
    int32_t i;
    for(i = 0;  i < numLineages; i++)
    {
        lineagesCons[i].group = (int32_t *)malloc(sizeof(int32_t) * lineagesCons[i].groupSize);
        CHECKPOINTER(lineagesCons[i].group);
        memcpy(lineagesCons[i].group, lineages[i].group, sizeof(int32_t) * lineages[i].groupSize);
        lineagesCons[i].groupSize = lineages[i].groupSize;
    }
    return;
}

void Rearrangement_enumerate_rearrs_recursive_(RearrGroup * lineages, int32_t * lineageIndivs, int32_t numLineages, int32_t curGen, int32_t targetGen, RearrContainer * rearrs, CoalPedigree * ped)
{
    int32_t i, j, k;

    if(curGen == targetGen)
    {
        int32_t maxGroupSizeOrLineageCount = numLineages;
        for(i = 0; i < numLineages; i++)
        {
            if(lineages[i].groupSize > maxGroupSizeOrLineageCount)
                maxGroupSizeOrLineageCount = lineages[i].groupSize;
        }
        // make rearrangement, add it to rearrs
        Rearrangement rearr;
        Rearrangement_init(&rearr, maxGroupSizeOrLineageCount, ped->numDemes);
        for(i = 0; i < numLineages; i++)
        {

            rearr.groups[i].groupSize = lineages[i].groupSize;
            for(j = 0; j < lineages[i].groupSize; j++)
                rearr.groups[i].group[j] = lineages[i].group[j];

            if(ped->numDemes == 1)
                rearr.groups[i].deme = 0;
            else
                rearr.groups[i].deme = (lineageIndivs[i]/2) / (ped->popnSize/ped->numDemes);
            // ^ divide lineageIndivs[i] by 2 because of branching of lastGen

            rearr.groups[i].countFactor = lineages[i].countFactor;
        }
        rearr.numGroups = numLineages;
        RearrContainer_tally_rearrangement(&rearr, rearrs);
        Rearrangement_free(&rearr);
        return;
    }

    int32_t * segregationPattern = (int32_t *)calloc((size_t)numLineages, sizeof(int32_t));
    CHECKPOINTER(segregationPattern);
    for(i = 0; i < numLineages; i++)
        segregationPattern[i] = 0;
    int32_t * parents = (int32_t *)malloc(sizeof(int32_t) * numLineages);
    CHECKPOINTER(parents);

    int32_t numLineagesCons; // number of consolidated lineages

    int32_t * parentsCons = (int32_t *)malloc(sizeof(int32_t) * numLineages);
    CHECKPOINTER(parentsCons);
    RearrGroup * lineagesCons = (RearrGroup *)malloc(sizeof(RearrGroup) * numLineages);
    CHECKPOINTER(lineagesCons);


    int32_t * lineageIndivsCons = (int32_t *)malloc(sizeof(int32_t) * numLineages);
    CHECKPOINTER(lineageIndivsCons);

    int32_t done = 0;
    while(!done)
    {
        Rearrangement_copy_lineages(lineages, lineagesCons, numLineages);
        for(i = 0; i < numLineages; i++)
            parents[i] = ped->relationships[curGen][lineageIndivs[i]][segregationPattern[i]];

        memcpy(lineageIndivsCons, lineageIndivs, sizeof(int32_t) * numLineages);

        memcpy(parentsCons, parents, sizeof(int32_t) * numLineages);

        numLineagesCons = numLineages;

        for(i = 0; i < numLineagesCons-1; i++)
        {
            for(j = i+1; j < numLineagesCons; j++)
            {
                if(lineageIndivsCons[i] == lineageIndivsCons[j] && parentsCons[i] == parentsCons[j])
                {
                    // coalesce
                    RearrGroup_coalesce_groups(&(lineagesCons[j]), &(lineagesCons[i]), curGen, targetGen); // order matters here

                    for(k = j; k < numLineagesCons-1; k++)
                    {
                        RearrGroup_copy_group(&(lineagesCons[k+1]), &(lineagesCons[k]));
                        lineageIndivsCons[k] = lineageIndivsCons[k+1];
                        parentsCons[k] = parentsCons[k+1];
                    }
                    numLineagesCons--;
                    j--;
                }
            }
        }

        Rearrangement_enumerate_rearrs_recursive_(lineagesCons, parentsCons, numLineagesCons, curGen+1, targetGen, rearrs, ped);

        done = Rearrangement_increment_segregation_(segregationPattern, numLineages);

        for(i = 0; i < numLineagesCons; i++)
            free(lineagesCons[i].group);
    }

    free(segregationPattern);
    free(parents);
    free(parentsCons);
    free(lineagesCons);
    free(lineageIndivsCons);

    return;
}

void Rearrangement_calculate_reconfiguration_probs_diploid(CoalPedigree * ped, int32_t * indivs, int32_t numIndivs, int32_t numGens, FILE * output)
{
	int32_t i, j, k;
	twoints * lastGen;
    if(ped->variablePopnSize)
        PERROR("Rearrangement_calculate_reconfiguration_probs_diploid() can only be used with constant-size pedigree");

    lastGen = (twoints *)malloc(sizeof(twoints) * (size_t)ped->popnSize);
    CHECKPOINTER(lastGen);
    for(i = 0; i < ped->popnSize; i++)
    {
        lastGen[i][0] = ped->relationships[numGens-1][i][0];
        lastGen[i][1] = ped->relationships[numGens-1][i][1];
        ped->relationships[numGens-1][i][0] = 2*i;
        ped->relationships[numGens-1][i][1] = 2*i + 1;
    }

    RearrGroup * lineages = (RearrGroup *)malloc(sizeof(RearrGroup) * 2*numIndivs);
	CHECKPOINTER(lineages);
    for(i = 0; i < 2*numIndivs; i++)
    {
        lineages[i].groupSize = 1;
        lineages[i].group = (int32_t *)malloc(sizeof(int32_t));
        lineages[i].group[0] = i;
        lineages[i].deme = 0; // will calculate later if relevant
        lineages[i].countFactor = 1;
    }

    int32_t * lineageIndivs = (int32_t *)malloc(sizeof(int32_t) * 2*numIndivs);
    CHECKPOINTER(lineageIndivs);
    for(i = 0; i < numIndivs; i++)
    {
        lineageIndivs[2*i] = ped->relationships[0][indivs[i]][0];
        lineageIndivs[2*i+1] = ped->relationships[0][indivs[i]][1];
    }

    RearrContainer rearrs;
    RearrContainer_init(&rearrs, 1000, 2*numIndivs, ped->numDemes);

    Rearrangement_enumerate_rearrs_recursive_(lineages, lineageIndivs, 2*numIndivs, 1, numGens, &rearrs, ped);

	// putting the relationships back
    for(i = 0; i < ped->popnSize; i++)
    {
        ped->relationships[numGens-1][i][0] = lastGen[i][0];
        ped->relationships[numGens-1][i][1] = lastGen[i][1];
    }

    double prob;
    for(i = 0; i < rearrs.curNumRearrs; i++)
    {
        if(ped->numDemes == 1)
            Rearrangement_print(rearrs.rearrs[i], output);
        else
            Rearrangement_print_twodemes(rearrs.rearrs[i], output);
        prob = (double)(rearrs.rearrs[i]->countFactor) / (double)(rearrs.totalCount);
		fprintf(output, "%f\n", prob);
    }

    //RearrContainer_print(&rearrs, output);

    RearrContainer_free(&rearrs);
    for(i = 0; i < 2*numIndivs; i++)
        free(lineages[i].group);
    free(lineages);
    free(lineageIndivs);
	free(lastGen);
	return;
}


/*
int32_t main()
{
	CoalPedigree ped;
	int32_t N = 20;
	CoalPedigree_init_singledeme(&ped, 5, 10, 10);

	int32_t indivs[2] = {1,0};
	int32_t numGens = 3;
	int32_t numIndivs = 2;
	int32_t numReps = 100000;

	randseed();
    
	ped.relationships[0][0][0] = 3;
	ped.relationships[0][0][1] = 5;
	ped.relationships[0][1][0] = 2;
	ped.relationships[0][1][1] = 5;

	ped.relationships[1][2][0] = 3;
	ped.relationships[1][2][1] = 9;
	ped.relationships[1][3][0] = 1;
	ped.relationships[1][3][1] = 5;
	ped.relationships[1][5][0] = 1;
	ped.relationships[1][5][1] = 6;

	Rearrangement_simulate_reconfiguration_probs_diploid(&ped, indivs, numIndivs, numGens, numReps, stdout);
    printf("\n");
	Rearrangement_calculate_reconfiguration_probs_diploid(&ped, indivs, numIndivs, numGens, stdout);

    //  simulation
    // |0|1 3|2|       0.063100
    // |0|1 2|3|       0.436610
    // |0|1 2 3|       0.063490
    // |0|1|2 3|       0.062010
    // |0|1|2|3|       0.374790

    //  enumeration
    // |0|1 2|3|       0.437500
    // |0|1|2 3|       0.062500
    // |0|1|2|3|       0.375000
    // |0|1 3|2|       0.062500
    // |0|1 2 3|       0.062500

	CoalPedigree_free(&ped);
	return 0;
}
*/

