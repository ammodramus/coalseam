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

#ifndef COALPEDIGREE_HEADER
#define COALPEDIGREE_HEADER

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include "definitions.h"
#include "options.h"

/*
 * CoalPedigree
 *
 * relationships[t][i][j]
 * 	t: generation under consideration (offspring's generation, back in time)
 * 	i: individual number (0 is female, 1 is male)
 * 	j: parent gender (0 is mother, 1 is father)
 */
/* A function for initializing CoalPedigree structs. Given numGens, numMales,
 * and numFemales, CoalPedigree_init allocates memory for an array of
 * dimensions t,i,j,k, where t is the number of generations, i (offspring
 * gender) = 2, j (parent gender) = 2, and k is the individual (either male
 * or female) number.
 *
 * Note I have removed the sweep functionality from this module.
*/

typedef struct coalpedigree
{
	int32_t numGenerations;
	int32_t numMales;
	int32_t numFemales;
	int32_t popnSize;

	int32_t numDemes;

	int32_t variablePopnSize; 		// == 0 if N is constant, 1 if N changes.
	int32_t * popnSizes;			// only malloc'ed if(variablePopnSize);

	twoints ** relationships;		// container for all the relationships in the pedigree
} CoalPedigree;

void CoalPedigree_init(CoalPedigree * ped, Options * opt);
void CoalPedigree_init_singledeme(CoalPedigree * ped, int32_t numGens, int32_t numFemales, int32_t numMales);
void CoalPedigree_init_ddemes(CoalPedigree * ped, int32_t numGens, int32_t numDemes, int32_t demeSize);
//void CoalPedigree_init(CoalPedigree * ped, int32_t numGens, int32_t numDemes, int32_t demeSize);
void CoalPedigree_init_varN(CoalPedigree * ped, int32_t numGens, int32_t * popnSizes);
void CoalPedigree_free(CoalPedigree * ped);
void CoalPedigree_shuffle(CoalPedigree * ped, Options * opt);
void CoalPedigree_shuffle_dioecious(CoalPedigree * ped);
void CoalPedigree_shuffle_dioecious_fixed_recent_relative(CoalPedigree * ped, int32_t * indivs, int32_t relAncGen);
void CoalPedigree_shuffle_dioecious_fixed_ibd(CoalPedigree * ped, int32_t ibdGen, int32_t clearnGen, int32_t numIndivs, int32_t ibdIndiv, int32_t * nonIBDindivs);
void CoalPedigree_shuffle_dioecious_shared_anc(CoalPedigree * ped, int32_t sharedAncGen, int32_t clearGen, int32_t * indivs, int32_t numIndivs);
void CoalPedigree_shuffle_twopop_nomigration(CoalPedigree * ped);
void CoalPedigree_shuffle_ddeme_nomigration(CoalPedigree * ped);
void CoalPedigree_shuffle_twopop_nomigration_interval(CoalPedigree * ped, int32_t start, int32_t end);
void CoalPedigree_shuffle_ddeme_diploid_stochastic(CoalPedigree * ped, double avgNumMigrants);
void CoalPedigree_shuffle_twopop_diploid_fully_stochastic(CoalPedigree * ped, double avgNumMigrants);
void CoalPedigree_shuffle_twopop_diploid_fully_stochastic_varM(CoalPedigree * ped, double * avgNumMigrants);
void CoalPedigree_shuffle_twopop_diploid_fixed(CoalPedigree * ped, int32_t numMigrants);
void CoalPedigree_shuffle_twopop_diploid_stochastic(CoalPedigree * ped, double avgNumMigrants);
void CoalPedigree_shuffle_twopop_diploid_stochastic_rm(CoalPedigree * ped, double avgNumMigrants);
void CoalPedigree_shuffle_twopop_diploid_stochastic_janet(CoalPedigree * ped, double avgNumMigrants);
void CoalPedigree_shuffle_twopop_diploid_interval(CoalPedigree * ped, int32_t interval, int32_t firstMigGen);
void CoalPedigree_shuffle_twopop_fixed_migrant(CoalPedigree * ped, double avgNumMigrants, int32_t * indivs, int32_t numIndivs, int32_t ancMigGen, int32_t clearGen);
void CoalPedigree_recursive_remove_ancestry_(CoalPedigree * ped, int32_t indiv, int32_t curGen, int32_t ancGen, int32_t anc);
void CoalPedigree_remove_ancestry(CoalPedigree * ped, int32_t indiv, int32_t ancGen, int32_t anc);
void CoalPedigree_shuffle_fixed_recent_relative(CoalPedigree * ped, int32_t * twoindivs, int32_t relAncGen);
int32_t CoalPedigree_find_migrant(CoalPedigree * ped, int32_t migAncGen, int32_t demeIdx);
double CoalPedigree_get_anc_weight(CoalPedigree * ped, int32_t indiv, int32_t target, int32_t ancGen, double * wc, double * wcp);
double CoalPedigree_get_mig_prob(CoalPedigree * ped, int32_t indiv, int32_t migrantGen, double * wc, double * wcp);
double CoalPedigree_shuffle_twopop_fixed_weight_single_migrant(CoalPedigree * ped, int32_t * indivs, int32_t numIndivsPerDeme, int32_t migAncGen, int32_t interval);
void CoalPedigree_shuffle_twopop_fixed_recent_relative(CoalPedigree * ped, int32_t * indivs, int32_t numIndivsPerDeme, int32_t relAncGen, int32_t interval);
void CoalPedigree_shuffle_twopop_haploid_fixed(CoalPedigree * ped, int32_t numMigrants);
void CoalPedigree_shuffle_twopop_haploid_interval(CoalPedigree * ped, int32_t interval);
void CoalPedigree_shuffle_twopop_haploid_stochastic(CoalPedigree * ped, double avgNumMigrants);
void CoalPedigree_shuffle_dioecious_varN(CoalPedigree * ped, int32_t loopGen);
void CoalPedigree_shuffle_monoecious(CoalPedigree * ped);
void CoalPedigree_shuffle_cyc_WF_dio(CoalPedigree * ped);
void CoalPedigree_shuffle_cyc_WF_monoecious(CoalPedigree * ped);
void CoalPedigree_update_type_probs_(double * tprobs, double * cTprobs, int32_t * numGenotypes, int32_t * gPerm, double s, double h);
void CoalPedigree_shuffle_sweep_dio(CoalPedigree * ped, int32_t ts, double s, double h, int32_t requireSweep);
void CoalPedigree_shuffle_sweep_dio_varN(CoalPedigree * ped, int32_t ts, double s, double h, int32_t requireSweep);
int32_t CoalPedigree_sim_coal_time(CoalPedigree * ped, int32_t * indiv);
int32_t CoalPedigree_sim_coal_time_loop(CoalPedigree * ped, int32_t * indiv, int32_t loopTime);
void CoalPedigree_calculate_num_migrants(CoalPedigree * ped, int32_t popn, int32_t * countBucket, int32_t * numMigrants, int32_t * numOffspring);
double CoalPedigree_calculate_migration_prob(CoalPedigree * ped, int32_t numGenerations, int32_t popn, int32_t * indiv, double * weights, double * weightsP);
double CoalPedigree_calculate_migration_prob_singleindiv(CoalPedigree * ped, int32_t numGenerations, int32_t indiv, double * weights, double * weightsP);
void CoalPedigree_calculate_migration_prob_cumulative(CoalPedigree * ped, int32_t numGenerations, int32_t popn, int32_t * indiv, double * weights, double * weightsP, double * cumProb);
void CoalPedigree_print_migrant_weights(CoalPedigree * ped, int32_t numGens, int32_t popn, double ** wc, double ** wcp, FILE * fout);
void CoalPedigree_print_migrant_weights_ddeme(CoalPedigree * ped, int32_t * indiv, int32_t numIndivs, int32_t numGens, double ** wc, double ** wcp, FILE * fout);
void CoalPedigree_calc_coal_probs(CoalPedigree * ped, int32_t * indivs, int32_t maxGen, double * weights, double * weightsP, double * probs);
void CoalPedigree_print_coal_probs(CoalPedigree * ped, int32_t * indivs, int32_t maxGen, double * weights, double * weightsP, FILE * fout);
void CoalPedigree_export_pedigree(CoalPedigree * ped, char * filename);
void CoalPedigree_import_pedigree(CoalPedigree * ped, char * filename);

#endif
