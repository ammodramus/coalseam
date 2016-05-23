#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <getopt.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "coalpedigree.h"
#include "definitions.h"
#include "rearrangement.h"
#include "genegenealogy.h"
#include "newicktree.h"
#include "random.h"

void get_init_indivs_pedsim(CoalPedigree * ped, Options * opt, int32_t * indivs)
{
    // shuffle routines that have particular individuals involved:
    // - fixed ibd
    // - fixed migrant

    if(opt->withinSet)
    {
        indivs[0] = runifd(0, ped->popnSize/ped->numDemes-1);
        do{
            indivs[1] = runifd(0, ped->popnSize/ped->numDemes-1);
        } while(indivs[0] == indivs[1]);
        return;
    }
    else // !opt->withinSet (between)
    {
        // first indiv is from first deme, second is from second, regardless of numDemes
        indivs[0] = runifd(0, ped->popnSize/ped->numDemes-1);
        indivs[1] = runifd(ped->popnSize/ped->numDemes, 2*ped->popnSize/ped->numDemes-1);
        return;
    }
}

/* 
 * print_coal_probs     a function to print the exact pairwise coalescence time
 *                      distribution for chromosomes sampled from two separate
 *                      individuals
 * 
 * ped                  initialized pedigree
 * indivs               int32_t vector of length two, containing the two sampled individuals
 * maxGen               stop after maxGen generations
 * fout                 file object for output
 */

void print_coal_probs(CoalPedigree * ped, int32_t * indivs, int32_t maxGen, FILE * fout)
{
	int32_t t, i, j, k, indivIdx, parentCombIdxs[4], popnSize = ped->popnSize, checkInterval, anyWeight;
	double prob;

    checkInterval = maxGen/30;

    double * weights = (double *)malloc(sizeof(double) * ped->popnSize * ped->popnSize);
    CHECKPOINTER(weights);

    double * weightsP = (double *)malloc(sizeof(double) * ped->popnSize * ped->popnSize);
    CHECKPOINTER(weightsP);

	for(i = 0; i < popnSize*popnSize; i++)
		weights[i] = weightsP[i] = 0.0;
	weights[popnSize*indivs[0] + indivs[1]] = 1.0;		// initial configuration.
	for(t = 0; t < maxGen-1; t++)		// note < not <=. Arbitrary.
	{
        prob = 0.0;
		// pass back the weight.
		for(i = 0; i < popnSize; i++)
		{
			for(j = 0; j < popnSize; j++)
			{
				indivIdx = i*popnSize + j;
				if(weights[indivIdx] > 0.0)
				{
					parentCombIdxs[0] = ped->relationships[t][i][0]*popnSize + ped->relationships[t][j][0];
					parentCombIdxs[1] = ped->relationships[t][i][0]*popnSize + ped->relationships[t][j][1];
					parentCombIdxs[2] = ped->relationships[t][i][1]*popnSize + ped->relationships[t][j][0];
					parentCombIdxs[3] = ped->relationships[t][i][1]*popnSize + ped->relationships[t][j][1];
					for(k = 0; k < 4; k++)
						weightsP[parentCombIdxs[k]] += 0.25*weights[indivIdx];
				}
			}
		}
		// trade weights and weightsP, set weightsP to zero
		for(i = 0; i < popnSize * popnSize; i++)
		{
			weights[i] = weightsP[i];
			weightsP[i] = 0.0;
		}
		// check for coalescence
		for(i = 0; i < popnSize; i++)
		{
			indivIdx = i*popnSize + i;
			if(weights[indivIdx] > 0.0)
			{
				prob += 0.5*weights[indivIdx];
				// give all of the non-coalescing weight to one vs. the other (equivalent) configuration
				weightsP[ped->relationships[t+1][i][0]*popnSize + ped->relationships[t+1][i][1]] += 0.5*weights[indivIdx];
				weights[indivIdx] = 0.0;
			}
		}
		fprintf(fout, "%i\t%g\n", t+1, prob);
		if(t % checkInterval == 0)
		{
			anyWeight = 0;
			for(i = 0; i < popnSize*popnSize; i++)
			{
				if(weights[i] > 0.0)
				{
					anyWeight = 1;
					break;
				}
			}
			if(!anyWeight)		// if there's no weight left (due to underflow), return
            {
                free(weights);
                free(weightsP);
				return;
            }
		}
	}
    free(weights);
    free(weightsP);
	return;
}

void print_sim_coal_times(CoalPedigree * ped, int32_t * indivs, int32_t numReps, FILE * fout)
{
    int32_t rep;
    for(rep = 0; rep < numReps; rep++)
        fprintf(fout, "%i\n", CoalPedigree_sim_coal_time(ped, indivs));
    return;
}

double simulate_fst(CoalPedigree * ped, int32_t * indivs, int32_t numCoalTimes)
{
    int32_t initIndivs[2], coalSum, i;
    double ETw0, ETw1, ETb, fst;

    // first calculate ET0_0
    coalSum = 0;
    for(i = 0; i < numCoalTimes; i++)
    {
        initIndivs[0] = initIndivs[1] = indivs[0];
        coalSum += CoalPedigree_sim_coal_time(ped, initIndivs);
    }
    ETw0 = (double)coalSum/(double)numCoalTimes;
    
    // then calculate ET0_1
    coalSum = 0;
    for(i = 0; i < numCoalTimes; i++)
    {
        initIndivs[0] = initIndivs[1] = indivs[1];
        coalSum += CoalPedigree_sim_coal_time(ped, initIndivs);
    }
    ETw1 = (double)coalSum/(double)numCoalTimes;
    
    // finally calculate ET1 (between)
    coalSum = 0;
    for(i = 0; i < numCoalTimes; i++)
    {
        initIndivs[0] = indivs[0];
        initIndivs[1] = indivs[1];
        coalSum += CoalPedigree_sim_coal_time(ped, initIndivs);
    }
    ETb = (double)coalSum/(double)numCoalTimes;

    fst = 1.0 - (ETw0 + ETw1)/2.0/ETb;
    
    return fst;
}

double simulate_mean_coal_time(CoalPedigree * ped, int32_t * indivs, int32_t numCoalTimes)
{
    int32_t rep;
    int32_t initIndivs[2];
    int32_t sum = 0;
    double mean;
    for(rep = 0; rep < numCoalTimes; rep++)
    {
        initIndivs[0] = indivs[0];
        initIndivs[1] = indivs[1];
        sum += CoalPedigree_sim_coal_time(ped, initIndivs);
    }
    mean = (double)sum/(double)numCoalTimes;
    return mean;
}

void print_many_means(CoalPedigree * ped, Options * opt, int32_t numCoalTimes, FILE * fout)
{
    int32_t demeSize = ped->popnSize/ped->numDemes;
    int32_t i, j, k, rep;
    int32_t numPossiblePairs = 0; 
    if(ped->numDemes == 1)
        numPossiblePairs = ped->popnSize*(ped->popnSize-1)/2;
    else if(ped->numDemes > 1 && opt->withinSet)
        numPossiblePairs = ped->popnSize/ped->numDemes*(ped->popnSize/ped->numDemes-1)/2;
    else if(ped->numDemes > 1 && !opt->withinSet)
        numPossiblePairs = ped->popnSize/ped->numDemes * ped->popnSize/ped->numDemes;
    int32_t indivs[2];
    if(opt->numPairs >= numPossiblePairs)
    {
        // print all possible pairs
        if(opt->withinSet || ped->numDemes == 1)
        {
            // print only within one deme
            for(i = 0; i < demeSize-1; i++)
            {
                indivs[0] = i;
                for(j = i+1; j < demeSize; j++)
                {
                    indivs[1] = j;
                    fprintf(fout, "%g\n", simulate_mean_coal_time(ped, indivs, numCoalTimes));
                }
            }
        }
        else // !opt->withinSet && ped->numDemes > 1 (between)
        {
            for(i = 0; i < demeSize-1; i++)
            {
                indivs[0] = i;
                for(j = demeSize; j < 2*demeSize-1; j++)
                {
                    indivs[1] = j;
                    fprintf(fout, "%g\n", simulate_mean_coal_time(ped, indivs, numCoalTimes));
                }
            }
        }
    }

    else if(opt->numPairs < numPossiblePairs)
    {
        int32_t * pairsIdxs = (int32_t *)malloc(sizeof(int32_t) * opt->numPairs);
        CHECKPOINTER(pairsIdxs);

        int32_t * pairsBucket = (int32_t *)malloc(sizeof(int32_t) * numPossiblePairs);

        equal_sample_noreplace(opt->numPairs, numPossiblePairs, pairsIdxs, pairsBucket);
        free(pairsBucket);

        if(opt->withinSet || ped->numDemes == 1)
        {
            int32_t curPairIdx = 0;
            for(k = 0; k < opt->numPairs; k++)
            {
                curPairIdx = 0;
                for(i = 0; i < demeSize-1; i++)
                {
                    for(j = 0; j < demeSize; j++)
                    {
                        if(pairsIdxs[k] == curPairIdx)
                        {
                            indivs[0] = i;
                            indivs[1] = j;
                            fprintf(fout, "%g\n", simulate_mean_coal_time(ped, indivs, numCoalTimes));
                        }
                        curPairIdx++;
                    }
                }
            }
        }

        else if(!opt->withinSet)
        {
            for(i = 0; i < opt->numPairs; i++)
            {
                indivs[0] = pairsIdxs[i] / demeSize;
                indivs[1] = pairsIdxs[i] % demeSize + demeSize;
                fprintf(fout, "%g\n", simulate_mean_coal_time(ped, indivs, numCoalTimes));
            }
        }
        free(pairsIdxs);
    }
    return;                    
}

void print_many_means_with_rearrs(CoalPedigree * ped, Options * opt, int32_t numCoalTimes, FILE * fout)
{
    int32_t demeSize = ped->popnSize/ped->numDemes;
    int32_t i, j, k, rep;
    int32_t numPossiblePairs = 0; 
    if(ped->numDemes == 1)
        numPossiblePairs = ped->popnSize*(ped->popnSize-1)/2;
    else if(ped->numDemes > 1 && opt->withinSet)
        numPossiblePairs = ped->popnSize/ped->numDemes*(ped->popnSize/ped->numDemes-1)/2;
    else if(ped->numDemes > 1 && !opt->withinSet)
        numPossiblePairs = ped->popnSize/ped->numDemes * ped->popnSize/ped->numDemes;
    int32_t indivs[2];

    if(opt->numPairs >= numPossiblePairs)
    {
        // print all possible pairs
        if(opt->withinSet || ped->numDemes == 1)
        {
            // print only within one deme
            for(i = 0; i < demeSize-1; i++)
            {
                indivs[0] = i;
                for(j = i+1; j < demeSize; j++)
                {
                    indivs[1] = j;
                    fprintf(fout, "%g\n", simulate_mean_coal_time(ped, indivs, numCoalTimes));
                    if(opt->numDemes == 1)
                        Rearrangement_calculate_reconfiguration_probs_diploid(ped, indivs, 2, opt->rearrGen, fout);
                    else if(opt->numDemes == 2)
                        Rearrangement_calculate_reconfiguration_probs_diploid(ped, indivs, 2, opt->rearrGen, fout);
                }
            }
        }
        else // !opt->withinSet && ped->numDemes > 1 (between)
        {
            for(i = 0; i < demeSize; i++)
            {
                indivs[0] = i;
                for(j = demeSize; j < 2*demeSize; j++)
                {
                    indivs[1] = j;
                    fprintf(fout, "%g\n", simulate_mean_coal_time(ped, indivs, numCoalTimes));
                    if(opt->numDemes == 1)
                        Rearrangement_calculate_reconfiguration_probs_diploid(ped, indivs, 2, opt->rearrGen, fout);
                    else if(opt->numDemes == 2)
                        Rearrangement_calculate_reconfiguration_probs_diploid(ped, indivs, 2, opt->rearrGen, fout);
                }
            }
        }
    }

    else if(opt->numPairs < numPossiblePairs)
    {
        int32_t * pairsIdxs = (int32_t *)malloc(sizeof(int32_t) * opt->numPairs);
        CHECKPOINTER(pairsIdxs);

        int32_t * pairsBucket = (int32_t *)malloc(sizeof(int32_t) * numPossiblePairs);

        equal_sample_noreplace(opt->numPairs, numPossiblePairs, pairsIdxs, pairsBucket);
        free(pairsBucket);

        int * indicesi = (int *)malloc(sizeof(int) * (size_t)numPossiblePairs);
        CHECKPOINTER(indicesi);
        int * indicesj = (int *)malloc(sizeof(int) * (size_t)numPossiblePairs);
        CHECKPOINTER(indicesj);


        if(opt->withinSet || ped->numDemes == 1)
        {
            int curPairIdx = 0;
            for(i = 0; i < demeSize-1; i++)
            {
                for(j = i+1; j < demeSize; j++)
                {
                    indicesi[curPairIdx] = i;
                    indicesj[curPairIdx] = j;
                    curPairIdx++;
                }
            }

            for(k = 0; k < opt->numPairs; k++)
            {
                indivs[0] = indicesi[pairsIdxs[k]];
                indivs[1] = indicesj[pairsIdxs[k]];
                fprintf(fout, "%g\n", simulate_mean_coal_time(ped, indivs, numCoalTimes));
                if(opt->numDemes == 1)
                    Rearrangement_calculate_reconfiguration_probs_diploid(ped, indivs, 2, opt->rearrGen, fout);
                else if(opt->numDemes == 2)
                    Rearrangement_calculate_reconfiguration_probs_diploid(ped, indivs, 2, opt->rearrGen, fout);
            }
        }

        else if(!opt->withinSet)
        {
            for(i = 0; i < opt->numPairs; i++)
            {
                indivs[0] = pairsIdxs[i] / demeSize;
                indivs[1] = pairsIdxs[i] % demeSize + demeSize;
                fprintf(fout, "%g\n", simulate_mean_coal_time(ped, indivs, numCoalTimes));
                if(opt->numDemes == 1)
                    Rearrangement_calculate_reconfiguration_probs_diploid(ped, indivs, 2, opt->rearrGen, fout);
                else if(opt->numDemes == 2)
                    Rearrangement_calculate_reconfiguration_probs_diploid(ped, indivs, 2, opt->rearrGen, fout);
            }
        }
        free(pairsIdxs);
        free(indicesi);
        free(indicesj);
    }
    return;                    
}

void print_many_distributions(CoalPedigree * ped, Options * opt, FILE * fout)
{
    int32_t demeSize = ped->popnSize/ped->numDemes;
    int32_t i, j, k, rep;
    int32_t numPossiblePairs = 0; 
    if(ped->numDemes == 1)
        numPossiblePairs = ped->popnSize*(ped->popnSize-1)/2;
    else if(ped->numDemes > 1 && opt->withinSet)
        numPossiblePairs = ped->popnSize/ped->numDemes*(ped->popnSize/ped->numDemes-1)/2;
    else if(ped->numDemes > 1 && !opt->withinSet)
        numPossiblePairs = ped->popnSize/ped->numDemes * ped->popnSize/ped->numDemes;
    int32_t indivs[2];
    if(opt->numPairs >= numPossiblePairs)
    {
        // print all possible pairs
        if(opt->withinSet || ped->numDemes == 1)
        {
            // print only within one deme
            for(i = 0; i < demeSize-1; i++)
            {
                indivs[0] = i;
                for(j = i+1; j < demeSize; j++)
                {
                    indivs[1] = j;
                    print_coal_probs(ped, indivs, ped->numGenerations, fout);
                    fprintf(fout, "---\n");
                }
            }
        }
        else // !opt->withinSet && ped->numDemes > 1 (between)
        {
            for(i = 0; i < demeSize-1; i++)
            {
                indivs[0] = i;
                for(j = demeSize; j < 2*demeSize-1; j++)
                {
                    indivs[1] = j;
                    print_coal_probs(ped, indivs, ped->numGenerations, fout);
                    fprintf(fout, "---\n");
                }
            }
        }
    }

    else if(opt->numPairs < numPossiblePairs)
    {
        int32_t * pairsIdxs = (int32_t *)malloc(sizeof(int32_t) * opt->numPairs);
        CHECKPOINTER(pairsIdxs);

        int32_t * pairsBucket = (int32_t *)malloc(sizeof(int32_t) * numPossiblePairs);

        equal_sample_noreplace(opt->numPairs, numPossiblePairs, pairsIdxs, pairsBucket);
        free(pairsBucket);

        if(opt->withinSet || ped->numDemes == 1)
        {
            int32_t curPairIdx = 0;
            for(k = 0; k < opt->numPairs; k++)
            {
                curPairIdx = 0;
                for(i = 0; i < demeSize-1; i++)
                {
                    for(j = 0; j < demeSize; j++)
                    {
                        if(pairsIdxs[k] == curPairIdx)
                        {
                            indivs[0] = i;
                            indivs[1] = j;
                            print_coal_probs(ped, indivs, ped->numGenerations, fout);
                            fprintf(fout, "---\n");
                        }
                        curPairIdx++;
                    }
                }
            }
        }

        else if(!opt->withinSet)
        {
            for(i = 0; i < opt->numPairs; i++)
            {
                indivs[0] = pairsIdxs[i] / demeSize;
                indivs[1] = pairsIdxs[i] % demeSize + demeSize;
                print_coal_probs(ped, indivs, ped->numGenerations, fout);
                fprintf(fout, "---\n");
            }
        }
        free(pairsIdxs);
    }
    return;                    
}

void print_many_distributions_with_rearrs(CoalPedigree * ped, Options * opt, FILE * fout)
{
    int32_t demeSize = ped->popnSize/ped->numDemes;
    int32_t i, j, k, rep;
    int32_t numPossiblePairs = 0; 
    if(ped->numDemes == 1)
        numPossiblePairs = ped->popnSize*(ped->popnSize-1)/2;
    else if(ped->numDemes > 1 && opt->withinSet)
        numPossiblePairs = ped->popnSize/ped->numDemes*(ped->popnSize/ped->numDemes-1)/2;
    else if(ped->numDemes > 1 && !opt->withinSet)
        numPossiblePairs = ped->popnSize/ped->numDemes * ped->popnSize/ped->numDemes;
    int32_t indivs[2];

    if(opt->numPairs >= numPossiblePairs)
    {
        // print all possible pairs
        if(opt->withinSet || ped->numDemes == 1)
        {
            // print only within one deme
            for(i = 0; i < demeSize-1; i++)
            {
                indivs[0] = i;
                for(j = i+1; j < demeSize; j++)
                {
                    indivs[1] = j;
                    print_coal_probs(ped, indivs, ped->numGenerations, fout);
                    if(opt->numDemes == 1)
                        Rearrangement_calculate_reconfiguration_probs_diploid(ped, indivs, 2, opt->rearrGen, fout);
                    else if(opt->numDemes == 2)
                        Rearrangement_calculate_reconfiguration_probs_diploid(ped, indivs, 2, opt->rearrGen, fout);
                    fprintf(fout, "---\n");
                }
            }
        }
        else // !opt->withinSet && ped->numDemes > 1 (between)
        {
            for(i = 0; i < demeSize; i++)
            {
                indivs[0] = i;
                for(j = demeSize; j < 2*demeSize; j++)
                {
                    indivs[1] = j;
                    print_coal_probs(ped, indivs, ped->numGenerations, fout);
                    if(opt->numDemes == 1)
                        Rearrangement_calculate_reconfiguration_probs_diploid(ped, indivs, 2, opt->rearrGen, fout);
                    else if(opt->numDemes == 2)
                        Rearrangement_calculate_reconfiguration_probs_diploid(ped, indivs, 2, opt->rearrGen, fout);
                    fprintf(fout, "---\n");
                }
            }
        }
    }

    else if(opt->numPairs < numPossiblePairs)
    {
        int32_t * pairsIdxs = (int32_t *)malloc(sizeof(int32_t) * opt->numPairs);
        CHECKPOINTER(pairsIdxs);

        int32_t * pairsBucket = (int32_t *)malloc(sizeof(int32_t) * numPossiblePairs);

        equal_sample_noreplace(opt->numPairs, numPossiblePairs, pairsIdxs, pairsBucket);
        free(pairsBucket);

        int * indicesi = (int *)malloc(sizeof(int) * (size_t)numPossiblePairs);
        CHECKPOINTER(indicesi);
        int * indicesj = (int *)malloc(sizeof(int) * (size_t)numPossiblePairs);
        CHECKPOINTER(indicesj);


        if(opt->withinSet || ped->numDemes == 1)
        {
            int curPairIdx = 0;
            for(i = 0; i < demeSize-1; i++)
            {
                for(j = i+1; j < demeSize; j++)
                {
                    indicesi[curPairIdx] = i;
                    indicesj[curPairIdx] = j;
                    curPairIdx++;
                }
            }

            for(k = 0; k < opt->numPairs; k++)
            {
                indivs[0] = indicesi[pairsIdxs[k]];
                indivs[1] = indicesj[pairsIdxs[k]];
                print_coal_probs(ped, indivs, ped->numGenerations, fout);
                if(opt->numDemes == 1)
                    Rearrangement_calculate_reconfiguration_probs_diploid(ped, indivs, 2, opt->rearrGen, fout);
                else if(opt->numDemes == 2)
                    Rearrangement_calculate_reconfiguration_probs_diploid(ped, indivs, 2, opt->rearrGen, fout);
                fprintf(fout, "---\n");
            }
        }

        else if(!opt->withinSet)
        {
            for(i = 0; i < opt->numPairs; i++)
            {
                indivs[0] = pairsIdxs[i] / demeSize;
                indivs[1] = pairsIdxs[i] % demeSize + demeSize;
                print_coal_probs(ped, indivs, ped->numGenerations, fout);
                if(opt->numDemes == 1)
                    Rearrangement_calculate_reconfiguration_probs_diploid(ped, indivs, 2, opt->rearrGen, fout);
                else if(opt->numDemes == 2)
                    Rearrangement_calculate_reconfiguration_probs_diploid(ped, indivs, 2, opt->rearrGen, fout);
                fprintf(fout, "---\n");
            }
        }
        free(pairsIdxs);
        free(indicesi);
        free(indicesj);
    }
    return;                    
}

void print_prob_other_deme(CoalPedigree * ped, int indiv)
{
    assert(ped->numDemes == 2);
    int t, curIndiv, i;

    const int demeSize = ped->popnSize / 2;
    const int popnSize = ped->popnSize;

    const int otherDemeIs1 = indiv < demeSize;

    double sum;

    double * probs = (double *)malloc(sizeof(double) * ped->popnSize);
    CHECKPOINTER(probs);
    double * probsP = (double *)malloc(sizeof(double) * ped->popnSize);
    CHECKPOINTER(probsP);

    for(i = 0; i < ped->popnSize; i++)
    {
        probs[i] = 0.0;
        probsP[i] = 0.0;
    }
    probs[indiv] = 1.0;

    for(t = 0; t < ped->numGenerations; t++)
    {
        sum = 0.0;
        if(otherDemeIs1)
        {
            for(i = demeSize; i < popnSize; i++)
                sum += probs[i];
        }
        else //  other deme is 0
        {
            for(i = 0; i < demeSize; i++)
                sum += probs[i];
        }
        printf("%i\t%e\n", t, sum);

        // now advance the probability
        for(i = 0; i < popnSize; i++)
        {
            probsP[ped->relationships[t][i][0]] += 0.5 * probs[i];
            probsP[ped->relationships[t][i][1]] += 0.5 * probs[i];
        }
        for(i = 0; i < popnSize; i++)
        {
            probs[i] = probsP[i];
            probsP[i] = 0.0;
        } 
    }
    free(probs);
    free(probsP);
    return;
}





int32_t main(int argc, char * argv[])
{
	CoalPedigree ped;

    Options opt;

	//char optString[] = "D:N:M:se";

	randseed();

    Options_parse_options_pedsim(argc, argv, &opt);

    int32_t indivs[2];

	CoalPedigree_init(&ped, &opt);

    get_init_indivs_pedsim(&ped, &opt, &(indivs[0]));
    opt.indivs = indivs;

	CoalPedigree_shuffle(&ped, &opt);

	if(opt.exportSet)
		CoalPedigree_export_pedigree(&ped, opt.exportFilename);

    if(opt.singleLineageSet)
    {
        assert(ped.numDemes == 2);
        const int demeSize = ped.popnSize / 2;
        int i, indiv;
        for(i = 0; i < opt.numSingleLineages-1; i++)
        {
            indiv = runifd(0, demeSize-1);
            print_prob_other_deme(&ped, indiv);
            printf("---\n");
        }
        indiv = runifd(0, demeSize-1);
        print_prob_other_deme(&ped, indiv);
        exit(0);
    }

	// rearrangement probabilities
	if(opt.rearrSet && !opt.multiplePairs)
	{
		FILE * rearrOut = fopen(opt.rearrFilename, "w");
		if(!rearrOut)
			PERROR("Could not open file to write rearrangement probability.");
		if(opt.numDemes == 1)
            Rearrangement_calculate_reconfiguration_probs_diploid(&ped, indivs, 2, opt.rearrGen, rearrOut);
        else if(opt.numDemes == 2)
            Rearrangement_calculate_reconfiguration_probs_diploid(&ped, indivs, 2, opt.rearrGen, rearrOut);
        else
            PERROR("Outputting rearrangements is supported for only 1 or 2 demes.");
		fclose(rearrOut);
	}

    // now simulate or calculate exact distribution
    if(!opt.fstSet && !opt.multiplePairs)
    {
        if(opt.exactSet) // --exact
            print_coal_probs(&ped, indivs, ped.numGenerations, stdout);
        else  //--means [default]
            print_sim_coal_times(&ped, indivs, opt.numCoalTimes, stdout);
    }
    else if(opt.fstSet && !opt.multiplePairs) // --fst
        printf("%g\n", simulate_fst(&ped, indivs, opt.numCoalTimes));
    else if(opt.multiplePairs && opt.numPairs > 1 && !opt.exactSet)
    {
        if(!opt.rearrSet)
            print_many_means(&ped, &opt, opt.numCoalTimes, stdout);
        else
            print_many_means_with_rearrs(&ped, &opt, opt.numCoalTimes, stdout);
    }
    else if(opt.multiplePairs && opt.numPairs > 1 && opt.exactSet)
    {
        if(!opt.rearrSet)
            print_many_distributions(&ped, &opt, stdout);
        else
            print_many_distributions_with_rearrs(&ped, &opt, stdout);
    }

	CoalPedigree_free(&ped);
	return 0;
}

