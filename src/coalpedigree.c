#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "definitions.h"
#include "random.h"
#include "options.h"
#include "coalpedigree.h"

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
 * ped          pointer to uninitialized pedigree
 * numGens      number of generations to simulate. enough memory is allocated to have numGens generations of relationships
 * numFemales   number of females in the population
 * numMales     number of males in the population
 *
*/
void CoalPedigree_init_singledeme(CoalPedigree * ped, int32_t numGens, int32_t numFemales, int32_t numMales)
{
	int32_t i, t;

	ped->numFemales = numFemales;
	ped->numMales = numMales;
	ped->numGenerations = numGens;
	ped->popnSize = numFemales + numMales;
	ped->variablePopnSize = 0;
	ped->numDemes = 1;
	ped->relationships = (twoints **)malloc(sizeof(twoints *) * (size_t)numGens);
	ped->popnSizes = NULL;
	for(t = 0; t < numGens; t++)
	{
		ped->relationships[t] = (twoints *)malloc(sizeof(twoints) * (size_t)(ped->popnSize));
		CHECKPOINTER(ped->relationships[t]);
	}
}

void CoalPedigree_init_ddemes(CoalPedigree * ped, int32_t numGens, int32_t numDemes, int32_t demeSize)
{
	int32_t i, t;

	if(demeSize % 2 != 0)
		PERROR("Odd-numbered deme size in CoalPedigree_init_ddemes. Must be even-numbered.");

	ped->numFemales = numDemes * demeSize / 2;
	ped->numMales = numDemes * demeSize / 2;
	ped->numGenerations = numGens;
	ped->popnSize = numDemes * demeSize;
	ped->variablePopnSize = 0;
	ped->popnSizes = NULL;
	ped->relationships = (twoints **)malloc(sizeof(twoints *) * (size_t)numGens);
	CHECKPOINTER(ped->relationships);
	for(t = 0; t < numGens; t++)
	{
		ped->relationships[t] = (twoints *)malloc(sizeof(twoints) * (size_t)(ped->popnSize));
		CHECKPOINTER(ped->relationships[t]);
	}
	ped->numDemes = numDemes;
}


/* Initializes a CoalPedigree struct when the population size varies
 *
 * ped         pointer to CoalPedigree struct to be initialized
 * numGens     number of generations spanned in the pedigree
 * popnSizes   size of the population in each generation, vector of length numGens
 */
void CoalPedigree_init_varN(CoalPedigree * ped, int32_t numGens, int32_t * popnSizes)
{
	int32_t i, t, maxPopnSize;

	ped->numDemes = 1;
	ped->numFemales = 0;
	ped->numMales = 0;
	ped->numGenerations = numGens;
	ped->variablePopnSize = 1;
	ped->popnSizes = (int32_t *)malloc(sizeof(int32_t) * numGens); // Allocate an array of population sizes.
	CHECKPOINTER(ped->popnSizes);
	maxPopnSize = 0;
	for(i = 0; i < numGens; i++)
	{
		ped->popnSizes[i] = popnSizes[i];
		if(popnSizes[i] > maxPopnSize)
			maxPopnSize = popnSizes[i];
	}
	ped->popnSize = maxPopnSize;
	ped->relationships = (twoints **)malloc(sizeof(twoints *) * (size_t)numGens);
	for(t = 0; t < numGens; t++)
		ped->relationships[t] = (twoints *)malloc(sizeof(twoints) * (size_t)(ped->popnSizes[t]));
}

void CoalPedigree_init(CoalPedigree * ped, Options * opt)
{
	int32_t numFemales, numMales, success;
	FILE * fin;
	char * line;

	// import a pedigree
	if(opt->importSet)
	{
		CoalPedigree_import_pedigree(ped, opt->importFilename);
		return;
	}

	// if opt->numGens wasn't specified
	if(opt->numGens == 0)
	{
		if(!(opt->multiplePopnSizes))
		{
			if(opt->numDemes == 1 && !(opt->multiplePopnSizes))
				opt->numGens = 20 * opt->demeSize;
			else if(opt->numDemes > 1)
				opt->numGens = 20 * opt->demeSize * opt->numDemes * (1.0+((double)opt->numDemes - 1.0)/(double)(opt->migRate * opt->numDemes));
		}
		else // opt->multiplePopnSizes == 1
		{
			char ch;
			fin = fopen(opt->popnSizeFile, "r");
			if(!fin)
			{
				fprintf(stderr, "Could not open %s for reading population sizes.\n", opt->popnSizeFile);
				exit(1);
			}
			// count number of lines / population sizes
			opt->numGens = 0;  // already set at zero...
			while(!feof(fin))
			{
				ch = fgetc(fin);
				if(ch == '\n')
					opt->numGens++;
			}
			fclose(fin);
		}
	}

	if(opt->numDemes == 1 && !(opt->multiplePopnSizes))
	{
		numFemales = opt->demeSize / 2 + opt->demeSize % 2; 
		numMales = opt->demeSize / 2;
		CoalPedigree_init_singledeme(ped, opt->numGens, numFemales, numMales);
	}

	else if(opt->numDemes == 1 && opt->multiplePopnSizes)
	{
		int32_t MAX_NUM_POPN_SIZES = 100000;
		int32_t DEFAULT_MAX_LINE_SIZE = 100;
		int32_t * popnSizesTemp = (int32_t *)malloc(sizeof(int32_t) * MAX_NUM_POPN_SIZES);
		CHECKPOINTER(popnSizesTemp);
		int32_t popnSizeIdx = 0;
		line = (char *)malloc(sizeof(char) * DEFAULT_MAX_LINE_SIZE);
		fin = fopen(opt->popnSizeFile, "r");
		if(!fin)
		{
			fprintf(stderr, "Could not open %s for reading population sizes.\n", opt->popnSizeFile);
			exit(1);
		}
		while(fgets(line, DEFAULT_MAX_LINE_SIZE, fin) != NULL)
		{
			success = sscanf(line, "%i\n", &(popnSizesTemp[popnSizeIdx++]));
			if(!success)
				PERROR("Invalid population size in --popn-size-file.");
		}
		CoalPedigree_init_varN(ped, opt->numGens, popnSizesTemp);
		free(popnSizesTemp);
		free(line);
		fclose(fin);
	}
	else
		CoalPedigree_init_ddemes(ped, opt->numGens, opt->numDemes, opt->demeSize);
	return;
}

void CoalPedigree_free(CoalPedigree * ped)
{
	int32_t t,i;
	for(t = 0; t < ped->numGenerations; t++)
		free(ped->relationships[t]);
	if(ped->popnSizes) // NULL if numDemes == 1
		free(ped->popnSizes);
	free(ped->relationships);
}

void CoalPedigree_shuffle(CoalPedigree * ped, Options * opt)
{
	if(ped->numDemes == 1)
	{
		if(opt->cyclical)
		{
			CoalPedigree_shuffle_cyc_WF_dio(ped);
			return;
		}
		if(opt->sweepSet)
		{
			CoalPedigree_shuffle_sweep_dio(ped, opt->sweep_ts, opt->sweep_s, opt->sweep_h, opt->requireFixation);
			return;
		}
		if(opt->monoecious)
		{
			CoalPedigree_shuffle_monoecious(ped);
			return;
		}
		if(opt->multiplePopnSizes)
		{
			CoalPedigree_shuffle_dioecious_varN(ped, opt->loopGen);
			return;
		}
        if(opt->fixedIBDset)
        {
            int32_t i;
            // first individual in opt->indivs is the IBD individual; the rest have no IBD
            opt->fixedIBDindiv = opt->indivs[0];
            int32_t * nonIBDindivs = (int32_t *)malloc(sizeof(int32_t) * opt->sampleSize-1);
            CHECKPOINTER(nonIBDindivs);
            for(i = 1; i < opt->sampleSize; i++)
                nonIBDindivs[i-1] = opt->indivs[i];
            int32_t clearGen = opt->fixedIBDgen;
            CoalPedigree_shuffle_dioecious_fixed_ibd(ped, opt->fixedIBDgen, clearGen, opt->sampleSize, opt->fixedIBDindiv, nonIBDindivs);
            free(nonIBDindivs);
            return;
        }
        if(opt->sharedAncSet)
        {
            opt->sharedAncIndivs[0] = opt->indivs[0];
            opt->sharedAncIndivs[1] = opt->indivs[1];
            int32_t clearGen = opt->sharedAncGen;
            CoalPedigree_shuffle_dioecious_shared_anc(ped, opt->sharedAncGen, clearGen, opt->indivs, opt->sampleSize);
            return;
        }
		// else...
		CoalPedigree_shuffle_dioecious(ped);
		return;
	}
	else // ped->numDemes != 1
	{
		double avgNumMigrants = opt->migRate/4.0;
		if(ped->numDemes > 2)
		{
			CoalPedigree_shuffle_ddeme_diploid_stochastic(ped, avgNumMigrants);
			return;
		}
		else if(ped->numDemes == 2)
		{
            if(opt->fixedMigrantSet)
            {
                opt->fixedMigrantIndiv = opt->indivs[0];
                //int32_t clearGen = 4;
                // now going to bifurcate pedigree only until the generation including the migrant ancestor
                int32_t clearGen = opt->fixedMigrantGen;
                
                CoalPedigree_shuffle_twopop_fixed_migrant(ped, avgNumMigrants, opt->indivs, opt->sampleSize*2, opt->fixedMigrantGen, clearGen);
            }
            else
                CoalPedigree_shuffle_twopop_diploid_fully_stochastic(ped, avgNumMigrants);
			return;
		}
		else
			PERROR("Invalid ped->numDemes in CoalPedigree_shuffle().");
	}
	PERROR("No shuffling in CoalPedigree_shuffle().");
	return;
}



void CoalPedigree_shuffle_dioecious(CoalPedigree * ped)
{
	int32_t t,i;
	const int32_t numGens = ped->numGenerations;
	const int32_t numFemales = ped->numFemales;
	const int32_t popnSize = ped->popnSize;
	const int32_t * popnSizes = ped->popnSizes;
	int32_t currentN, parentN, parentNumFemales, parentNumMales;

	if(!(ped->variablePopnSize))
	{
		for(t = 0; t < numGens; t++)
		{
			for(i = 0; i < popnSize; i++)
			{
				ped->relationships[t][i][0] = runifd(0,numFemales-1);
				ped->relationships[t][i][1] = runifd(numFemales,popnSize-1);
			}
		}
		return;
	}
	else
	{
		for(t = 0; t < numGens-1; t++)
		{
			currentN = popnSizes[t];
			parentN = popnSizes[t+1];
			parentNumFemales = parentN / 2;
			parentNumMales = parentN / 2;
			/* Pick an integer randomly numF times from 0, ... , numM, and fill the
			 * given array with those numbers. */
			for(i = 0; i < currentN; i++)
			{
				// runifd is inclusive
				ped->relationships[t][i][0] = runifd(0,parentNumFemales-1);
				ped->relationships[t][i][1] = runifd(parentNumFemales,parentN-1);
			}
		}
		// Assuming the population size doesn't change between numGens-1 and numGens, since we don't have info for generation numGens.
		for(i = 0; i < currentN; i++)
		{
			currentN = parentN = popnSizes[numGens-1];
			parentNumFemales = parentN / 2;
			parentNumMales = parentN / 2;
			ped->relationships[numGens-1][i][0] = runifd(0,parentNumFemales-1);
			ped->relationships[numGens-1][i][1] = runifd(parentNumFemales,parentN-1);
		}
	}
	return;
}

/* gives the two individuals in indivs an ancestor in relAncGen */
void CoalPedigree_shuffle_dioecious_fixed_recent_relative(CoalPedigree * ped, int32_t * indivs, int32_t relAncGen)
{
	int32_t t,i,j, ancestor, ancestorP, ancestor1, ancestorP1, sex;
	const int32_t indiv = indivs[0];
	const int32_t N = ped->popnSize;
	const int32_t firstFemale = 0;
	const int32_t lastFemale = ped->popnSize/2-1;
	const int32_t firstMale = ped->popnSize/2;
	const int32_t lastMale = ped->popnSize-1;
	double migWeight;
	CoalPedigree_shuffle_dioecious(ped);
	ancestor = ped->relationships[0][indiv][1];		// the 1-chromosome is the inbred/migrant one.
	for(t = 1; t < relAncGen; t++)
	{
		ancestorP = ped->relationships[t][ancestor][rbern(0.5)];
		ancestor = ancestorP;
	}
	ancestor1 = indivs[1];
	for(t = 0; t < relAncGen-1; t++)
	{
		ancestorP1 = ped->relationships[t][ancestor1][rbern(0.5)];
		ancestor1 = ancestorP1;
	}
	sex = (ancestor >= firstMale);
	ped->relationships[relAncGen-1][ancestor1][sex] = ancestor;
	return;
}

void CoalPedigree_recursive_binary_tree_(CoalPedigree * ped, int32_t curGen, int32_t targetGen, int32_t curIndiv, int32_t parentIdx)
{
    if(curGen == targetGen)
        return;
    const int32_t firstMale = ped->popnSize/2;
    ped->relationships[curGen][curIndiv][0] = parentIdx;
    ped->relationships[curGen][curIndiv][1] = firstMale + parentIdx;
    CoalPedigree_recursive_binary_tree_(ped, curGen+1, targetGen, parentIdx, parentIdx*2);
    CoalPedigree_recursive_binary_tree_(ped, curGen+1, targetGen, firstMale+parentIdx, parentIdx*2+1);
    return;
}

/* shuffles a pedigree to give ibdIndiv's parents a relationship in generation
 * ibdGen, while the individuals in nonIBDindivs have no identity by descent
 * (perfect binary tree pedigrees)
 *
 * ped            initialized single-deme pedigree
 * ibdGen         generation in which ibdIndiv's parents share an ancestor (number of generations before present)
 * clearGen       generation until which all the individuals' ancestry is a perfect binary tree (no inbreeding)
 * numIndivs      number of individuals in teh sample
 * ibdIndiv       individual that will get the IBD
 * nonIBDindivs   other individuals in the sample, which need perfect binary tree
 * */
void CoalPedigree_shuffle_dioecious_fixed_ibd(CoalPedigree * ped, int32_t ibdGen, int32_t clearGen, int32_t numIndivs, int32_t ibdIndiv, int32_t * nonIBDindivs)
{
    // make sure there are enough individuals in the pedigree
    if(numIndivs * pow((double)2, (double)ibdGen)  > (double)(ped->popnSize))
        PERROR("Too many individuals for perfect binary tree in CoalPedigree_shuffle_dioecious_fixed_ibd()");
    CoalPedigree_shuffle_dioecious(ped);
    int32_t i, t;
    
    // Make a perfect binary tree for all of the individuals, including ibdIndiv
    // ibdIndiv will be the one whose parents are the first females and males
    // in the first generation back in time.
    CoalPedigree_recursive_binary_tree_(ped, 0, clearGen, ibdIndiv, 0);
    for(i = 1; i < numIndivs; i++)
        CoalPedigree_recursive_binary_tree_(ped, 0, clearGen, nonIBDindivs[i-1], i);

    // add the inbred relationship for ibdIndiv
    int32_t leftIndiv = ped->relationships[0][ibdIndiv][0];
    int32_t rightIndiv = ped->relationships[0][ibdIndiv][1];
    for(t = 1; t < ibdGen-1; t++)
    {
        leftIndiv = ped->relationships[t][leftIndiv][0];
        rightIndiv = ped->relationships[t][rightIndiv][0];
    }
    ped->relationships[ibdGen-1][rightIndiv][0] = ped->relationships[ibdGen-1][leftIndiv][0]; 
    return;
}

// indivs[0] and indivs[1] are the two that are going to shared the ancestor
void CoalPedigree_shuffle_dioecious_shared_anc(CoalPedigree * ped, int32_t sharedAncGen, int32_t clearGen, int32_t * indivs, int32_t numIndivs)
{
    // make sure there are enough individuals in the pedigree
    if(numIndivs * pow((double)2, (double)clearGen)  > (double)(ped->popnSize))
        PERROR("Too many individuals for perfect binary tree in CoalPedigree_shuffle_dioecious_fixed_ibd()");
    CoalPedigree_shuffle_dioecious(ped);
    int32_t i, t;
    
    for(i = 0; i < numIndivs; i++)
        CoalPedigree_recursive_binary_tree_(ped, 0, clearGen, indivs[i], i);

    // add the inbred relationship for indivs[0] and indivs[1]
    int32_t leftIndiv = indivs[0];
    int32_t rightIndiv = indivs[1];
    for(t = 0; t < sharedAncGen-1; t++)
    {
        leftIndiv = ped->relationships[t][leftIndiv][0];
        rightIndiv = ped->relationships[t][rightIndiv][0];
    }
    ped->relationships[sharedAncGen-1][rightIndiv][0] = ped->relationships[sharedAncGen-1][leftIndiv][0]; 
    return;
}


void CoalPedigree_shuffle_twopop_nomigration(CoalPedigree * ped)
{
	int32_t t,i;
	const int32_t numGens = ped->numGenerations;
	const int32_t numFemales = ped->numFemales;
	const int32_t popnSize = ped->popnSize;
	const int32_t demeSize = popnSize/2;

	if(popnSize % 4 != 0)
	{
		fprintf(stderr, "Warning: Shuffling two-population pedigree with total population size that \
				is not a multiple of 4. (in CoalPedigree_shuffle_twopop)\n");
	}

	const int32_t first0Female = 0;
	const int32_t last0Female = demeSize/2-1;
	const int32_t first0Male = demeSize/2;
	const int32_t last0Male = demeSize-1;
	const int32_t first1Female = demeSize;
	const int32_t last1Female = demeSize + demeSize/2-1;
	const int32_t first1Male = demeSize + demeSize/2;
	const int32_t last1Male = popnSize-1;

	for(t = 0; t < numGens; t++)
	{
		/* Pick an integer randomly numF times from 0, ... , numM, and fill the
		 * given array with those numbers. */
		for(i = 0; i < demeSize; i++)
		{
			ped->relationships[t][i][0] = runifd(first0Female,last0Female);
			ped->relationships[t][i][1] = runifd(first0Male,last0Male);
		}
		for(i = first1Female; i < popnSize; i++)
		{
			ped->relationships[t][i][0] = runifd(first1Female,last1Female);
			ped->relationships[t][i][1] = runifd(first1Male,last1Male);
		}
	}
	return;
}

void CoalPedigree_shuffle_ddeme_nomigration(CoalPedigree * ped)
{
	int32_t t, i, j, firstFemale, lastFemale, firstMale, lastMale;
	const int32_t numGens = ped->numGenerations;
	const int32_t popnSize = ped->popnSize;
	const int32_t numDemes = ped->numDemes;
	const int32_t demeSize = popnSize/numDemes;
	const int32_t numFemales = demeSize / 2;
	const int32_t numMales = demeSize / 2;
	if(popnSize/numDemes % 2 != 0)
	{
		fprintf(stderr, "Warning: Shuffling d-deme pedigree with total population size that is not a multiple of 2*numDemes.\n");
	}
	for(t = 0; t < numGens; t++)
	{
		for(i = 0; i < numDemes; i++)
		{
			firstFemale = demeSize * i;
			lastFemale = firstFemale + numFemales - 1;
			firstMale = firstFemale + numFemales;
			lastMale = firstMale + numMales - 1;
			for(j = firstFemale; j <= lastMale; j++)
			{
				ped->relationships[t][j][0] = runifd(firstFemale,lastFemale);
				ped->relationships[t][j][1] = runifd(firstMale,lastMale);
			}
		}
	}
	return;
}

// The following function shuffles two subpopulations without migration over the *inclusive* interval
// [start, end].
void CoalPedigree_shuffle_twopop_nomigration_interval(CoalPedigree * ped, int32_t start, int32_t end)
{
	int32_t t,i;
	const int32_t popnSize = ped->popnSize;
	const int32_t demeSize = popnSize/2;
	if(popnSize % 4 != 0)
	{
		fprintf(stderr, "Warning: Shuffling two-population pedigree with total population size that \
				is not a multiple of 4. (in CoalPedigree_shuffle_twopop)\n");
	}

	const int32_t first0Female = 0;
	const int32_t last0Female = demeSize/2-1;
	const int32_t first0Male = demeSize/2;
	const int32_t last0Male = demeSize-1;
	const int32_t first1Female = demeSize;
	const int32_t last1Female = demeSize + demeSize/2-1;
	const int32_t first1Male = demeSize + demeSize/2;
	const int32_t last1Male = popnSize-1;

	for(t = start; t <= end; t++)
	{
		/* Pick an integer randomly numF times from 0, ... , numM, and fill the
		 * given array with those numbers. */
		for(i = first0Female; i <= last0Male; i++)
		{
			// runifd is inclusive
			ped->relationships[t][i][0] = runifd(first0Female,last0Female);
			ped->relationships[t][i][1] = runifd(first0Male,last0Male);
		}
		for(i = first1Female; i <= last1Male; i++)
		{
			// runifd is inclusive
			ped->relationships[t][i][0] = runifd(first1Female,last1Female);
			ped->relationships[t][i][1] = runifd(first1Male,last1Male);
		}
	}
	return;
}

void CoalPedigree_shuffle_ddeme_diploid_stochastic(CoalPedigree * ped, double avgNumMigrants)
{
	int32_t t, i, j, relationship, numFemaleMigrants, numMaleMigrants, originDeme, firstFemale, lastFemale, firstMale, lastMale;
	const int32_t numGens = ped->numGenerations;
	const int32_t popnSize = ped->popnSize;
	const int32_t numDemes = ped->numDemes;
	const int32_t demeSize = popnSize/numDemes;
	const double m = avgNumMigrants/(double)(demeSize);

	CoalPedigree_shuffle_ddeme_nomigration(ped);
	for(t = 0; t < numGens; t++)
	{
		for(i = 0; i < numDemes; i++)
		{
			firstFemale = i*demeSize;
			lastFemale = i*demeSize+demeSize/2-1;
			firstMale = i*demeSize+demeSize/2;
			lastMale = (i+1)*demeSize-1;
			numFemaleMigrants = gsl_ran_binomial(demeSize/2, m);
			numMaleMigrants = gsl_ran_binomial(demeSize/2, m);
			for(j = firstFemale; j < firstFemale+numFemaleMigrants; j++)
			{
				do
				{
					originDeme = runifd(0, numDemes-1);
				} while(originDeme == i);		// note the migrant's deme of origin cannot be it's own deme.
				ped->relationships[t][j][0] = runifd(originDeme*demeSize, originDeme*demeSize + demeSize/2-1);
				ped->relationships[t][j][1] = runifd(originDeme*demeSize+demeSize/2, (originDeme+1)*demeSize-1);
			}
			for(j = firstMale; j < firstMale+numMaleMigrants; j++)
			{
				do
				{
					originDeme = runifd(0, numDemes-1);
				} while(originDeme == i);
				ped->relationships[t][j][0] = runifd(originDeme*demeSize, originDeme*demeSize + demeSize/2-1);
				ped->relationships[t][j][1] = runifd(originDeme*demeSize+demeSize/2, (originDeme+1)*demeSize-1);
			}
		}
	}
	return;
}

void CoalPedigree_shuffle_twopop_diploid_fully_stochastic(CoalPedigree * coalped, double avgNumMigrants)
{
	int32_t t,i;
	const int32_t numGens = coalped->numGenerations;
	const int32_t numFemales = coalped->numFemales;
	const int32_t popnSize = coalped->popnSize;

	if(popnSize % 4 != 0)
	{
		fprintf(stderr, "Warning: Shuffling two-population pedigree with total population size that \
				is not a multiple of 4. (in CoalPedigree_shuffle_twopop)\n");
	}

	const int32_t first0Female = 0;
	const int32_t last0Female = popnSize/4 - 1;
	const int32_t first0Male = popnSize/4;
	const int32_t last0Male = popnSize/2 - 1;
	const int32_t first1Female = popnSize/2;
	const int32_t last1Female = 3*popnSize/4 - 1;
	const int32_t first1Male = 3*popnSize/4;
	const int32_t last1Male = popnSize - 1;

	double m = avgNumMigrants/(double)(popnSize/2);
	int32_t numMigrants0, numMigrants1;

	// have to decide how to split up the migrants amongst males and females.
	int32_t numFemaleMigrants0, numFemaleMigrants1;
	int32_t numMaleMigrants0, numMaleMigrants1;
	int32_t relationship;


	/* For each generation, pick a father and a mother randomly for all sons and
	 * daughters */
	CoalPedigree_shuffle_twopop_nomigration(coalped);
	for(t = 0; t < numGens; t++)
	{
		numFemaleMigrants0 = gsl_ran_binomial_tpe(popnSize/4, m);
		numMaleMigrants0 = gsl_ran_binomial_tpe(popnSize/4, m);
		for(i = first0Female; i <= last0Male; i++)
		{
			if(coalped->relationships[t][i][0] < first0Female + numFemaleMigrants0)
				coalped->relationships[t][i][0] = runifd(first1Female, last1Female);
			if(coalped->relationships[t][i][1] < first0Male + numMaleMigrants0)
				coalped->relationships[t][i][1] = runifd(first1Male, last1Male);
		}
		numFemaleMigrants1 = gsl_ran_binomial_tpe(popnSize/4, m);
		numMaleMigrants1 = gsl_ran_binomial_tpe(popnSize/4, m);
		for(i = first1Female; i <= last1Male; i++)
		{
			if(coalped->relationships[t][i][0] < first1Female + numFemaleMigrants1)
				coalped->relationships[t][i][0] = runifd(first0Female, last0Female);
			if(coalped->relationships[t][i][1] < first1Male + numMaleMigrants1)
				coalped->relationships[t][i][0] = runifd(first0Male, last0Male);
		}
	}
	return;
}

/* shuffles a two-deme pedigree with constant
 * population sizes and variable migration rates. The
 * migration rates are given in the avgNumMigrants vector,
 * which is of length numGens (== ped->numGenerations). 
 *
 * ped				initialized two-deme pedigree
 * avgNumMigrants	vector containing average number of migrants (Nm) in each generation (length ped->numGens)
 *
 * */
void CoalPedigree_shuffle_twopop_diploid_fully_stochastic_varM(CoalPedigree * ped, double * avgNumMigrants)
{
	int32_t t,i;
	const int32_t numGens = ped->numGenerations;
	const int32_t numFemales = ped->numFemales;
	const int32_t popnSize = ped->popnSize;

	if(popnSize % 4 != 0)
	{
		fprintf(stderr, "Warning: Shuffling two-population pedigree with total population size that \
				is not a multiple of 4. (in CoalPedigree_shuffle_twopop)\n");
	}

	const int32_t first0Female = 0;
	const int32_t last0Female = popnSize/4 - 1;
	const int32_t first0Male = popnSize/4;
	const int32_t last0Male = popnSize/2 - 1;
	const int32_t first1Female = popnSize/2;
	const int32_t last1Female = 3*popnSize/4 - 1;
	const int32_t first1Male = 3*popnSize/4;
	const int32_t last1Male = popnSize - 1;

	double m;
	int32_t numMigrants0, numMigrants1;

	// have to decide how to split up the migrants amongst males and females.
	int32_t numFemaleMigrants0, numFemaleMigrants1;
	int32_t numMaleMigrants0, numMaleMigrants1;
	int32_t relationship;


	/* For each generation, pick a father and a mother randomly for all sons and
	 * daughters */
	CoalPedigree_shuffle_twopop_nomigration(ped);
	for(t = 0; t < numGens; t++)
	{
		m = avgNumMigrants[t]/(double)(popnSize/2);
		numFemaleMigrants0 = gsl_ran_binomial_tpe(popnSize/4, m);
		numMaleMigrants0 = gsl_ran_binomial_tpe(popnSize/4, m);
		for(i = first0Female; i <= last0Male; i++)
		{
			if(ped->relationships[t][i][0] < first0Female + numFemaleMigrants0)
				ped->relationships[t][i][0] = runifd(first1Female, last1Female);
			if(ped->relationships[t][i][1] < first0Male + numMaleMigrants0)
				ped->relationships[t][i][1] = runifd(first1Male, last1Male);
		}
		numFemaleMigrants1 = gsl_ran_binomial_tpe(popnSize/4, m);
		numMaleMigrants1 = gsl_ran_binomial_tpe(popnSize/4, m);
		for(i = first1Female; i <= last1Male; i++)
		{
			if(ped->relationships[t][i][0] < first1Female + numFemaleMigrants1)
				ped->relationships[t][i][0] = runifd(first0Female, last0Female);
			if(ped->relationships[t][i][1] < first1Male + numMaleMigrants1)
				ped->relationships[t][i][0] = runifd(first0Male, last0Male);
		}
	}
	return;
}


void CoalPedigree_shuffle_twopop_diploid_fixed(CoalPedigree * ped, int32_t numMigrants)
{
	int32_t t,i;
	const int32_t numGens = ped->numGenerations;
	const int32_t numFemales = ped->numFemales;
	const int32_t popnSize = ped->popnSize;
	const int32_t demeSize = popnSize/2;

	if(popnSize % 4 != 0)
	{
		fprintf(stderr, "Warning: Shuffling two-population pedigree with total population size that \
				is not a multiple of 4. (in CoalPedigree_shuffle_twopop)\n");
	}

	const int32_t first0Female = 0;
	const int32_t last0Female = demeSize/2-1;
	const int32_t first0Male = demeSize/2;
	const int32_t last0Male = demeSize-1;
	const int32_t first1Female = demeSize;
	const int32_t last1Female = demeSize + demeSize/2-1;
	const int32_t first1Male = demeSize + demeSize/2;
	const int32_t last1Male = popnSize-1;

	int32_t numFemaleMigrants = numMigrants/2;
	int32_t numMaleMigrants = numMigrants - numFemaleMigrants;
	int32_t relationship;

	/* For each generation, pick a father and a mother randomly for all sons and
	 * daughters */
	for(t = 0; t < numGens; t++)
	{
		/* Pick an integer randomly numF times from 0, ... , numM, and fill the
		 * given array with those numbers. */
		for(i = first0Female; i <= last0Male; i++)
		{
			// runifd is inclusive
			relationship = runifd(first0Female, last0Female);
			if(relationship < numFemaleMigrants)
				relationship += demeSize;		// take it over to the second population
			ped->relationships[t][i][0] = relationship;

			relationship = runifd(first0Male, last0Male);
			if(relationship < first0Male + numMaleMigrants)
				relationship += demeSize;
			ped->relationships[t][i][1] = relationship;
		}
		for(i = first1Female; i <= last1Male; i++)
		{
			// runifd is inclusive
			relationship = runifd(first1Female, last1Female);
			if(relationship < first1Female + numFemaleMigrants)
				relationship -= demeSize;
			ped->relationships[t][i][0] = relationship;

			relationship = runifd(last1Female+1, popnSize-1);
			if(relationship < first1Male + numMaleMigrants)
				relationship -= demeSize;
			ped->relationships[t][i][1] = relationship;
		}
	}
	return;
}

void CoalPedigree_shuffle_twopop_diploid_stochastic(CoalPedigree * ped, double avgNumMigrants)
{
	int32_t t,i;
	const int32_t numGens = ped->numGenerations;
	const int32_t numFemales = ped->numFemales;
	const int32_t popnSize = ped->popnSize;
	const int32_t demeSize = popnSize/2;
	if(popnSize % 4 != 0)
	{
		fprintf(stderr, "Warning: Shuffling two-population pedigree with total population size that \
				is not a multiple of 4. (in CoalPedigree_shuffle_twopop)\n");
	}

	const int32_t first0Female = 0;
	const int32_t last0Female = demeSize/2-1;
	const int32_t first0Male = demeSize/2;
	const int32_t last0Male = demeSize-1;
	const int32_t first1Female = demeSize;
	const int32_t last1Female = demeSize + demeSize/2-1;
	const int32_t first1Male = demeSize + demeSize/2;
	const int32_t last1Male = popnSize-1;

	double m = avgNumMigrants/(double)(popnSize/2);
	int32_t numFemaleMigrants;
	int32_t numMaleMigrants;
	int32_t relationship;

	/* For each generation, pick a father and a mother randomly for all sons and
	 * daughters */
	CoalPedigree_shuffle_twopop_nomigration(ped);
	for(t = 0; t < numGens; t++)
	{
		numFemaleMigrants = gsl_ran_binomial_tpe(demeSize/2, m);
		numMaleMigrants = gsl_ran_binomial_tpe(demeSize/2, m);
		// n.b. migration is amongst the first numMigrants parents in each generation
		for(i = first0Female; i <= last0Male; i++)
		{
			if(ped->relationships[t][i][0] < first0Female + numFemaleMigrants)
				ped->relationships[t][i][0] += demeSize;
			if(ped->relationships[t][i][1] < first0Male + numMaleMigrants)
				ped->relationships[t][i][1] += demeSize;
		}
		for(i = first1Female; i <= last1Male; i++)
		{
			if(ped->relationships[t][i][0] < first1Female + numFemaleMigrants)
				ped->relationships[t][i][0] -= demeSize;
			if(ped->relationships[t][i][1] < first1Male + numMaleMigrants)
				ped->relationships[t][i][1] -= demeSize;
		}
	}
	return;
}

/* Note that this is reproduction followed by migration, forward in time. */
void CoalPedigree_shuffle_twopop_diploid_stochastic_rm(CoalPedigree * ped, double avgNumMigrants)
{
	int32_t t,i;
	const int32_t numGens = ped->numGenerations;
	const int32_t numFemales = ped->numFemales;
	const int32_t popnSize = ped->popnSize;
	const int32_t demeSize = popnSize/2;

	if(popnSize % 4 != 0)
	{
		fprintf(stderr, "Warning: Shuffling two-population pedigree with total population size that \
				is not a multiple of 4. (in CoalPedigree_shuffle_twopop)\n");
	}

	const int32_t first0Female = 0;
	const int32_t last0Female = demeSize/2-1;
	const int32_t first0Male = demeSize/2;
	const int32_t last0Male = demeSize-1;
	const int32_t first1Female = demeSize;
	const int32_t last1Female = demeSize + demeSize/2-1;
	const int32_t first1Male = demeSize + demeSize/2;
	const int32_t last1Male = popnSize-1;

	double m = avgNumMigrants/(double)(popnSize/2);

	// have to decide how to split up the migrants amongst males and females.
	int32_t numFemaleMigrants;
	int32_t numMaleMigrants;
	int32_t relationship;
	twoints migRelationships;


	/* For each generation, pick a father and a mother randomly for all sons and
	 * daughters */
	CoalPedigree_shuffle_twopop_nomigration(ped);
	for(t = 0; t < numGens; t++)
	{
		numFemaleMigrants = gsl_ran_binomial(demeSize/2, m);
		numMaleMigrants = gsl_ran_binomial(demeSize/2, m);
		// swap parent-relationships with an individual in the other deme.
		for(i = 0; i < numFemaleMigrants; i++)
		{
			migRelationships[0] = ped->relationships[t][first0Female+i][0];
			migRelationships[1] = ped->relationships[t][first0Female+i][1];
			ped->relationships[t][first0Female+i][0] = ped->relationships[t][first1Female+i][0];
			ped->relationships[t][first0Female+i][1] = ped->relationships[t][first1Female+i][1];
			ped->relationships[t][first1Female+i][0] = migRelationships[0];
			ped->relationships[t][first1Female+i][1] = migRelationships[1];
		}
		for(i = 0; i < numMaleMigrants; i++)
		{
			migRelationships[0] = ped->relationships[t][i+first0Male][0];
			migRelationships[1] = ped->relationships[t][i+first0Male][1];
			ped->relationships[t][i+first0Male][0] = ped->relationships[t][i+first1Male][0];
			ped->relationships[t][i+first0Male][1] = ped->relationships[t][i+first1Male][1];
			ped->relationships[t][i+first1Male][0] = migRelationships[0];
			ped->relationships[t][i+first1Male][1] = migRelationships[1];
		}
	}
	return;
}

/* To be removed */
void CoalPedigree_shuffle_twopop_diploid_stochastic_janet(CoalPedigree * ped, double avgNumMigrants)
{
	int32_t t,i;
	const int32_t numGens = ped->numGenerations;
	const int32_t numFemales = ped->numFemales;
	const int32_t popnSize = ped->popnSize;
	const int32_t demeSize = popnSize/2;
	if(popnSize % 4 != 0){
		fprintf(stderr, "Warning: Shuffling two-population pedigree with total population size that \
				is not a multiple of 4. (in CoalPedigree_shuffle_twopop)\n");
	}
	const int32_t first0Female = 0;
	const int32_t last0Female = demeSize/2-1;
	const int32_t first0Male = demeSize/2;
	const int32_t last0Male = demeSize-1;
	const int32_t first1Female = demeSize;
	const int32_t last1Female = demeSize + demeSize/2-1;
	const int32_t first1Male = demeSize + demeSize/2;
	const int32_t last1Male = popnSize-1;
	double m = avgNumMigrants/(double)(demeSize);
	int32_t numFemaleMigrants;
	int32_t numMaleMigrants;
	int32_t relationship;
	/* For each generation, pick a father and a mother randomly for all sons and
	 * daughters */
	CoalPedigree_shuffle_twopop_nomigration(ped);
	for(t = 0; t < numGens; t++)
	{
		numFemaleMigrants = gsl_ran_binomial_tpe(demeSize/2, m);
		numMaleMigrants = gsl_ran_binomial_tpe(demeSize/2, m);
		for(i = first0Female; i <= last0Male; i++)
		{
			if(ped->relationships[t][i][0] < first0Female + numFemaleMigrants)
				ped->relationships[t][i][0] += demeSize;
			if(ped->relationships[t][i][1] < first0Male + numMaleMigrants)
				ped->relationships[t][i][1] += demeSize;
		}
		for(i = first1Female; i <= last1Male; i++)
		{
			if(ped->relationships[t][i][0] < first1Female + numFemaleMigrants)
				ped->relationships[t][i][0] -= demeSize;
			if(ped->relationships[t][i][1] < first1Male + numMaleMigrants)
				ped->relationships[t][i][1] -= demeSize;
		}
	}
	return;
}

/* Creates a two-deme pedigree where in each `interval` generations there is
 * one migrant individual */
void CoalPedigree_shuffle_twopop_diploid_interval(CoalPedigree * ped, int32_t interval, int32_t firstMigGen)
{
	int32_t t,i;
	const int32_t numGens = ped->numGenerations;
	const int32_t numFemales = ped->numFemales;
	const int32_t popnSize = ped->popnSize;
	const int32_t demeSize = popnSize/2;
	if(popnSize % 4 != 0)
	{
		fprintf(stderr, "Warning: Shuffling two-population pedigree with total population size that \
				is not a multiple of 4. (in CoalPedigree_shuffle_twopop)\n");
	}
	int32_t femaleMigrant;
	int32_t relationship;
	const int32_t first0Female = 0;
	const int32_t last0Female = demeSize/2-1;
	const int32_t first0Male = demeSize/2;
	const int32_t last0Male = demeSize-1;
	const int32_t first1Female = demeSize;
	const int32_t last1Female = demeSize + demeSize/2-1;
	const int32_t first1Male = demeSize + demeSize/2;
	const int32_t last1Male = popnSize-1;
	CoalPedigree_shuffle_twopop_nomigration(ped);
	for(t = firstMigGen; t < numGens; t += interval)
	{
		femaleMigrant = rbern(0.5);
		if(femaleMigrant)
		{
			ped->relationships[t][first0Female][0] = runifd(first1Female, last1Female);
			ped->relationships[t][first0Female][1] = runifd(first1Male, last1Male);
			ped->relationships[t][first1Female][0] = runifd(first0Female, last0Female);
			ped->relationships[t][first1Female][1] = runifd(first0Male, last0Male);
		}
		else
		{
			ped->relationships[t][first0Male][0] = runifd(first1Female, last1Female);
			ped->relationships[t][first0Male][1] = runifd(first1Male, last1Male);
			ped->relationships[t][first1Male][0] = runifd(first0Female, last0Female);
			ped->relationships[t][first1Male][1] = runifd(first0Male, last0Male);
		}
	}
	return;
}

void CoalPedigree_recursive_remove_ancestry_(CoalPedigree * ped, int32_t indiv, int32_t curGen, int32_t ancGen, int32_t anc)
{
	const int32_t demeSize = ped->popnSize / ped->numDemes; 
	const int32_t deme = indiv / demeSize;
	const int32_t male = (indiv - demeSize*deme) / (demeSize/2);
	if(curGen == ancGen-1)
	{
		while(ped->relationships[curGen][indiv][0] == anc)
			ped->relationships[curGen][indiv][0] = runifd(deme*demeSize, deme*demeSize+demeSize/2-1);
		while(ped->relationships[curGen][indiv][1] == anc)
			ped->relationships[curGen][indiv][1] = runifd(deme*demeSize+demeSize/2, deme*(demeSize+1)-1);
		return;
	}
	CoalPedigree_recursive_remove_ancestry_(ped, ped->relationships[curGen][indiv][0], curGen+1, ancGen, anc);
	CoalPedigree_recursive_remove_ancestry_(ped, ped->relationships[curGen][indiv][1], curGen+1, ancGen, anc);
	return;
}

void CoalPedigree_remove_ancestry(CoalPedigree * ped, int32_t indiv, int32_t ancGen, int32_t anc)
{
	CoalPedigree_recursive_remove_ancestry_(ped, indiv, 0, ancGen, anc);
	return;
}


/* Shuffles a panmictic population and gives two individuals a shared ancestor
 * in generation relAncGen. Good to combine with CoalPedigree_remove_ancestry().
 * */
void CoalPedigree_shuffle_fixed_recent_relative(CoalPedigree * ped, int32_t * twoindivs, int32_t relAncGen)
{
	int32_t t,i,j, ancestor, ancestorP, ancestor1, sex;
	CoalPedigree_shuffle_dioecious(ped);
	// the 1-chromosome of twoindivs[0] is the inbred/migrant one.
	ancestor = ped->relationships[0][twoindivs[0]][1];		
	for(t = 1; t < relAncGen; t++)
	{
		ancestorP = ped->relationships[t][ancestor][rbern(0.5)];
		ancestor = ancestorP;
	}
	ancestor1 = twoindivs[1];
	for(t = 0; t < relAncGen-1; t++)
	{
		ancestorP = ped->relationships[t][ancestor1][rbern(0.5)];
		ancestor1 = ancestorP;
	}
	sex = (ancestor >= ped->popnSize/2);
	ped->relationships[relAncGen-1][ancestor1][sex] = ancestor;
	return;
}

/*
 * This function returns the index of the first migrant in generation
 * migAncGen. A migrant is defined as any individual with two parents
 * coming from the deme not containing in the individual. If there
 * are multiple such migrants, only the first (numerically) is 
 * returned. Assumes that ped->numDemes was set (and is > 1). May not
 * be compatible with shuffling functions that produce migrants who
 * have one migrant parent and one non-migrant parent. (e.g., haploid
 * shuffling functions.) Returns -1 when no migrants are found. */
int32_t CoalPedigree_find_migrant(CoalPedigree * ped, int32_t migAncGen, int32_t demeIdx)
{
	int32_t i;
	const int32_t demeSize = ped->popnSize/ped->numDemes;
	const int32_t firstDemeIndiv = demeIdx * demeSize;
	const int32_t lastDemeIndiv = (demeIdx+1)*demeSize - 1;
	for(i = firstDemeIndiv; i <= lastDemeIndiv; i++)
	{
		if( NOTININTERVAL(ped->relationships[migAncGen][i][0], firstDemeIndiv, lastDemeIndiv) &&
				NOTININTERVAL(ped->relationships[migAncGen][i][1], firstDemeIndiv, lastDemeIndiv))
			return(i);
	}
	return -1;
}

/* 
 * This function returns the weight/probability that an individual 'indiv' 
 * has in a single individual 'target' in generation migrantGen. wc and 
 * wcp must be popnSize-length double arrays. ancGen is 0-based.
 */
double CoalPedigree_get_anc_weight(CoalPedigree * ped, int32_t indiv, int32_t target, int32_t ancGen, double * wc, double * wcp)
{
	int32_t t, i;
	const int32_t popnSize = ped->popnSize;
	for(i = 0; i < popnSize; i++)
		wc[i] = wcp[i] = 0.0;
	wc[indiv] = 1.0;
	for(t = 0; t < ancGen; t++)		// note that this is < rather than <=.
	{									// (ped->relationships extends to the previous gen.)
		for(i = 0; i < popnSize; i++)
		{
			if(wc[i] > 0.0)
			{
				wcp[ped->relationships[t][i][0]] += 0.5*wc[i];
				wcp[ped->relationships[t][i][1]] += 0.5*wc[i];
			}
		}
		for(i = 0; i < popnSize; i++)
		{
			wc[i] = wcp[i];
			wcp[i] = 0.0;
		}
	}
	return wc[target];
}

/* 
 This function returns the weight/probability that an individual 'indiv' 
 has in a single migrant in generation migrantGen. There must be exactly one
 migrant in that generation in order for this to make sense. This also requires
 a two-deme pedigree. wc and wcp must be demeSize-length double arrays.

 Accomodates indivs in the first or second deme.

TODO: Determine whether this is the same as combining the above two functions
 (namely CoalPedigree_get_anc_weight() and CoalPedigree_find_migrant()).
 */
double CoalPedigree_get_mig_prob(CoalPedigree * ped, int32_t indiv, int32_t migrantGen, double * wc, double * wcp)
{
	int32_t t, i, migrantIdx, N = ped->popnSize/2;
	if(indiv < N)
	{
		for(i = 0; i < N; i++)
		{
			if(ped->relationships[migrantGen][i][0] >= N)
			{
				migrantIdx = i;
				break;
			}
		}
		for(i = 0; i < N; i++)
			wc[i] = wcp[i] = 0.0;
		wc[indiv] = 1.0;
		for(t = 0; t < migrantGen; t++)		// note that this is < rather than <=.
		{									// (ped->relationships extends to the previous gen.)
			for(i = 0; i < N; i++)
			{
				wcp[ped->relationships[t][i][0]] += 0.5*wc[i];
				wcp[ped->relationships[t][i][1]] += 0.5*wc[i];
			}
			for(i = 0; i < N; i++)
			{
				wc[i] = wcp[i];
				wcp[i] = 0.0;
			}
		}
		return wc[migrantIdx];
	}
	else
	{
		for(i = N; i < 2*N; i++)
		{
			if(ped->relationships[migrantGen][i][0] < N)
			{
				migrantIdx = i;
				break;
			}
		}
	
		for(i = 0; i < N; i++)
			wc[i] = wcp[i] = 0.0;
		wc[indiv-N] = 1.0;
		for(t = 0; t < migrantGen; t++)		// note that this is < rather than <=.
		{									// (ped->relationships extends to the previous gen.)
			for(i = N; i < 2*N; i++)
			{
				wcp[ped->relationships[t][i][0]-N] += 0.5*wc[i-N];
				wcp[ped->relationships[t][i][1]-N] += 0.5*wc[i-N];
			}
			for(i = 0; i < N; i++)
			{
				wc[i] = wcp[i];
				wcp[i] = 0.0;
			}
		}
		return wc[migrantIdx-N];
	}
}

void CoalPedigree_recursive_binary_tree_twopop_(CoalPedigree * ped, int32_t curGen, int32_t targetGen, int32_t curIndiv, int32_t parentIdx)
{
    if(curGen == targetGen)
        return;
	const int32_t demeSize = ped->popnSize/2;
	const int32_t deme = (curIndiv >= demeSize);
	const int32_t first0Female = 0;
	const int32_t last0Female = demeSize/2-1;
	const int32_t first0Male = demeSize/2;
	const int32_t last0Male = demeSize-1;
	const int32_t first1Female = demeSize;
	const int32_t last1Female = demeSize + demeSize/2-1;
	const int32_t first1Male = demeSize + demeSize/2;
	const int32_t last1Male = ped->popnSize-1;
    ped->relationships[curGen][curIndiv][0] = first1Female*deme + parentIdx;
    ped->relationships[curGen][curIndiv][1] = first1Female*deme + first0Male + parentIdx;
    CoalPedigree_recursive_binary_tree_twopop_(ped, curGen+1, targetGen, first1Female*deme + parentIdx, parentIdx*2);
    CoalPedigree_recursive_binary_tree_twopop_(ped, curGen+1, targetGen, first1Female*deme + first0Male+parentIdx, parentIdx*2+1);
    return;
}

/* Shuffles a pedigree such that each individual in indivs has a pedigree that
 * is a perfect binary tree with no migration, except for one individual (the
 * first in indivs), who has a parent from the other deme migAncGen generations
 * ago.
 *
 * ped         initialized two-deme pedigree
 * indivs      vector of indices of sampled individuals. first individual will have migrant ancestry, others will have no migrant ancestry
 * numIndivs   total number of individuals
 * migAncGen   number of generations ago in which the migrant ancestor is found (migrant in the sense that they are in the other deme)
 * clearGen    generation until which each individual's ancestry is a perfect binary tree (inclusive)
 *
 */

void CoalPedigree_shuffle_twopop_fixed_migrant(CoalPedigree * ped, double avgNumMigrants, int32_t * indivs, int32_t numIndivs, int32_t migAncGen, int32_t clearGen)
{
    // background shuffling
    CoalPedigree_shuffle_twopop_diploid_fully_stochastic(ped, avgNumMigrants);

    // perfect binary trees
    int32_t count0 = 0, count1 = 0, idx, i, demeSize = ped->popnSize/2;
    for(i = 0; i < numIndivs; i++)
    {
        if(indivs[i] < demeSize)
            idx = count0++;
        else
            idx = count1++;
        CoalPedigree_recursive_binary_tree_twopop_(ped, 0, clearGen, indivs[i], idx);
    }

    // add the migrant relationship for the first individual
    int32_t curDeme = (indivs[0] >= demeSize);
    int32_t first0Female = 0;
    int32_t last0Female = demeSize/2-1;
    int32_t first1Female = demeSize;
    int32_t last1Female = demeSize+demeSize/2-1;
    int32_t firstAvailableFemaleAnc, newMaternalAncestor;
    int32_t t, curIndiv = indivs[0];
    for(t = 0; t < migAncGen-1; t++)
        curIndiv = ped->relationships[t][curIndiv][0];
    if(curDeme == 0)
    {
        firstAvailableFemaleAnc = first1Female + count1*pow(2,migAncGen);
        if(firstAvailableFemaleAnc > last1Female)
            PERROR("Too small deme size in CoalPedigree_shuffle_twopop_fixed_migrant()");
        newMaternalAncestor = runifd(firstAvailableFemaleAnc, last1Female);
        ped->relationships[migAncGen-1][curIndiv][0] = newMaternalAncestor;
    }
    else // curDeme == 1
    {
        firstAvailableFemaleAnc = first0Female + count0*pow(2,migAncGen);
        if(firstAvailableFemaleAnc > last0Female)
            PERROR("Too small deme size in CoalPedigree_shuffle_twopop_fixed_migrant()");
        newMaternalAncestor = runifd(firstAvailableFemaleAnc, last0Female);
        ped->relationships[migAncGen-1][curIndiv][0] = newMaternalAncestor;
    }
    return;
}

/* 
 This function returns the weight/probability that an individual 'indiv' 
 has in a single migrant in generation migrantGen. There must be exactly one
 migrant in that generation in order for this to make sense. This also requires
 a two-deme pedigree. wc and wcp must be demeSize-length double arrays.

 Accomodates indivs in the first or second deme.
 */
double CoalPedigree_shuffle_twopop_fixed_weight_single_migrant(CoalPedigree * ped, int32_t * indivs, int32_t numIndivsPerDeme, int32_t migAncGen, int32_t interval)
{
	int32_t t,i,j, ancestor, ancestorP;
	const int32_t indiv = indivs[0];
	const int32_t popnSize = ped->popnSize;
	const int32_t demeSize = popnSize/2;
	const int32_t pop = (indiv >= demeSize);
	const int32_t resumeGen = 2*(int32_t)(log2((double)demeSize));
	const int32_t first0Female = 0;
	const int32_t last0Female = demeSize/2-1;
	const int32_t first0Male = demeSize/2;
	const int32_t last0Male = demeSize-1;
	const int32_t first1Female = demeSize;
	const int32_t last1Female = demeSize + demeSize/2-1;
	const int32_t first1Male = demeSize + demeSize/2;
	const int32_t last1Male = popnSize-1;
	double * wc, * wcp;
	double migWeight;
	CoalPedigree_shuffle_twopop_diploid_interval(ped, interval, resumeGen);
	// Previously the migrant was on a random lineage. Now it's always paternal.
	//ancestor = ped->relationships[0][indiv][rbern(0.5)];
	ancestor = ped->relationships[0][indiv][1];
	for(t = 1; t < migAncGen; t++)
	{
		// this is still random, though, no need to make it always paternal.
		ancestorP = ped->relationships[t][ancestor][rbern(0.5)];
		ancestor = ancestorP;
	}
	// note that there is no reciprocal migration from the other deme into pop.
	if(pop == 0)
	{
		ped->relationships[migAncGen][ancestor][0] = runifd(first1Female, last1Female);
		ped->relationships[migAncGen][ancestor][1] = runifd(first1Male, last1Male);
	}
	else
	{
		ped->relationships[migAncGen][ancestor][0] = runifd(first0Female, last0Female);
		ped->relationships[migAncGen][ancestor][1] = runifd(first0Male, last0Male);
	}
	// at this point, ancestor is the migrant ancestor in generation migAncGen.
	// for each of the remaining individuals in this deme in indivs, remove migrant ancestry in those individuals.
	wc = (double *)malloc(sizeof(double) * (size_t)demeSize);
	wcp = (double *)malloc(sizeof(double) * (size_t)demeSize);
	CHECKPOINTER(wc);
	CHECKPOINTER(wcp);
	for(i = 1; i < numIndivsPerDeme; i++)
		CoalPedigree_remove_ancestry(ped, indivs[i], migAncGen, ancestor);
	migWeight = CoalPedigree_get_mig_prob(ped, indiv, migAncGen, wc, wcp);
	free(wc);
	free(wcp);
	return migWeight;
}

/* The following function shuffles a twopop CoalPedigree and gives two
 * individuals 'indivs' a shared ancestor in generation 'relAncGen'.
 * Migration starts at generation 2*log(demeSize). *This function 
 * assumes that both individuals come from deme 0.* */
void CoalPedigree_shuffle_twopop_fixed_recent_relative(CoalPedigree * ped, int32_t * indivs, int32_t numIndivsPerDeme, int32_t relAncGen, int32_t interval)
{
	int32_t t,i,j, ancestor, ancestorP, ancestor1, ancestorP1, sex;
	const int32_t demeSize = ped->popnSize/2;
	const int32_t resumeGen = 2*(int32_t)(log2((double)demeSize));
	CoalPedigree_shuffle_twopop_diploid_interval(ped, interval, resumeGen);
	ancestor = ped->relationships[0][indivs[0]][1];
	for(t = 1; t < relAncGen; t++)
	{
		ancestorP = ped->relationships[t][ancestor][rbern(0.5)];
		ancestor = ancestorP;
	}
	ancestor1 = indivs[1];
	for(t = 0; t < relAncGen-1; t++)
	{
		ancestorP1 = ped->relationships[t][ancestor1][rbern(0.5)];
		ancestor1 = ancestorP1;
	}
	sex = (ancestor >= demeSize/2);
	ped->relationships[relAncGen-1][ancestor1][sex] = ancestor;
	for(i = 2; i < numIndivsPerDeme; i++)
		CoalPedigree_remove_ancestry(ped, indivs[i], relAncGen, ancestor);
	return;
}

/* numMigrants is now the number of offspring in each generation that have one
 * migrant parent. haploid refers to haploid migration */
void CoalPedigree_shuffle_twopop_haploid_fixed(CoalPedigree * ped, int32_t numMigrants)
{
	int32_t t,i;
	const int32_t numGens = ped->numGenerations;
	const int32_t popnSize = ped->popnSize;
	const int32_t demeSize = popnSize/2;

	if(popnSize % 4 != 0)
	{
		fprintf(stderr, "Warning: Shuffling two-population pedigree with total population size that \
				is not a multiple of 4. (in CoalPedigree_shuffle_twopop)\n");
	}

	const int32_t first0Female = 0;
	const int32_t last0Female = demeSize/2-1;
	const int32_t first0Male = demeSize/2;
	const int32_t last0Male = demeSize-1;
	const int32_t first1Female = demeSize;
	const int32_t last1Female = demeSize + demeSize/2-1;
	const int32_t first1Male = demeSize + demeSize/2;
	const int32_t last1Male = popnSize-1;
	int32_t relationship;

	/* haploid shuffling assumes that the migrant gamete is always male (pollen). 
	 * this means that num[Fem|M]aleMigrants is the number of offspring
	 * sired by migrant pollen/gametes. */
	int32_t numFemaleMigrants = numMigrants/2;
	int32_t numMaleMigrants = numMigrants - numFemaleMigrants;

	CoalPedigree_shuffle_twopop_nomigration(ped);
	for(t = 0; t < numGens; t++)
	{
		// relationships for migrant females in pop 0.
		for(i = first0Female; i < first0Female + numFemaleMigrants; i++)
		{
			// runifd is inclusive
			relationship = runifd(first1Male, last1Male);
			ped->relationships[t][i][1] = relationship;
		}
		// relationships for migrant males in pop 0.
		for(i = first0Male; i < first0Male + numMaleMigrants; i++)
		{
			relationship = runifd(first1Male, last1Male);
			ped->relationships[t][i][1] = relationship;
		}
		// relationships for migrant females in pop 1
		for(i = first1Female; i < first1Female + numFemaleMigrants; i++)
		{
			relationship = runifd(first0Male, last0Male);
			ped->relationships[t][i][1] = relationship;
		}
		// relationships for migrant males in pop 1.
		for(i = first1Male; i < first1Male + numMaleMigrants; i++)
		{
			relationship = runifd(first0Male, last0Male);
			ped->relationships[t][i][1] = relationship;
		}
	}
	return;
}

void CoalPedigree_shuffle_twopop_haploid_interval(CoalPedigree * ped, int32_t interval)
{
	int32_t t,i;
	const int32_t numGens = ped->numGenerations;
	const int32_t numFemales = ped->numFemales;
	const int32_t popnSize = ped->popnSize;
	const int32_t demeSize = popnSize/2;

	if(popnSize % 4 != 0)
	{
		fprintf(stderr, "Warning: Shuffling two-population pedigree with total population size that \
				is not a multiple of 4. (in CoalPedigree_shuffle_twopop)\n");
	}

	const int32_t first0Female = 0;
	const int32_t last0Female = demeSize/2-1;
	const int32_t first0Male = demeSize/2;
	const int32_t last0Male = demeSize-1;
	const int32_t first1Female = demeSize;
	const int32_t last1Female = demeSize + demeSize/2-1;
	const int32_t first1Male = demeSize + demeSize/2;
	const int32_t last1Male = popnSize-1;

	int32_t femaleMigrant;

	CoalPedigree_shuffle_twopop_nomigration(ped);

	for(t = 0; t < numGens; t += interval)
	{
		// assuming it's a male migrant parent (see ..._haploid_fixed)
		// (a random 0 individual has a random 1 father)
		ped->relationships[t][runifd(first0Female,last0Male)][1] = runifd(first1Male,last1Male);
		// (a random 1 individual has a random 0 father)
		ped->relationships[t][runifd(first1Female,last1Male)][1] = runifd(first0Male,last0Male);
	}
	return;
}

void CoalPedigree_shuffle_twopop_haploid_stochastic(CoalPedigree * ped, double avgNumMigrants)
{
	int32_t t,i;
	const int32_t numGens = ped->numGenerations;
	const int32_t popnSize = ped->popnSize;
	const int32_t demeSize = popnSize/2;

	if(popnSize % 4 != 0)
	{
		fprintf(stderr, "Warning: Shuffling two-population pedigree with total population size that \
				is not a multiple of 4. (in CoalPedigree_shuffle_twopop)\n");
	}
	const int32_t first0Female = 0;
	const int32_t last0Female = demeSize/2-1;
	const int32_t first0Male = demeSize/2;
	const int32_t last0Male = demeSize-1;
	const int32_t first1Female = demeSize;
	const int32_t last1Female = demeSize + demeSize/2-1;
	const int32_t first1Male = demeSize + demeSize/2;
	const int32_t last1Male = popnSize-1;

	int32_t numFemaleMigrants0;
	int32_t numMaleMigrants0;
	int32_t numFemaleMigrants1;
	int32_t numMaleMigrants1;

	double m = avgNumMigrants/((double)demeSize);

	int32_t relationship;

	/* For each generation, pick a father and a mother randomly for all sons and
	 * daughters */
	CoalPedigree_shuffle_twopop_nomigration(ped);
	for(t = 0; t < numGens; t++)
	{
		numFemaleMigrants0 = gsl_ran_binomial(demeSize/2, m);
		numMaleMigrants0 = gsl_ran_binomial(demeSize/2, m);
		numFemaleMigrants1 = gsl_ran_binomial(demeSize/2, m);
		numMaleMigrants1 = gsl_ran_binomial(demeSize/2, m);
		/* Pick an integer randomly numF times from 0, ... , numM, and fill the
		 * given array with those numbers. */
		//for(i = 0; i <= last0Indiv; i++)
		// relationships for migrant females in pop 0
		for(i = first0Female; i < first0Female + numFemaleMigrants0; i++)
		{
			// runifd is inclusive
			relationship = runifd(first0Female, last0Female);
			ped->relationships[t][i][0] = relationship;

			relationship = runifd(first1Male, last1Male);
			ped->relationships[t][i][1] = relationship;
		}
		// relationships for migrant males in pop 0.
		for(i = first0Male; i < first0Male + numMaleMigrants0; i++)
		{
			// runifd is inclusive
			relationship = runifd(0, last0Female);
			ped->relationships[t][i][0] = relationship;

			relationship = runifd(first1Male, last1Male);
			ped->relationships[t][i][1] = relationship;
		}
		// relationships for migrant females in pop 1
		for(i = first1Female; i < first1Female + numFemaleMigrants1; i++)
		{
			// runifd is inclusive
			relationship = runifd(first1Female, last1Female);
			ped->relationships[t][i][0] = relationship;

			relationship = runifd(first0Male, last0Male);
			ped->relationships[t][i][1] = relationship;
		}
		// relationships for migrant males in pop 1.
		for(i = first1Male; i < first1Male + numMaleMigrants1; i++)
		{
			// runifd is inclusive
			relationship = runifd(first1Female, last1Female);
			ped->relationships[t][i][0] = relationship;

			relationship = runifd(first0Male, last0Male);
			ped->relationships[t][i][1] = relationship;
		}
	}
	return;
}

/*  loopGen		generation looped back to after going through pedigree
 *
 */
void CoalPedigree_shuffle_dioecious_varN(CoalPedigree * ped, int32_t loopGen)
{
	int32_t t,i;
	const int32_t numGens = ped->numGenerations;
	const int32_t numFemales = ped->numFemales;
	const int32_t popnSize = ped->popnSize;
	const int32_t * popnSizes = ped->popnSizes;
	int32_t currentN, parentN, parentNumFemales, parentNumMales;

	if(!(ped->variablePopnSize))
		PERROR("Attempted shuffling of non-varN pedigree in CoalPedigree_shuffle_dioecious_varN");

	// give parents of final generation the same size as loopGen
	for(i = 0; i < popnSizes[numGens-1]; i++)
	{
		ped->relationships[numGens-1][i][0] = runifd(0,popnSizes[loopGen]/2-1);
		ped->relationships[numGens-1][i][1] = runifd(popnSizes[loopGen]/2,popnSizes[loopGen]-1);
	}

	for(t = 0; t < numGens-1; t++)
	{
		currentN = popnSizes[t];
		parentN = popnSizes[t+1];
		parentNumFemales = parentN / 2;
		parentNumMales = parentN / 2;
		/* Pick an integer randomly numF times from 0, ... , numM, and fill the
		 * given array with those numbers. */
		for(i = 0; i < currentN; i++)
		{
			// runifd is inclusive
			ped->relationships[t][i][0] = runifd(0,parentNumFemales-1);
			ped->relationships[t][i][1] = runifd(parentNumFemales,parentN-1);
		}
	}
	/*
	// Notice how the loops stops at numGens-2 above. That's because we don't provide the size of the
	// last (back in time) generation's parents. The loop below is just to guard against bugs.
	for(i = 0; i < currentN; i++)
	{
		ped->relationships[numGens-1][i][0] = runifd(0,parentNumFemales-1);
		ped->relationships[numGens-1][i][1] = runifd(parentNumFemales,parentN-1);
	}
	*/
	return;
}

/* shuffles a monoecious pedigree
 *
 * ped	pointer to initialized pedigree
 */
void CoalPedigree_shuffle_monoecious(CoalPedigree * ped)
{
	int32_t t,i;
	/* Not sure whether numFemales is ever used in the coal. time calculations. */
	const int32_t numGens = ped->numGenerations;
	const int32_t popnSize = ped->popnSize;

	/* For each generation, pick a father and a mother randomly for all sons and
	 * daughters */
	for(t = 0; t < numGens; t++)
	{
		for(i = 0; i < popnSize; i++)
		{
			ped->relationships[t][i][0] = runifd(0,popnSize-1);
			ped->relationships[t][i][1] = runifd(0,popnSize-1);
		}
	}
	return;
}

/* shuffles a "cyclical" dioecious pedigree, where the same relationships are
 * repeated in every generation
 *
 * ped	pointer to initialized pedigree
 */
void CoalPedigree_shuffle_cyc_WF_dio(CoalPedigree * ped)
{
	int32_t t,i;
	const int32_t numFemales = ped->numFemales;
	const int32_t numGens = ped->numGenerations;
	const int32_t popnSize = ped->popnSize;

	for(i = 0; i < popnSize; i++)
	{
		ped->relationships[0][i][0] = runifd(0,numFemales-1);
		ped->relationships[0][i][1] = runifd(numFemales,popnSize-1);
	}
	for(t = 1; t < numGens; t++)
	{
		for(i = 0; i < popnSize; i++)
		{
			ped->relationships[t][i][0] = ped->relationships[0][i][0];
			ped->relationships[t][i][1] = ped->relationships[0][i][1];
		}
	}
	return;
}

/* shuffles a "cyclical" monoecious pedigree, where the same relationships are
 * repeated in every generation
 *
 * ped	pointer to initialized pedigree
 */
void CoalPedigree_shuffle_cyc_WF_monoecious(CoalPedigree * ped)
{
	int32_t t,i;
	const int32_t numGens = ped->numGenerations;
	const int32_t popnSize = ped->popnSize;

	for(i = 0; i < popnSize; i++)
	{
		ped->relationships[0][i][0] = runifd(0,popnSize-1);
		ped->relationships[0][i][1] = runifd(0,popnSize-1);
	}
	for(t = 1; t < numGens; t++)
	{
		ped->relationships[t][i][0] = ped->relationships[t-1][i][0];
		ped->relationships[t][i][1] = ped->relationships[t-1][i][1];
	}
	return;
}

/* (used in CoalPedigree_shuffle_sweep_dio()) */
inline void CoalPedigree_update_type_probs_(double * tprobs, double * cTprobs, int32_t * numGenotypes, int32_t * gPerm, double s, double h)
{
	int32_t i, itemp, swapped;
	double dtemp;
	double total = 0.0;
	cTprobs[0] = 0.0;
	//double max[2] = {0.0,0.0};

	tprobs[0] = (double)numGenotypes[0];
	total += tprobs[0];
	tprobs[1] = ((double)numGenotypes[1] * (1+s*h));
	total += tprobs[1];
	tprobs[2] = ((double)numGenotypes[2] * (1+s));
	total += tprobs[2];
	
	tprobs[0] /= total;
	tprobs[1] /= total;
	tprobs[2] /= total;

	for(i = 0; i < 3; i++)
	{
		gPerm[i] = i;
		cTprobs[i] = tprobs[i];
	}

	// Bubble sort the three cTprobs into descending order.
	while(1)
	{
		swapped = 0;
		for(i = 0; i < 2; i++)
		{
			// < because we want to sort from highest to lowest.
			if(cTprobs[i] < cTprobs[i+1])
			{
				dtemp = cTprobs[i+1];
				cTprobs[i+1] = cTprobs[i];
				cTprobs[i] = dtemp;
				itemp = gPerm[i];
				gPerm[i] = gPerm[i+1];
				gPerm[i+1] = itemp;
				swapped = 1;
			}
		}
		if(!swapped)
			break;
	}

	// Now, make them cumulative.
	cTprobs[1] += cTprobs[0];
	cTprobs[2] += cTprobs[1];
	return;
}

/* Simulates pedigrees with selective sweep in history.
 *
 * ped            pedigree
 * ts             number of generations ago that sweep started
 * s              selective advantage of beneficial allele (fitness of heterozygous is 1+s*h, homozygote is 1+s)
 * h              dominance coefficient                                          0 <= h <= 1
 * requireSweep   indicator for whether to require fixation before the present
 *
 */

void CoalPedigree_shuffle_sweep_dio(CoalPedigree * ped, int32_t ts, double s, double h, int32_t requireSweep)
{
	if(ts > ped->numGenerations)
		PERROR("Selective sweep to start farther back in time than CoalPedigree extends.");
	if(ts < 0)
		PERROR("Negative ts value in CoalPedigree_shuffle_sweep_dio\nts = %i\n",ts);
	int32_t t,i,j;
	const int32_t numGens = ped->numGenerations;
	const int32_t numFemales = ped->numFemales;
	const int32_t numMales = ped->numMales;
	const int32_t popnSize = ped->popnSize;

	// numGenotypes[i][j] is the number of individuals of sex i and 
	// genotype j in the generation currently under consideration.
	int32_t numGenotypes[2][3] = {{numFemales,0,0},{numMales,0,0}};
	// genotypePerm[i][j] is the jth most common genotype amongst
	// individuals of sex i in the generation currently under consideration.
	// It's updated each generation by CoalPedigree_update_type_probs_.
	int32_t genotypePerm[2][3];
	// typeProbs[i][j] gives the probability that an individual chooses a parent
	// (of sex i) of genotype genotypePerm[i][j]
	double typeProbs[2][3]; 
	// cTypeProbs is the cumulative version of typeProbs.
	double cTypeProbs[2][3];
	// initialPos is the initial individual in which the beneficial mutation
	// occurs in generation ts.
	int32_t initialPos;
	double rU;
	int fixed[2] = {0,0};

	// parentGenotypes[i][j] gives the genotype of parent-of-sex-j of
	// individual i in the generation currently under consideration.
	//fprintf(stderr,"popnSize = %i\n",popnSize);
	twoints * parentGenotypes = (twoints *)malloc(sizeof(twoints) * popnSize);
	CHECKPOINTER(parentGenotypes);

	//Shuffle neutrally until the mutation arises.
	//for(t = numGens-1; t > ts; t--)
	for(t = numGens-1; t >= ts; t--)
	{
		for(i = 0; i < popnSize; i++)
		{
			ped->relationships[t][i][0] = runifd(0,numFemales-1);
			ped->relationships[t][i][1] = runifd(numFemales,popnSize-1);
			//ped->genotypes[t][i] = 0;
		}
	}

	// genotypeIdxs[s][g][i] is the index (\in {0,...,popnSize-1}) of the ith
	// individual of genotype g (\in {0,1,2}) amongst individuals of sex s.
	static int32_t *** genotypeIdxs; // Note that resizing popnSize in the middle of
									 // the running of the program might cause a segfault.
	if(genotypeIdxs == NULL)
	{
		genotypeIdxs = (int32_t ***)malloc(sizeof(int32_t **) * 2);
		CHECKPOINTER(genotypeIdxs);
		for(i = 0; i < 2; i++)
		{
			genotypeIdxs[i] = (int32_t **)malloc(sizeof(int32_t *) * 3);
			CHECKPOINTER(genotypeIdxs[i]);
			for(j = 0; j < 3; j++)
			{
				genotypeIdxs[i][j] = (int32_t *)malloc(sizeof(int32_t) * popnSize);
				CHECKPOINTER(genotypeIdxs[i][j]);
			}
		}
	}
	
	// Deciding where the mutation initially arises...
	initialPos = runifd(0,popnSize-1);
	if(initialPos < numFemales) // it arises in a female...
	{
		numGenotypes[0][1]++;
		numGenotypes[0][0]--;
	}
	else // it arises in a male...
	{
		(numGenotypes[1][1])++;
		(numGenotypes[1][0])--;
	}


	int curPos = 0;
	for(i = 0; i < numFemales; i++)
	{
		if(initialPos != i)
		{
			genotypeIdxs[0][0][curPos++] = i;
			//ped->genotypes[ts][i] = 0;
		}
		else
		{
			genotypeIdxs[0][1][0] = i;
			//ped->genotypes[ts][i] = 1;
		}
	}

	curPos = 0;
	for(i = numFemales; i < popnSize; i++)
	{
		if(initialPos != i)
		{
			genotypeIdxs[1][0][curPos++] = i;
			//ped->genotypes[ts][i] = 0;
		}
		else
		{
			genotypeIdxs[1][1][0] = i;
			//ped->genotypes[ts][i] = 1;
		}
	}


	int curGenotype;
	int32_t * curGenotypes = (int32_t *)malloc(sizeof(int32_t) * popnSize);
	CHECKPOINTER(curGenotypes);

	
	int finalGenotype;
	for(t = ts-1; t >= 0; t--) 
	{						// genotypes and produces the next [(t-1)'s] generation's genotypes.
		CoalPedigree_update_type_probs_(typeProbs[0],cTypeProbs[0],numGenotypes[0],genotypePerm[0],s,h);
		CoalPedigree_update_type_probs_(typeProbs[1],cTypeProbs[1],numGenotypes[1],genotypePerm[1],s,h);
		for(i = 0; i < popnSize; i++)
		{
			rU = runif();
			for(j = 0; j < 3; j++)
			{
				if(rU <= cTypeProbs[0][j])
				{
					parentGenotypes[i][0] = genotypePerm[0][j];
					break;
				}
			}
			rU = runif();
			for(j = 0; j < 3; j++)
			{
				if(rU <= cTypeProbs[1][j])
				{
					parentGenotypes[i][1] = genotypePerm[1][j];
					break;
				}
			}
			ped->relationships[t][i][0] = genotypeIdxs[0][parentGenotypes[i][0]][runifd(0,numGenotypes[0][parentGenotypes[i][0]]-1)];
			ped->relationships[t][i][1] = genotypeIdxs[1][parentGenotypes[i][1]][runifd(0,numGenotypes[1][parentGenotypes[i][1]]-1)];

			curGenotypes[i] = 0;
			for(j = 0; j < 2; j++)
			{
				switch(parentGenotypes[i][j])
				{
					case 0:
						break;
					case 2:
						curGenotypes[i] += 1;
						break;
					case 1:
						curGenotypes[i] += rbern(0.5);
						break;
					default:
						PERROR("Bad curGenotypes[i]");
				}
			}
			//ped->genotypes[t][i] = curGenotypes[i];
		}

		// Now, update the genotypes.
		for(i = 0; i < 3; i++)
		{
			numGenotypes[0][i] = 0;
			numGenotypes[1][i] = 0;
		}

		for(i = 0; i < numFemales; i++)
		{
			curGenotype = curGenotypes[i];
			genotypeIdxs[0][curGenotype][numGenotypes[0][curGenotype]] = i;
			(numGenotypes[0][curGenotype])++;
		}
		for(i = numFemales; i < popnSize; i++)
		{
			curGenotype = curGenotypes[i];
			genotypeIdxs[1][curGenotype][numGenotypes[1][curGenotype]] = i;
			(numGenotypes[1][curGenotype])++;
		}

		if((fixed[1] = (numGenotypes[0][2] == numFemales && numGenotypes[1][2] == popnSize-numFemales)) || 
				(fixed[0] = (numGenotypes[0][0] == numFemales && numGenotypes[1][0] == popnSize-numFemales)))
		{
			if(requireSweep && fixed[0])
			{
				//fprintf(stderr,"Lost allele.\n");
				numGenotypes[0][0] = numFemales;
				numGenotypes[0][1] = 0;
				numGenotypes[0][2] = 0;

				numGenotypes[1][0] = popnSize - numFemales;
				numGenotypes[1][1] = 0;
				numGenotypes[1][2] = 0;
				if(initialPos < numFemales) // it arises in a female...
				{
					numGenotypes[0][1]++;
					numGenotypes[0][0]--;
				}
				else // it arises in a male...
				{
					(numGenotypes[1][1])++;
					(numGenotypes[1][0])--;
				}
				curPos = 0;
				for(i = 0; i < numFemales; i++)
				{
					if(initialPos != i)
					{
						genotypeIdxs[0][0][curPos++] = i;
						//ped->genotypes[ts][i] = 0; // Should be unnecessary.
					}
					else
					{
						genotypeIdxs[0][1][0] = i;
						//ped->genotypes[ts][i] = 1;
					}
				}

				curPos = 0;
				for(i = numFemales; i < popnSize; i++)
				{
					if(initialPos != i)
					{
						genotypeIdxs[1][0][curPos++] = i;
						//ped->genotypes[ts][i] = 0;
					}
					else
					{
						genotypeIdxs[1][1][0] = i;
						//ped->genotypes[ts][i] = 1;
					}
				}
				t = ts; // To be deincremented by the loop to ts-1 to start over.
				//printf("t = ts... starting over.\n");
			}
			else // If it fixed at 1 or it fixed and we don't care whether it fixed at 1 or 0...
			{
				int tp;
				finalGenotype = 2*fixed[1];
				for(tp = t-1; tp >= 0; tp--)
				{
					for(i = 0; i < popnSize; i++)
					{
						ped->relationships[tp][i][0] = runifd(0,numFemales-1);
						ped->relationships[tp][i][1] = runifd(numFemales,popnSize-1);
						//ped->genotypes[tp][i] = finalGenotype;	
					}
				}
				break; // breaks from the t loop.
			}
		}
	}
	if(!fixed[0] && !fixed[1])
		fprintf(stderr, "WARNING: selective sweep did not fix.\n");

	free(parentGenotypes);	
	free(curGenotypes);	

	return;
}

/* Same as above except with a variable population size. */
void CoalPedigree_shuffle_sweep_dio_varN(CoalPedigree * ped, int32_t ts, double s, double h, int32_t requireSweep)
{
	if(ts > ped->numGenerations)
		PERROR("Selective sweep to start farther back in time than CoalPedigree extends.");
	if(ts < 0)
		PERROR("Negative ts value in CoalPedigree_shuffle_sweep_dio\nts = %i\n",ts);
	int32_t t,i,j;
	const int32_t numGens = ped->numGenerations;
	const int32_t * popnSizes = ped->popnSizes;
	int32_t maxPopnSize = 0;
	for(i = 0; i < numGens; i++)
	{
		if(maxPopnSize < popnSizes[i])
			maxPopnSize = popnSizes[i];
	}
	int32_t currentPopnSize = popnSizes[0];
	int32_t parentPopnSize;

	// numGenotypes[i][j] is the number of individuals of sex i and 
	// genotype j in the generation currently under consideration.
	//int32_t numGenotypes[2][3] = {{numFemales,0,0},{numMales,0,0}};
	//int32_t numGenotypes[2][3] = {{popnSizes[0]/2,0,0},{popnSizes[0]/2,0,0}};
	int32_t numGenotypes[2][3];
	// genotypePerm[i][j] is the jth most common genotype amongst
	// individuals of sex i in the generation currently under consideration.
	// It's updated each generation by CoalPedigree_update_type_probs_.
	int32_t genotypePerm[2][3];
	// typeProbs[i][j] gives the probability that an individual chooses a parent
	// (of sex i) of genotype genotypePerm[i][j]
	double typeProbs[2][3]; 
	// cTypeProbs is the cumulative version of typeProbs.
	double cTypeProbs[2][3];
	// initialPos is the initial individual in which the beneficial mutation
	// occurs in generation ts.
	int32_t initialPos;
	double rU;
	int fixed[2] = {0,0};

	// parentGenotypes[i][j] gives the genotype of parent-of-sex-j of
	// individual i in the generation currently under consideration.
	twoints * parentGenotypes = (twoints *)malloc(sizeof(twoints) * maxPopnSize);
	CHECKPOINTER(parentGenotypes);

	//Shuffle neutrally until the mutation arises.

	// We don't have a population size for the first generation's parents, so assume it's the same.
	// Hopefully prevents some of the segfaults that were popping up.
	parentPopnSize = currentPopnSize = popnSizes[numGens-1];
	for(i = 0; i < currentPopnSize; i++)
	{
		ped->relationships[numGens-1][i][0] = runifd(0,parentPopnSize/2-1);
		ped->relationships[numGens-1][i][1] = runifd(parentPopnSize/2,parentPopnSize-1);
	}

	for(t = numGens-2; t >= ts; t--)
	{
		parentPopnSize = popnSizes[t+1];
		currentPopnSize = popnSizes[t];
		for(i = 0; i < currentPopnSize; i++)
		{
			ped->relationships[t][i][0] = runifd(0,parentPopnSize/2-1);
			ped->relationships[t][i][1] = runifd(parentPopnSize/2,parentPopnSize-1);
		}
	}

	// genotypeIdxs[s][g][i] is the index (\in {0,...,popnSize-1}) of the ith
	// individual of genotype g (\in {0,1,2}) amongst individuals of sex s.
	static int32_t *** genotypeIdxs; // Note that resizing popnSize in the middle of
									 // the running of the program might cause a segfault.
	if(genotypeIdxs == NULL)
	{
		genotypeIdxs = (int32_t ***)malloc(sizeof(int32_t **) * 2);
		CHECKPOINTER(genotypeIdxs);
		for(i = 0; i < 2; i++)
		{
			genotypeIdxs[i] = (int32_t **)malloc(sizeof(int32_t *) * 3);
			CHECKPOINTER(genotypeIdxs[i]);
			for(j = 0; j < 3; j++)
			{
				genotypeIdxs[i][j] = (int32_t *)malloc(sizeof(int32_t) * maxPopnSize);
				CHECKPOINTER(genotypeIdxs[i][j]);
			}
		}
	}
	
	currentPopnSize = popnSizes[ts]; 		// Working.
	//printf("currentPopnSize = %i (when the mutation arises)\nt = %i\nts = %i\n",currentPopnSize,t,ts);

 	numGenotypes[0][0] = numGenotypes[1][0] = currentPopnSize/2;
	numGenotypes[0][1] = numGenotypes[0][2] = numGenotypes[1][1] = numGenotypes[1][2] = 0;
 	
	// Deciding where the mutation initially arises...
	initialPos = runifd(0,currentPopnSize-1);
	if(initialPos < currentPopnSize/2) // it arises in a female...
	{
		numGenotypes[0][1]++;
		numGenotypes[0][0]--;
	}
	else // it arises in a male...
	{
		(numGenotypes[1][1])++;
		(numGenotypes[1][0])--;
	}


	int curPos = 0;
	for(i = 0; i < currentPopnSize/2; i++)
	{
		if(initialPos != i)
			genotypeIdxs[0][0][curPos++] = i;
		else
			genotypeIdxs[0][1][0] = i;
	}

	curPos = 0;
	for(i = currentPopnSize/2; i < currentPopnSize; i++)
	{
		if(initialPos != i)
			genotypeIdxs[1][0][curPos++] = i;
		else
			genotypeIdxs[1][1][0] = i;
	}

	int curGenotype;
	int32_t * curGenotypes = (int32_t *)malloc(sizeof(int32_t) * maxPopnSize);
	CHECKPOINTER(curGenotypes);

	
	double freq;
	int32_t numSelectedAlleles;
	int32_t finalGenotype;
	const int32_t printFrequencies = 0;
	const int32_t printLastPedigree = 0;


	/*
	int32_t * numsOutOfRange = (int32_t *)malloc(sizeof(int32_t) * (size_t)ts);
	int32_t numCompleted = 0;
	*/

	/*
	FILE * freqs;
	if(printFrequencies)
	{
		freqs = fopen("lactase_frequencies.txt","w");
	}
	*/
	//for(t = ts; t > 0; t--) // Yes, > 0 (and not >= 0) because each iteration through the loop takes the current
	FILE * pedout;
	char fileName[1000];
	//int32_t maxRelationships, maxRelIdx;
	int32_t sumNumGenotypes;
	int32_t whichParent;
	FILE * max = fopen("ERROR_max_relationships.txt","w");
	for(t = ts-1; t >= 0; t--) 		// Already starts at ts-1, but that's A-OK there champ.
	{						// genotypes and produces the next [(t-1)'s] generation's genotypes.
		parentPopnSize = ped->popnSizes[t+1];
		currentPopnSize = ped->popnSizes[t];
		//numCompleted++;
		/*
		if(printFrequencies)
		{
			numSelectedAlleles = numGenotypes[0][1] + 2*numGenotypes[0][2] + numGenotypes[1][1] + 2*numGenotypes[1][2];
			freq = ((double)numSelectedAlleles)/(2*currentPopnSize);
			fprintf(freqs,"%f\n",freq);
		}
		*/

		/*
		sumNumGenotypes = 0;
		for(i = 0; i < 2; i++)
			for(j = 0; j < 3; j++)
				sumNumGenotypes += numGenotypes[i][j];
		if(sumNumGenotypes != popnSizes[t+1])
		{
			fprintf(stderr,"sumNumGenotypes = %i\npopnSizes[t+1] = %i\nt = %i\n",sumNumGenotypes,popnSizes[t+1],t);
			PERROR("Bad sumNumGenotypes!");
		}
		*/
		
		CoalPedigree_update_type_probs_(typeProbs[0],cTypeProbs[0],numGenotypes[0],genotypePerm[0],s,h);
		CoalPedigree_update_type_probs_(typeProbs[1],cTypeProbs[1],numGenotypes[1],genotypePerm[1],s,h);

		if(printLastPedigree && t < 10)
		{
			sprintf(fileName,"peds/relationships_%i",t);
			printf("t<10!\n");
			printf("%s\n",fileName);
			pedout = fopen(fileName,"w");
		}

		for(i = 0; i < currentPopnSize; i++)
		{
			rU = runif();
			for(j = 0; j < 3; j++)
			{
				if(rU <= cTypeProbs[0][j])
				{
					parentGenotypes[i][0] = genotypePerm[0][j];
					break;
				}
			}
			rU = runif();
			for(j = 0; j < 3; j++)
			{
				if(rU <= cTypeProbs[1][j])
				{
					parentGenotypes[i][1] = genotypePerm[1][j];
					break;
				}
			}

			for(j = 0; j < 2; j++)
			{
				whichParent = runifd(0,numGenotypes[j][parentGenotypes[i][j]]-1);
				ped->relationships[t][i][j] = genotypeIdxs[j][parentGenotypes[i][j]][whichParent];
				if(ped->relationships[t][i][j] > parentPopnSize-1)
				{
					fprintf(stderr,"bad relationships assignment.\n");
				}
			}

			if(printLastPedigree && t < 10)
			{
				fprintf(pedout, "%i\t%i\n", ped->relationships[t][i][0],ped->relationships[t][i][1]);
				fflush(pedout);
			}


			curGenotypes[i] = 0;
			for(j = 0; j < 2; j++)
			{
				switch(parentGenotypes[i][j])
				{
					case 0:
						break;
					case 2:
						curGenotypes[i] += 1;
						break;
					case 1:
						curGenotypes[i] += rbern(0.5);
						break;
					default:
						PERROR("Bad curGenotypes[i]");
				}
			}
		}

		/*	
		maxRelationships = 0;
		maxRelIdx = 0;
		for(i = 0; i < currentPopnSize; i++)
		{
			if(ped->relationships[t][i][1] > maxRelationships)
			{
				maxRelationships = ped->relationships[t][i][1];
				maxRelIdx = i;
			}
		}
		fprintf(max,"%i\n",maxRelationships);
		if(maxRelationships > popnSizes[t+1])
		{
			fprintf(stderr, "Bad maxRelationships\n");
			fprintf(stderr,"popnSizes[t-1] = %i\nmaxRelationships = %i\nt = %i\n",popnSizes[t-1],maxRelationships, t);
		}
		*/

		if(printLastPedigree && t < 10)
			fclose(pedout);

		// Now, update the genotypes.
		for(i = 0; i < 3; i++)
		{
			numGenotypes[0][i] = 0;
			numGenotypes[1][i] = 0;
		}

		for(i = 0; i < currentPopnSize/2; i++)
		{
			curGenotype = curGenotypes[i];
			genotypeIdxs[0][curGenotype][numGenotypes[0][curGenotype]] = i;
			(numGenotypes[0][curGenotype])++;
		}
		for(i = currentPopnSize/2; i < currentPopnSize; i++)
		{
			curGenotype = curGenotypes[i];
			genotypeIdxs[1][curGenotype][numGenotypes[1][curGenotype]] = i;
			(numGenotypes[1][curGenotype])++;
		}

		if((fixed[1] = (numGenotypes[0][2] == currentPopnSize/2 && numGenotypes[1][2] == currentPopnSize/2)) || 
				(fixed[0] = (numGenotypes[0][0] == currentPopnSize/2 && numGenotypes[1][0] == currentPopnSize/2)))
		{
			if(requireSweep && fixed[0])
			{
				/*
				if(printFrequencies)
				{
					fclose(freqs);
					freqs = fopen("lactase_frequencies.txt","w");
				}
				*/
				//currentPopnSize = popnSizes[ts-1];
				currentPopnSize = popnSizes[ts]; 		// TODO: Check this.
				numGenotypes[0][0] = currentPopnSize/2;
				numGenotypes[0][1] = 0;
				numGenotypes[0][2] = 0;

				numGenotypes[1][0] = currentPopnSize/2;
				numGenotypes[1][1] = 0;
				numGenotypes[1][2] = 0;
				if(initialPos < currentPopnSize/2) // it arises in a female...
				{
					numGenotypes[0][1]++;
					numGenotypes[0][0]--;
				}
				else // it arises in a male...
				{
					(numGenotypes[1][1])++;
					(numGenotypes[1][0])--;
				}
				curPos = 0;
				for(i = 0; i < currentPopnSize/2; i++)
				{
					if(initialPos != i)
					{
						genotypeIdxs[0][0][curPos++] = i;
					}
					else
					{
						genotypeIdxs[0][1][0] = i;
					}
				}

				curPos = 0;
				for(i = currentPopnSize/2; i < currentPopnSize; i++)
				{
					if(initialPos != i)
					{
						genotypeIdxs[1][0][curPos++] = i;
					}
					else
					{
						genotypeIdxs[1][1][0] = i;
					}
				}
				t = ts; // To be deincremented by the loop to ts-1 to start over.
				//printf("t = ts... starting over.\n");
			}
			else 		// If it fixed at 1, or it fixed at 0 and we don't care whether it fixed at 1 or 0...
			{
				int tp;
				finalGenotype = 2*fixed[1];
				for(tp = t-1; tp >= 0; tp--)
				{
					parentPopnSize = popnSizes[tp+1];
					currentPopnSize = popnSizes[tp];
					// This was a bug: parentPopnSize was currentPopnSize instead. Ugh.
					for(i = 0; i < currentPopnSize; i++)
					{
						ped->relationships[tp][i][0] = runifd(0,parentPopnSize/2-1);
						ped->relationships[tp][i][1] = runifd(parentPopnSize/2,parentPopnSize-1);
						//ped->genotypes[tp][i] = finalGenotype;	
					}
				}
				fclose(max); //TOREMOVE
				break; // breaks from the t loop.
			}
		}
	}
	free(parentGenotypes);	
	free(curGenotypes);	
	return;
}

/* Excluding a handful of "sweep" pedigree-shuffling functions for embedding a selective
 * sweep in the pedigree and then simulating segregation at a locus linked to the
 * selected site. */


/* simulates a coalescence time for chromosomes taken randomly from the
 * individuals in indiv.
 * ped		initialized pedigree
 * indiv	individuals from which the two chromosomes are sampled (length 2)
 */
int32_t CoalPedigree_sim_coal_time(CoalPedigree * ped, int32_t * indivs)
{
    int32_t indivsP[2];
    indivsP[0] = indivs[0];
    indivsP[1] = indivs[1];

	int32_t t=0,i;
	/* Cycle through the pedigree however many times until coalescence happens...
	 * But what happens if no pedigree common ancestor happens?
	 * This would be extremely rare. */
	while(1)	
	{
		for(i = 0; i < ped->numGenerations-1; i++) // Note numGenerations-1! This is because there often 
		{										   // isn't any information for (numGenerations-1)'s parents.
			if(indivsP[0] == indivsP[1])
			{
				/* If the two chromosomes are in the same parent, with probability
				 * 1/2 they go to the same parent, and otherwise they go to the two
				 * parents of that individual. */
				if(rbern(0.5) && t != 0) // if t == 0, we're sampling two chromosomes from the same individual, and we forbid them from coalescing
					return t+1;
				else
				{
					indivsP[0] = ped->relationships[i][indivsP[0]][0];
					indivsP[1] = ped->relationships[i][indivsP[1]][1];
					t++;
				}
			}
			else
			{
				/* If they're not in the same parent, go back one generation, picking a
				 * a parent randomly for each individual and increment t. */
				indivsP[0] = ped->relationships[i][indivsP[0]][rbern(0.5)];
				indivsP[1] = ped->relationships[i][indivsP[1]][rbern(0.5)];
				t++;
			}
		}
	}
}

// Just like CoalPedigree_sim_coal_time, but loops back through only a certain part of the
// pedigree, useful when say you have an initial period of exponential decline
// that you don't want to loop through many times. loopTime is the time
// (back in time) that the process should go back to when it's made it all the
// way through the pedigree.
int32_t CoalPedigree_sim_coal_time_loop(CoalPedigree * ped, int32_t * indiv, int32_t loopTime)
{
	int32_t t=0,i;
	for(i = 0; i < ped->numGenerations-1; i++) // Note numGenerations-1! This is because there often 
	{										   // isn't any information for (numGenerations-1)'s parents.
		if(indiv[0] == indiv[1])
		{
			/* If the two chromosomes are in the same parent, with probability
			 * 1/2 they go to the same parent, and otherwise they go to the two
			 * parents of that individual. */
			if(rbern(0.5))
				return t+1;
			else
			{
				indiv[0] = ped->relationships[i][indiv[0]][0];
				indiv[1] = ped->relationships[i][indiv[0]][1];
				t++;
			}
		}
		else
		{
			/* If they're not in the same parent, go back one generation, picking a
			 * a parent randomly for each individual and increment t. */
			indiv[0] = ped->relationships[i][indiv[0]][rbern(0.5)];
			indiv[1] = ped->relationships[i][indiv[1]][rbern(0.5)];
			t++;
		}
	}

	while(1)	
	{
		for(i = loopTime; i < ped->numGenerations-1; i++) // Note numGenerations-1! This is because there often 
		{										   // isn't any information for (numGenerations-1)'s parents.
			if(indiv[0] == indiv[1])
			{
				/* If the two chromosomes are in the same parent, with probability
				 * 1/2 they go to the same parent, and otherwise they go to the two
				 * parents of that individual. */
				if(rbern(0.5))
					return t+1;
				else
				{
					indiv[0] = ped->relationships[i][indiv[0]][0];
					indiv[1] = ped->relationships[i][indiv[0]][1];
					t++;
				}
			}
			else
			{
				/* If they're not in the same parent, go back one generation, picking a
				 * a parent randomly for each individual and increment t. */
				indiv[0] = ped->relationships[i][indiv[0]][rbern(0.5)];
				indiv[1] = ped->relationships[i][indiv[1]][rbern(0.5)];
				t++;
			}
		}
	}
}

void CoalPedigree_calculate_num_migrants(CoalPedigree * ped, int32_t popn, int32_t * countBucket, int32_t * numMigrants, int32_t * numOffspring)
{
	int32_t numGenerations = ped->numGenerations;
	int32_t demeSize = ped->popnSize/2;
	int32_t t,i,j;

	for(i = 0; i < numGenerations; i++)
	{
		numMigrants[i] = 0;
		numOffspring[i] = 0;
	}

	if(!popn)
	{
		for(t = 0; t < numGenerations; t++)
		{
			for(i = 0; i < demeSize; i++)
			{
				if(ped->relationships[t][i][0] >= demeSize)
				{
					numMigrants[t]++;
					if(t > 0)
					{
						if(i < demeSize/2)
						{
							for(j = 0; j < demeSize; j++)
							{
								if(ped->relationships[t-1][j][0] == i)
									numOffspring[t]++;
							}
						}
						else
						{
							for(j = 0; j < demeSize; j++)
							{
								if(ped->relationships[t-1][j][1] == i)
									numOffspring[t]++;
							}

						}
					}
				}
			}
		}
	}
	else if(popn == 1)
	{
		for(t = 0; t < numGenerations; t++)
		{
			for(i = demeSize; i < 2*demeSize; i++)
			{
				if(ped->relationships[t][i][0] < demeSize)
				{
					numMigrants[t]++;
					if(t > 0)
					{
						if(i < 3*demeSize/2)
						{
							for(j = demeSize; j < 2*demeSize-1; j++)
							{
								if(ped->relationships[t-1][j][0] == i)
									numOffspring[t]++;
							}
						}
						else
						{
							for(j = demeSize; j < 2*demeSize-1; j++)
							{
								if(ped->relationships[t-1][j][1] == i)
									numOffspring[t]++;
							}

						}
					}
				}
			}
		}
	}
	else
		PERROR("non 0/1 argument to popn in CoalPedigree_calculate_num_migrants");
	return;
}

/* This function returns the probability of escaping from a given subpopulation
 * within numGenerations generations given that the lineage starts in the two
 * individuals in indiv.
 * weights and weightsP are both pointers to double arrays of length N (deme size)
 *
 * Update: this should be redundant with
 * reconfiguration-probability-calculating functions.
 * */
double CoalPedigree_calculate_migration_prob(CoalPedigree * ped, int32_t numGenerations, int32_t popn, int32_t * indiv, double * weights, double * weightsP)
{
	int32_t i,t;
	int32_t N = ped->popnSize / 2;
	double totalProb = 0.0;
	for(i = 0; i < N; i++)
		weights[i] = weightsP[i] = 0.0;
	weights[indiv[0]] += 0.5;
	weights[indiv[1]] += 0.5;

	if(!popn)
	{
		for(t = 0; t < numGenerations; t++)
		{
			for(i = 0; i < N; i++)
			{
				if(weights[i] > 0.0)		// check to see if that weight is an escape, if so, ..., if not ...
				{
					if(ped->relationships[t][i][0] >= N) // both relationships will be >= N, so can just check mother relationship.
					{
						totalProb += weights[i];
						weights[i] = 0.0;
					}
					else
					{
						weightsP[ ped->relationships[t][i][0] ] += 0.5*weights[i];
						weightsP[ ped->relationships[t][i][1] ] += 0.5*weights[i];
					}
				}
			}
			for(i = 0; i < N; i++)
			{
				weights[i] = weightsP[i];
				weightsP[i] = 0.0;
			}
		}
	}
	else if(popn == 1)
	{
		for(t = 0; t < numGenerations; t++)
		{
			for(i = 0; i < N; i++)
			{
				if(weights[i] > 0.0)		// check to see if that weight is an escape, if so, ..., if not ...
				{
					if(ped->relationships[t][N+i][0] < N) // both relationships will be < N, so can just check mother relationship.
					{
						totalProb += weights[i];
						weights[i] = 0.0;
					}
					else
					{
						weightsP[ ped->relationships[t][N+i][0] ] += 0.5*weights[i];
						weightsP[ ped->relationships[t][N+i][1] ] += 0.5*weights[i];
					}
				}
			}
			for(i = 0; i < N; i++)
			{
				weights[i] = weightsP[i];
				weightsP[i] = 0.0;
			}
		}
	}
	else
		PERROR("non 0/1 subpopulation specification in CoalPedigree_calculate_prob_escape.");
	return totalProb;
}

// note still only accomodates one migrant.
double CoalPedigree_calculate_migration_prob_singleindiv(CoalPedigree * ped, int32_t numGenerations, int32_t indiv, double * weights, double * weightsP)
{
	int32_t i,t;
	int32_t N = ped->popnSize/2;
	int32_t popn = indiv >= N;
	double totalProb = 0.0;
	for(i = 0; i < N; i++)
		weights[i] = weightsP[i] = 0.0;
	weights[indiv] += 1.0;
	// TODO: make this able to accomodate tests in the second subpopulation
	if(!popn)
	{
		for(t = 0; t < numGenerations; t++)
		{
			for(i = 0; i < N; i++)
			{
				if(weights[i] > 0.0)		// check to see if that weight is an escape, if so, ..., if not ...
				{
					if(ped->relationships[t][i][0] >= N) // both relationships will be >= N, so can just check mother relationship.
					{
						totalProb += weights[i];
						weights[i] = 0.0;
					}
					else
					{
						weightsP[ ped->relationships[t][i][0] ] += 0.5*weights[i];
						weightsP[ ped->relationships[t][i][1] ] += 0.5*weights[i];
					}
				}
			}
			for(i = 0; i < N; i++)
			{
				weights[i] = weightsP[i];
				weightsP[i] = 0.0;
			}
		}
	}
	else if(popn == 1)
	{
		for(t = 0; t < numGenerations; t++)
		{
			for(i = 0; i < N; i++)
			{
				if(weights[i] > 0.0)		// check to see if that weight is an escape, if so, ..., if not ...
				{
					if(ped->relationships[t][N+i][0] < N) // both relationships will be < N, so can just check mother relationship.
					{
						totalProb += weights[i];
						weights[i] = 0.0;
					}
					else
					{
						weightsP[ ped->relationships[t][N+i][0] ] += 0.5*weights[i];
						weightsP[ ped->relationships[t][N+i][1] ] += 0.5*weights[i];
					}
				}
			}
			for(i = 0; i < N; i++)
			{
				weights[i] = weightsP[i];
				weightsP[i] = 0.0;
			}
		}
	}
	else
		PERROR("non 0/1 subpopulation specification in CoalPedigree_calculate_prob_escape.");
	return totalProb;
}

void CoalPedigree_calculate_migration_prob_cumulative(CoalPedigree * ped, int32_t numGenerations, int32_t popn, int32_t * indiv, double * weights, double * weightsP, double * cumProb)
{
	int32_t i,t;
	int32_t N = ped->popnSize / 2;
	//double totalProb = 0.0;
	for(i = 0; i < N; i++)
		weights[i] = weightsP[i] = 0.0;
	for(t = 0; t < numGenerations; t++)
		cumProb[t] = 0.0;

	weights[indiv[0]] += 0.5;
	weights[indiv[1]] += 0.5;

	// TODO: make this able to accomodate tests in the second subpopulation
	if(!popn)
	{
		for(t = 0; t < numGenerations; t++)
		{
			for(i = 0; i < N; i++)
			{
				if(weights[i] > 0.0)		// check to see if that weight is an escape, if so, ..., if not ...
				{
					if(ped->relationships[t][i][0] >= N) // both relationships will be >= N, so can just check mother relationship.
					{
						//totalProb += weights[i];
						cumProb[t] += weights[i];
						weights[i] = 0.0;
					}
					else
					{
						weightsP[ ped->relationships[t][i][0] ] += 0.5*weights[i];
						weightsP[ ped->relationships[t][i][1] ] += 0.5*weights[i];
					}
				}
			}
			for(i = 0; i < N; i++)
			{
				weights[i] = weightsP[i];
				weightsP[i] = 0.0;
			}
			if(t < numGenerations-1)
				cumProb[t+1] += cumProb[t];
		}
	}
	else if(popn == 1)
	{
		for(t = 0; t < numGenerations; t++)
		{
			for(i = 0; i < N; i++)
			{
				if(weights[i] > 0.0)		// check to see if that weight is an escape, if so, ..., if not ...
				{
					if(ped->relationships[t][N+i][0] < N) // both relationships will be < N, so can just check mother relationship.
					{
						//totalProb += weights[i];
						cumProb[t] += weights[i];
						weights[i] = 0.0;
					}
					else
					{
						weightsP[ ped->relationships[t][N+i][0] ] += 0.5*weights[i];
						weightsP[ ped->relationships[t][N+i][1] ] += 0.5*weights[i];
					}
				}
			}
			for(i = 0; i < N; i++)
			{
				weights[i] = weightsP[i];
				weightsP[i] = 0.0;
			}
			if(t < numGenerations-1)
				cumProb[t+1] += cumProb[t];
		}
	}
	else
		PERROR("non 0/1 subpopulation specification in CoalPedigree_calculate_prob_escape.");
	return;
}

void CoalPedigree_print_migrant_weights(CoalPedigree * ped, int32_t numGens, int32_t popn, double ** wc, double ** wcp, FILE * fout)
{
	int32_t i, j, t;
	int32_t N = ped->popnSize / 2;
	double weight = 0.0;

	// all the weight starts in the initial individual.
	for(i = 0; i < N; i++)
	{
		for(j = 0; j < N; j++)
		{
			wc[i][j] = 0.0;
			wcp[i][j] = 0.0;
		}
		wc[i][i] = 1.0;
	}
	
	if(popn == 0)
	{
		for(t = 0; t < numGens; t++)
		{
			// for each gen., first check migrant status for i, then update weights
			for(i = 0; i < N; i++)
			{
				if(ped->relationships[t][i][0] >= N)
				{
					for(j = 0; j < N; j++)		// note the switching of j and i, here.
					{
						weight += wc[j][i];
						if(wc[j][i] > 0.0)		// remove all of the weight going to the migrant.
							wc[j][i] = 0.0;
					}
					fprintf(fout, "%i\t%f\n", t, weight);
					weight = 0.0;
				}
				for(j = 0; j < N; j++)
				{
					wcp[i][ ped->relationships[t][j][0] ] += 0.5*(wc[i][j]);
					wcp[i][ ped->relationships[t][j][1] ] += 0.5*(wc[i][j]);
				}
				for(j = 0; j < N; j++)
				{
					wc[i][j] = wcp[i][j];
					wcp[i][j] = 0.0;
				}
			}
		}
	}
	else if(popn == 1)
	{
		for(t = 0; t < numGens; t++)
		{
			// for each gen., first check migrant status for i, then update weights
			for(i = 0; i < N; i++)
			{
				if(ped->relationships[t][N+i][0] < N)
				{
					for(j = 0; j < N; j++)		// note the switching of j and i, here.
					{
						weight += wc[j][i];
						if(wc[j][i] > 0.0)		// remove all of the weight going to the migrant.
							wc[j][i] = 0.0;
					}
					fprintf(fout, "%i\t%f\n", t, weight);
					weight = 0.0;
				}
				for(j = 0; j < N; j++)
				{
					wcp[i][ ped->relationships[t][j][0] ] += 0.5*(wc[i][j]);
					wcp[i][ ped->relationships[t][j][1] ] += 0.5*(wc[i][j]);
				}
				for(j = 0; j < N; j++)
				{
					wc[i][j] = wcp[i][j];
					wcp[i][j] = 0.0;
				}
			}
		}
	}
	else
		PERROR("non-0/1 value for popn in CoalPedigree_print_migrant_weights");
	return;
}

/*	This function takes a D-deme pedigree and a set of individuals in the same
 *	deme and determines the total migrant weight for those individuals in the
 *	first numGens generations. 
 *
 * */
void CoalPedigree_print_migrant_weights_ddeme(CoalPedigree * ped, int32_t * indiv, int32_t numIndivs, int32_t numGens, double ** wc, double ** wcp, FILE * fout)
{
	int32_t i, j, t, demeIdx, firstDemeFemale, lastDemeFemale;
	int32_t N = ped->popnSize / ped->numDemes;
	double * weight = (double *)malloc(sizeof(double) * numIndivs);
	CHECKPOINTER(weight);
	// all the weight starts in the initial individuals.
	demeIdx = indiv[0] / N;				// again, assuming all indivs in the same deme.
	firstDemeFemale = N*demeIdx;
	lastDemeFemale = firstDemeFemale + N/2;

	for(j = 0; j < numIndivs; j++)
	{
		for(i = 0; i < N; i++)
		{
			wc[j][i] = 0.0;
			wcp[j][i] = 0.0;
		}
		wc[j][ indiv[j]-firstDemeFemale ] = 1.0;		// for the jth set of weights, all the weight starts in indiv[j].
		weight[j] = 0.0;
	}
	
	for(t = 0; t < numGens; t++)
	{
		// for each gen., first check migrant status for i, then update weights
		for(i = 0; i < N; i++)
		{
			if(ped->relationships[t][firstDemeFemale+i][0] < firstDemeFemale || ped->relationships[t][firstDemeFemale+i][0] > lastDemeFemale)
			{
				for(j = 0; j < numIndivs; j++)		// note the switching of j and i, here.
				{
					if(wc[j][i] > 0.0)		// remove all of the weight going to the migrant.
					{
						weight[j] += wc[j][i];
						wc[j][i] = 0.0;
					}
				}
			}
		}
		for(j = 0; j < numIndivs; j++)
		{
			for(i = 0; i < N; i++)
			{
				wcp[j][ ped->relationships[t][firstDemeFemale+i][0]-firstDemeFemale ] += 0.5*(wc[j][i]);
				wcp[j][ ped->relationships[t][firstDemeFemale+i][1]-firstDemeFemale ] += 0.5*(wc[j][i]);
			}
			for(i = 0; i < N; i++)
			{
				wc[j][i] = wcp[j][i];
				wcp[j][i] = 0.0;
			}
		}
	}
	for(j = 0; j < numIndivs; j++)
		fprintf(fout, "%f\t", weight[j]);
	fprintf(fout, "\n");		// likely this will be called again with a different set of individuals, for a different deme and a different line.
	free(weight);
	return;
}

/* This function calculates the probability of a pair coalescing during each of
 * the first maxGen generations.  weights and weightsP must be
 * popnSize^2-length array of doubles, since positions need to be dealt with
 * jointly. probs is a maxGen-length double vector into which the coalescence
 * probabilities go.*/
void CoalPedigree_calc_coal_probs(CoalPedigree * ped, int32_t * indivs, int32_t maxGen, double * weights, double * weightsP, double * probs)
{
	int32_t t, i, j, k, indivIdx, parentCombIdxs[4], popnSize = ped->popnSize, checkInterval = maxGen/30, anyWeight;
	double prob = 0.0;
	for(i = 0; i < popnSize*popnSize; i++)
		weights[i] = weightsP[i] = 0.0;
	weights[popnSize*indivs[0] + indivs[1]] = 1.0;		// initial configuration.
	for(t = 0; t < maxGen-1; t++)		// note < not <=. Arbitrary.
	{
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
		probs[t] = prob;
		prob = 0.0;
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
			if(!anyWeight)		// if there's no weight left (due to the limits of double variables), return
				return;
		}
			

	}
	return;
}
/* This function prints the probability of a pair coalescing in the first maxGen generations.
 * weights and weightsP must be popnSize^2-length array of doubles, since positions need to
 * be dealt with jointly. */
void CoalPedigree_print_coal_probs(CoalPedigree * ped, int32_t * indivs, int32_t maxGen, double * weights, double * weightsP, FILE * fout)
{
	int32_t t, i, j, k, indivIdx, parentCombIdxs[4], popnSize = ped->popnSize, checkInterval = maxGen/30, anyWeight;
	double prob = 0.0;
	for(i = 0; i < popnSize*popnSize; i++)
		weights[i] = weightsP[i] = 0.0;
	weights[popnSize*indivs[0] + indivs[1]] = 1.0;		// initial configuration.
	for(t = 0; t < maxGen-1; t++)		// note < not <=. Arbitrary.
	{
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
		fprintf(fout, "%i\t%f\n", t+1, prob);
		prob = 0.0;
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
			if(!anyWeight)		// if there's no weight left (due to the limits of double variables), return
				return;
		}
	}
	return;
}

/* Exports a pedigree into a .csv file, with the first line being various
 * parameters of the pedigree. First a space-separated header is printed, then
 * the relationships are printed, started with the relationships of generation
 * 0 (present generation) back to its parents, generation 1, and so on back in
 * time.
 *
 * Currently only supports constant-sized populations. A 30*1000-generation
 * pedigree for a population of size 1000 produces a pedigree file that is
 * approximately 225 MB. Any larger (say, 10000 individuals for 30*10000
 * generations) will probably not be feasible. But that wouldn't really be
 * feasible with memory anyway.
 *
 * ped    initialized, shuffled pedigree
 * fout   output file
 *
 * */
void CoalPedigree_export_pedigree(CoalPedigree * ped, char * filename)
{
	int32_t t, i;
	FILE * fout = fopen(filename, "w");
	if(!fout)
	{
		fprintf(stderr, "Could not open %s to export pedigree.\n", filename);
		exit(1);
	}

	// header 
	fprintf(fout, "#");		// common comment symbol
	fprintf(fout, " %i %i %i %i %i\n", ped->numGenerations, ped->numMales, ped->numFemales, ped->popnSize, ped->numDemes);

	// relationships
	for(t = 0; t < ped->numGenerations-1; t++)
	{
		for(i = 0; i < ped->popnSize-1; i++)
			fprintf(fout, "%i,%i,", ped->relationships[t][i][0], ped->relationships[t][i][1]);
		fprintf(fout, "%i,%i\n", ped->relationships[t][ped->popnSize-1][0], ped->relationships[t][ped->popnSize-1][1]);
	}
	fclose(fout);
	return;
}

/* quick note: ped should be *uninitialized*. */
void CoalPedigree_import_pedigree(CoalPedigree * ped, char * filename)
{
	int32_t t, i;
	char * line = malloc(sizeof(char) * 100000);
	char * linecopy = malloc(sizeof(char) * 100000);
	CHECKPOINTER(line);
	FILE * fin = fopen(filename, "r");
	if(!fin)
	{
		fprintf(stderr, "Could not locate file %s for importing.\n", filename);
		exit(1); // requires stdlib.h
	}
	line = fgets(line, 100000, fin);

	// for reference (order):
	//fprintf(fout, " %i %i %i %i %i\n", ped->numGenerations, ped->numMales, ped->numFemales, ped->popnSize, ped->numDemes);
	
	// first line
	char * curTok;
	curTok = strtok(line, " ");    // #
	// numGenerations
	curTok = strtok(NULL, " ");
	ped->numGenerations = strtol(curTok, NULL, 10);
	// numMales
	curTok = strtok(NULL, " ");
	ped->numMales = strtol(curTok, NULL, 10);
	// numFemales
	curTok = strtok(NULL, " ");
	ped->numFemales = strtol(curTok, NULL, 10);
	// popnSize
	curTok = strtok(NULL, " ");
	ped->popnSize = strtol(curTok, NULL, 10);
	// numDemes
	curTok = strtok(NULL, " ");
	ped->numDemes = strtol(curTok, NULL, 10);

	CoalPedigree_init_singledeme(ped, ped->numGenerations, ped->numFemales, ped->numMales);

	t = 0;
	while(fgets(line, 100000, fin))
	{
		curTok = strtok(line, ",");
		ped->relationships[t][0][0] = strtol(curTok, NULL, 10);
		curTok = strtok(NULL, ",");
		ped->relationships[t][0][1] = strtol(curTok, NULL, 10);
		for(i = 1; i < ped->popnSize; i++)
		{
			curTok = strtok(NULL, ",");
			ped->relationships[t][i][0] = strtol(curTok, NULL, 10);
			curTok = strtok(NULL, ",");
			ped->relationships[t][i][1] = strtol(curTok, NULL, 10);
		}
		t++;
	}
	// relationships
	/*
	while(fgets(line, 100000, fin) != NULL)
	{
	}
	*/
	fclose(fin);
	return;
}

/*
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
*/


