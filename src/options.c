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
#include <string.h>
#include <getopt.h>
#include "definitions.h"
#include "options.h"


static struct option long_options[] =
{
	{"help", no_argument, 0, 'h'},
	{"theta",  required_argument, 0, 't'},
	{"migration-rate", required_argument, 0, 'M'},
	{"num-demes", required_argument, 0, 'D'},
	{"popn-size", required_argument, 0, 'N'},
	{"deme-size", required_argument, 0, 'N'},
	{"sample-size", required_argument, 0, 'n'},
	{"popn-size-file", required_argument, 0, 1},
	{"num-loci", required_argument, 0, 'k'},
	{"genealogies", no_argument, 0, 'g'},
	{"ms-out", required_argument, 0, 3},
	{"genepop-out", required_argument, 0, 4},
	{"sweep", no_argument, 0, 11},
	{"sweep-s", required_argument, 0, 13}, // selective advantage of advantageous allele
	{"sweep-h", required_argument, 0, 14}, // dominance of advantageous allele
	{"sweep-start-gen", required_argument, 0, 15}, // generation in which sweep starts
	{"require-fixation", no_argument, 0, 6},    // does the selective sweep have to fix?
	{"cyclical", no_argument, 0, 7},
	{"monoecious", no_argument, 0, 8},
	{"export-pedigree", required_argument, 0, 9},
	{"import-pedigree", required_argument, 0, 10},
	{"num-gens", required_argument, 0, 12},
	{"popn-size-file", required_argument, 0, 16},
	{"loop-generation", required_argument, 0, 17},
	{"rearrangement-probabilities", required_argument, 0, 18},
	{"rearrangement-generation", required_argument, 0, 19},
	{"rearrangement-reps", required_argument, 0, 20},
	{"fixed-ibd", required_argument, 0, 21},
	{"fixed-migrant", required_argument, 0, 22},
    {"num-coal-times", required_argument, 0, 23},
    {"exact", no_argument, 0, 24},
    {"within", no_argument, 0, 25},
    {"between", no_argument, 0, 26},
    {"fst", no_argument, 0, 27},
    {"num-pairs", required_argument, 0, 28},
    {"shared-ancestor", required_argument, 0, 29},
    {"single-lineage", no_argument, 0, 30},
    {"num-single-lineages", required_argument, 0, 31},
    {"means", no_argument, 0, 32},
    {"num-prob-gens", required_argument, 0, 33},
	{0, 0, 0, 0}
};

char usage[] = "\ncoalped [OPTIONS]\n\
    \n\
    -t, --theta                     mutation parameter 4*N*mu, where N is (diploid) population or\n\
                                    deme size and mu is the per-generation mutation probability\n\
    -M, --migration-rate            mutation parameter 4*N*m, where N is (diploid) population or deme size\n\
    -D, --num-demes                 number of demes or subpopulations\n\
    -N, --popn-size, --deme-size    size of each deme (or size of population when number of demes is 1)\n\
    -n, --sample-size               number of diploid individuals to sample from each deme\n\
    -k, --num-loci                  number of replicate loci to simulate\n\
    -g, --genealogies               display genealogies in Newick format in output\n\
    --popn-size-file                file containing population sizes for variable population sizes\n\
    --genepop-out                   filename for genepop output\n\
    --sweep                         indicate that a selective sweep is to be simulated (requires --sweep-s, --sweep-h, --sweep-start-gen)\n\
    --sweep-s                       selective advantage for beneficial allele (homozygote fitness: 1+s; heterozygote: 1+s*h)\n\
    --sweep-h                       dominance of positively selected allele\n\
    --sweep-start-gen               number of generations before present when selected allele arose\n\
    --require-fixation              require fixation of selected sweep before present (default: false)\n\
    --cyclical                      simulate a cyclical pedigree\n\
    --monoecious                    simulate a monoecious pedigree\n\
    --import-pedigree               filename for file containing pedigree to import\n\
    --export-pedigree               filename for file to contain exported pedigree\n\
    --num-gens                      number of generations to simulate pedigree (default 30 * deme size)\n\
    --loop-generation               generation to loop to when running through pedigree (to avoid sweep or change in deme size)\n\
    --rearrangement-probabilities   filename to which the rearrangement probabilities will be written\n\
    --rearrangement-generation      number of generations to keep track of when simulating rearrangement probabilities\n\
    --rearrangement-reps            number of repetitions to simulate rearrangement probabilities\n\
    --fixed-ibd x                   have a single inbred individual whose\n\
                                    parents parents are related by a single\n\
                                    ancestor x generations before the present\n\
                                    (i.e., generation x, 0-indexed)\n\
    --fixed-migrant x               have a single individual with recent\n\
                                    migrant ancestry, who had one ancestor who\n\
                                    was a migrant x generations before the present\n\
    --shared-ancestor x             include a pair of individuals who share an ancestor\n\
                                    x generations before the present.\n\
	\n";
void Options_print_help_statement()
{
	printf("%s", &usage[0]);
	return;
}

void Options_set_defaults(Options * opt)
{
	opt->migRate = -1.0;
	opt->theta = -1.0;
	opt->numDemes = 1;
	opt->demeSize = -1;
	opt->sampleSize = -1;
	opt->importSet = 0;
	opt->exportSet = 0;
	opt->monoecious = 0;
	opt->cyclical = 0;
	opt->multiplePopnSizes = 0;
	opt->genealogiesSet = 0;
	opt->msSet = 0;
	opt->numLoci = -1;
	opt->genepopSet = 0;
	opt->infAllelesSet = 0;
	opt->sweepSet = 0;
	opt->requireFixation = 0;
	opt->sweep_s = -1.0;
	opt->s_set = 0;
	opt->h_set = 0;
	opt->ts_set = 0;
	opt->numGens = 0;
	opt->loopGen = 0;
	opt->rearrSet = 0;
	opt->rearrGen = 3;
	opt->rearrReps = 100000;
    opt->fixedIBDset = 0;
    opt->fixedIBDgen = -1;
    opt->fixedIBDindiv = -1;
    opt->fixedMigrantSet = 0;
    opt->fixedMigrantGen = -1;
    opt->fixedMigrantIndiv = -1;
    opt->exactSet = 0;
    opt->numCoalTimes = 1000;
    opt->fstSet = 0;
    opt->numPairs = 1;
    opt->multiplePairs = 0;
    opt->sharedAncSet = 0;
    opt->sharedAncGen = -1;
    opt->sharedAncIndivs[0] = opt->sharedAncIndivs[1] = -1;
    opt->withinSet = 1;
    opt->fstSet = 0;
    opt->numPairs = 1;
    opt->singleLineageSet = 0;
    opt->numSingleLineages = 1;
    opt->pedsimNumProbGens = 0;
    return;
}



void Options_parse_options(int argc, char ** argv, Options * opt)
{
	int32_t optionIndex, success;
	char c;
	// set default values
    Options_set_defaults(opt);

    if(argc == 1)
    {
        Options_print_help_statement();
        exit(0);
    }


	while(1)
	{
		c = getopt_long(argc, argv, "N:hi:D:t:M:o:gn:k:", long_options, &optionIndex);
		if(c == -1)
			break;
		switch(c)
		{
			case 'h':
				Options_print_help_statement();
                exit(0);
				break;
			case 't':
				success = sscanf(optarg, "%lf", &(opt->theta));
				if(!success)
					PERROR("Could not read theta (--theta, -t).\n");
				break;
			case 'M':
				success = sscanf(optarg, "%lf", &(opt->migRate));
				if(!success)
					PERROR("Could not read migration rate (--migrate, -M).\n");
				break;
			case 'D':
				success = sscanf(optarg, "%i", &(opt->numDemes));
				if(!success)
					PERROR("Could not read number of demes (--num-demes, -D).\n");
				break;
			case 'N':
				success = sscanf(optarg, "%i", &(opt->demeSize));
				if(!success)
					PERROR("Could not read population / deme size (--popn-size, --deme-size, -N).\n");
				break;
			case 'n':
				success = sscanf(optarg, "%i", &(opt->sampleSize));
				if(!success)
					PERROR("Could not read sample size.\n");
				break;
			case 'k':
				success = sscanf(optarg, "%i", &(opt->numLoci));
				if(!success)
					PERROR("Could not read number of loci.\n");
				break;
			case 1: // --popn-size-file
				success = sscanf(optarg, "%s", opt->popnSizeFile);
				opt->multiplePopnSizes = 1;
				if(!success)
					PERROR("Invalid population-size file.\n");
				break;
			case 'g': // --genealogies
				opt->genealogiesSet = 1;
				break;
			case 4: // --genepop-out
				opt->genepopSet = 1;
				strncpy(opt->geneFilename, optarg, 200);
				break;
			case 11: // --sweep
				opt->sweepSet = 1;
				break;
			case 13: // --sweep-s
				opt->s_set = 1;
				success = sscanf(optarg, "%lf", &(opt->sweep_s)); 
				if(!success)
					PERROR("Invalid selective coefficient."); 
				break;
			case 14: // --sweep-s
				opt->h_set = 1;
				success = sscanf(optarg, "%lf", &(opt->sweep_h)); 
				if(!success)
					PERROR("Invalid dominance of advantageous allele (--sweep-h)."); 
				break;
			case 15: // --sweep-start-gen
				opt->ts_set = 1;
				success = sscanf(optarg, "%i", &(opt->sweep_ts)); 
				if(!success)
					PERROR("Invalid starting generation for selective sweep");
				break;
			case 6: //--require-fixation
				opt->requireFixation = 1;
				break;
			case 7: // --cyclical
				opt->cyclical = 1;
				break;
			case 8: // --monoecious
				opt->monoecious = 1;
				break;
			case 9: // --export-pedigree
				opt->exportSet = 1;
				strncpy(opt->exportFilename, optarg, 200);
				break;
			case 10: // --import-pedigree
				opt->importSet = 1;
				strncpy(opt->importFilename, optarg, 200);
				break;
			case 12: // --num-gens
				success = sscanf(optarg, "%i", &(opt->numGens)); 
				if(!success)
					PERROR("Invalid number of generations (--num-gens)."); 
                break;
			case 17: // --loop-generation
				success = sscanf(optarg, "%i", &(opt->loopGen));
				if(!success)
					PERROR("Invalid loop generation (--loop-gen).");
                break;
			case 18: // --rearrangement-probabilities
				opt->rearrSet = 1;
				success = sscanf(optarg, "%s", &(opt->rearrFilename[0]));
				if(!success)
					PERROR("Invalid rearrangement filename (--rearrangement-probabilities).");
				break;
			case 19: // --rearrangement-generation
				success = sscanf(optarg, "%i", &(opt->rearrGen));
				if(!success)
					PERROR("Invalid rearrangement generation (--rearrangement-generation).");
				break;
			case 20: // --rearrangement-reps
				success = sscanf(optarg, "%i", &(opt->rearrReps));
				if(!success)
					PERROR("Invalid rearrangement reps (--rearrangement-reps).");
				break;
			case 21: // --fixed-ibd
                opt->fixedIBDset = 1;
				success = sscanf(optarg, "%i", &(opt->fixedIBDgen));
				if(!success)
					PERROR("Invalid IBD generation (--fixed-ibd).");
				break;
			case 22: // --fixed-migrant
                opt->fixedMigrantSet = 1;
				success = sscanf(optarg, "%i", &(opt->fixedMigrantGen));
				if(!success)
					PERROR("Invalid fixed-migrant generation (--fixed-migrant).");
                break;
			case 29: // --shared-ancestor
                opt->sharedAncSet = 1;
				success = sscanf(optarg, "%i", &(opt->sharedAncGen));
				if(!success)
					PERROR("Invalid shared-ancestor generation (--shared-ancestor).");
                break;
			default:
				Options_print_help_statement();
				exit(1);
				break;
		}
	}

	/* check validity/compatibility of options */
    if(opt->numLoci < 0)
        PERROR("--num-loci (-k) must be specified and > 0");
	if(opt->numDemes <= 0)
		PERROR("Number of demes (--num-demes, -D) must be greater than or equal to zero");
	if(opt->demeSize <= 0)
		PERROR("Deme size (--deme-size, -N) must be set to a positive number");
	if(opt->sampleSize <= 0)
		PERROR("Sample size (--sample-size, -n) must be set to a positive number");
	if(opt->numLoci <= 0 && (opt->msSet || opt->infAllelesSet || opt->genealogiesSet || opt->genepopSet))
		PERROR("Must specify number of loci");
	if(opt->sampleSize == -1)
		PERROR("Must specify a sample size (--sample-size, -n)");
	if((opt->msSet && opt->genepopSet) || (opt->msSet && opt->infAllelesSet) || (opt->genepopSet && opt->infAllelesSet))
		PERROR("Must specify just one option for type of locus to output: --ms-out, --genepop-out, --inf-alleles-out");
	if(opt->numDemes > 1 && opt->migRate == -1.0)
		PERROR("Must specify --migration-rate (-M) when --num-demes > 1");
	if((opt->msSet || opt->genepopSet || opt->infAllelesSet) && opt->theta <= 0.0)
		PERROR("Must specify --theta (-t) > 0");
    if(opt->fixedMigrantSet && opt->fixedMigrantGen < 1)
		PERROR("--fixed-migrant argument (generation) must be >= 1");
    if(opt->fixedIBDset && opt->fixedIBDgen < 2)
        PERROR("--fixed-ibd generation must be >= 2.");
    if(opt->sharedAncSet && opt->sharedAncGen < 1)
        PERROR("--shared-ancestor generation must be >= 1.");
    if(opt->sharedAncSet && opt->sampleSize < 2)
        PERROR("--shared-ancestor requires --sample-size (-n) >= 2");
	// (exclusive options)
	int32_t numExclusiveOptions = 0;
	numExclusiveOptions += opt->monoecious;
	numExclusiveOptions += opt->cyclical;
	numExclusiveOptions += opt->multiplePopnSizes;
	numExclusiveOptions += (opt->numDemes > 1);
	numExclusiveOptions += opt->sweepSet;
	numExclusiveOptions += opt->fixedIBDset;
    numExclusiveOptions += opt->sharedAncSet;
	if(numExclusiveOptions > 1)
		PERROR("The following options are mutually exclusive:\n\
                --sweep, --cyclical, --popn-size-file, --monoecious, --fixed-ibd, --shared-ancestor, and --num-demes > 1");
    numExclusiveOptions = 0;
	numExclusiveOptions += opt->monoecious;
	numExclusiveOptions += opt->cyclical;
	numExclusiveOptions += opt->multiplePopnSizes;
	numExclusiveOptions += opt->sweepSet;
	numExclusiveOptions += opt->fixedMigrantSet;
    if(numExclusiveOptions > 1)
		PERROR("The following options are mutually exclusive: --sweep, --cyclical, --popn-size-file, --monoecious, --fixed-migrant");
    if(opt->fixedMigrantSet && opt->numDemes != 2)
        PERROR("--num-demes must be 2 when using --fixed-migrant");
	if(opt->numGens < 0)
		PERROR("Number of generations (--num-gens) must be greater than zero.");
	// sweep settings:
	if(opt->sweepSet)
	{
		if(!opt->s_set || !opt->h_set || !opt->ts_set)
			PERROR("Must specify selective advantage (--sweep-s), dominance (--sweep-h), and starting generation (--sweep-start-gen) to call selective sweep.");
	}
	return;
}

char usage_pedsim[] = "\npedsim [OPTIONS]\n\
	\n\
    --exact                         calculate coalescence time distribution exactly\n\
    --means                         calculate means instead of distributions\n\
    --num-coal-times                number of coalescence times to simulate\n\
    -M, --migration-rate            mutation parameter 4*N*m, where N is (diploid) population or deme size\n\
    -D, --num-demes                 number of demes or subpopulations\n\
    --within                        simulate coalescence times within demes if D>1\n\
    --between                       simulate coalescence times between demes if D>1\n\
    --fst                           simulate FST = 1-mean(within)/mean(between)\n\
    -N, --popn-size, --deme-size    size of each deme (or size of population when number of demes is 1)\n\
    --popn-size-file                file containing population sizes for variable population sizes\n\
    --sweep                         indicate that a selective sweep is to be simulated (requires --sweep-s, --sweep-h, --sweep-start-gen)\n\
    --sweep-s                       selection coefficient for positively selected allele\n\
    --sweep-h                       dominance of positively selected allele\n\
    --sweep-start-gen               number of generations before present when selected allele arose\n\
    --require-fixation              require fixation of selected sweep before present (default: false)\n\
    --cyclical                      simulate a cyclical pedigree\n\
    --monoecious                    simulate a monoecious pedigree\n\
    --import-pedigree               filename for file containing pedigree to import\n\
    --export-pedigree               filename for file to contain exported pedigree\n\
    --num-gens                      number of generations to simulate pedigree (default 30 * deme size)\n\
    --loop-generation               generation to loop to when running through pedigree (to avoid sweep or change in deme size)\n\
    --rearrangement-probabilities   filename to which the rearrangement probabilities will be written\n\
    --rearrangement-generation      number of generations to keep track of when simulating rearrangement probabilities\n\
    --rearrangement-reps            number of repetitions to simulate rearrangement probabilities\n\
    --fixed-ibd x                   have a single inbred individual whose parents parents are related by a single ancestor x generations ago\n\
    --fixed-migrant x               have a single individual with recent migrant ancestry, who had one ancestor who was a migrant x generations ago\n\
    --num-pairs x                   simulate x random pairs of individuals\n\
    --num-prob-gens x               number of generations for which coalescence probabilities are calculated\n\
	\n";

void Options_print_help_statement_pedsim()
{
	printf("%s", &usage_pedsim[0]);
	exit(0);
	return;
}

void Options_parse_options_pedsim(int argc, char ** argv, Options * opt)
{
	int32_t optionIndex, success;
	char c;
	// set default values
    Options_set_defaults(opt);

	while(1)
	{
		c = getopt_long(argc, argv, "N:hi:D:t:M:", long_options, &optionIndex);
		if(c == -1)
			break;
		switch(c)
		{
			case 'h':
				Options_print_help_statement_pedsim();
				break;
			case 'M':
				success = sscanf(optarg, "%lf", &(opt->migRate));
				if(!success)
					PERROR("Could not read migration rate (--migrate, -M).\n");
				break;
			case 'D':
				success = sscanf(optarg, "%i", &(opt->numDemes));
				if(!success)
					PERROR("Could not read number of demes (--num-demes, -D).\n");
				break;
			case 'N':
				success = sscanf(optarg, "%i", &(opt->demeSize));
				if(!success)
					PERROR("Could not read population / deme size (--popn-size, --deme-size, -N).\n");
				break;
			case 1: // --popn-size-file
				success = sscanf(optarg, "%s", opt->popnSizeFile);
				opt->multiplePopnSizes = 1;
				if(!success)
					PERROR("Invalid population-size file.\n");
				break;
			case 11: // --sweep
				opt->sweepSet = 1;
                opt->sweepSet = 1;
				break;
			case 13: // --sweep-s
				opt->s_set = 1;
                opt->sweepSet = 1;
				success = sscanf(optarg, "%lf", &(opt->sweep_s)); 
				if(!success)
					PERROR("Invalid selective coefficient."); 
				break;
			case 14: // --sweep-s
				opt->h_set = 1;
                opt->sweepSet = 1;
				success = sscanf(optarg, "%lf", &(opt->sweep_h)); 
				if(!success)
					PERROR("Invalid dominance of advantageous allele (--sweep-h)."); 
				break;
			case 15: // --sweep-start-gen
				opt->ts_set = 1;
                opt->sweepSet = 1;
				success = sscanf(optarg, "%i", &(opt->sweep_ts)); 
				if(!success)
					PERROR("Invalid starting generation for selective sweep");
				break;
			case 6: //--require-fixation
				opt->requireFixation = 1;
                opt->sweepSet = 1;
				break;
			case 7: // --cyclical
				opt->cyclical = 1;
				break;
			case 8: // --monoecious
				opt->monoecious = 1;
				break;
			case 9: // --export-pedigree
				opt->exportSet = 1;
				strncpy(opt->exportFilename, optarg, 200);
				break;
			case 10: // --import-pedigree
				opt->importSet = 1;
				strncpy(opt->importFilename, optarg, 200);
				break;
			case 12: // --num-gens
				success = sscanf(optarg, "%i", &(opt->numGens)); 
				if(!success)
					PERROR("Invalid number of generations (--num-gens)."); 
                break;
			case 17: // --loop-generation
				success = sscanf(optarg, "%i", &(opt->loopGen));
				if(!success)
					PERROR("Invalid loop generation (--loop-gen).");
                break;
			case 18: // --rearrangement-probabilities
				opt->rearrSet = 1;
				success = sscanf(optarg, "%s", &(opt->rearrFilename[0]));
				if(!success)
					PERROR("Invalid rearrangement filename (--rearrangement-probabilities).");
				break;
			case 19: // --rearrangement-generation
				success = sscanf(optarg, "%i", &(opt->rearrGen));
				if(!success)
					PERROR("Invalid rearrangement generation (--rearrangement-generation).");
				break;
			case 20: // --rearrangement-reps
				success = sscanf(optarg, "%i", &(opt->rearrReps));
				if(!success)
					PERROR("Invalid rearrangement reps (--rearrangement-reps).");
				break;
			case 21: // --fixed-ibd
                opt->fixedIBDset = 1;
				success = sscanf(optarg, "%i", &(opt->fixedIBDgen));
				if(!success)
					PERROR("Invalid IBD generation (--fixed-ibd).");
				break;
			case 22: // --fixed-migrant
                opt->fixedMigrantSet = 1;
				success = sscanf(optarg, "%i", &(opt->fixedMigrantGen));
				if(!success)
					PERROR("Invalid fixed-migrant generation (--fixed-migrant).");
                break;
			case 23: // --num-simulated-times
				success = sscanf(optarg, "%i", &(opt->numCoalTimes));
				if(!success)
					PERROR("Invalid fixed-migrant generation (--fixed-migrant).");
                break;
			case 24: // --exact
                opt->exactSet = 1;
                break;
			case 25: // --within
                opt->withinSet = 1;
                // (within by default, also)
                break;
			case 26: // --between
                opt->withinSet = 0;
                // (within by default)
                break;
			case 27: // --fst
                opt->fstSet = 1;
                opt->withinSet = 0;
                break;
			case 28: // --num-pairs
				success = sscanf(optarg, "%i", &(opt->numPairs));
                if(opt->numPairs > 1)
                    opt->multiplePairs = 1;
                if(!success)
                    PERROR("Invalid number of pairs (--num-pairs)");
                if(opt->numPairs <= 0)
                    PERROR("--num-pairs must be >= 1");
                break;
			case 29: // --shared-ancestor
                opt->sharedAncSet = 1;
				success = sscanf(optarg, "%i", &(opt->sharedAncGen));
				if(!success)
					PERROR("Invalid shared-ancestor generation (--shared-ancestor).");
                break;
            case 30: // --single-lineage
                opt->singleLineageSet = 1;
                break;
            case 31: // --num-single-lineages
                success = sscanf(optarg, "%i", &(opt->numSingleLineages));
                if(!success)
                    PERROR("Invalid number of single lineages (--num-single-lineages)");
                break;
            case 32: // --means
                opt->exactSet = 0;
                break;
            case 33:
                success = sscanf(optarg, "%i", &(opt->pedsimNumProbGens));
                if(!success)
                    PERROR("invalid --num-prob-gens");
                break;
			default:
                if(optarg)
                    fprintf(stderr, "\nInvalid option: %s\n\n", optarg);
				Options_print_help_statement_pedsim();
				exit(1);
				break;
		}
	}

	/* check validity/compatibility of options */
	if(opt->numDemes <= 0)
		PERROR("Number of demes (--num-demes, -D) must be greater than or equal to zero.");
	if(opt->demeSize <= 0)
		PERROR("Deme size (--deme-size, -N) must be set to a positive integer.");
	if(opt->numDemes > 1 && opt->migRate == -1.0)
		PERROR("Must specify --migration-rate (-M) when --num-demes > 1.");
	// (exclusive options)
	int32_t numExclusiveOptions = 0;
	numExclusiveOptions += opt->monoecious;
	numExclusiveOptions += opt->cyclical;
	numExclusiveOptions += opt->multiplePopnSizes;
	numExclusiveOptions += (opt->numDemes > 1);
	numExclusiveOptions += opt->sweepSet;
	numExclusiveOptions += opt->fixedIBDset;
	if(numExclusiveOptions > 1)
		PERROR("The following options are mutually exclusive: --sweep, --cyclical, --popn-size-file, --monoecious, --fixed-ibd, and --num-demes > 1");
    numExclusiveOptions = 0;
	numExclusiveOptions += opt->monoecious;
	numExclusiveOptions += opt->cyclical;
	numExclusiveOptions += opt->multiplePopnSizes;
	numExclusiveOptions += opt->sweepSet;
	numExclusiveOptions += opt->fixedMigrantSet;
    if(numExclusiveOptions > 1)
		PERROR("The following options are mutually exclusive: --sweep, --cyclical, --popn-size-file, --monoecious, --fixed-migrant");
    if(opt->fixedMigrantSet && opt->numDemes != 2)
        PERROR("--num-demes must be 2 when using --fixed-migrant");
	if(opt->numGens < 0)
		PERROR("Number of generations (--num-gens) must be greater than zero.");
	// sweep settings:
	if(opt->sweepSet)
	{
		if(!opt->s_set || !opt->h_set || !opt->ts_set)
			PERROR("Must specify selective advantage (--sweep-s), dominance (--sweep-h), and starting generation (--sweep-start-gen) to specify selective sweep.");
	}
    if(opt->fstSet && opt->numDemes < 2)
        PERROR("--num-demes (-D) must be >= 2 for --fst");
    if(opt->fstSet && opt->exactSet)
        PERROR("--fst and --exact are mutually exclusive");
    if(opt->multiplePairs && opt->fstSet)
        PERROR("--num-pairs > 1 and --fst are mutually exclusive");
    if(opt->singleLineageSet && opt->numDemes != 2)
        PERROR("--num-demes (-D) must be 2 when using --single-lineage");
    if(opt->pedsimNumProbGens < 0)
        PERROR("invalid --num-prob-gens (must be > 0)");
	return;
}

