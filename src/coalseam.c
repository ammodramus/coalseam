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

/* outputs results in ms(-like) format
 *
 * muts         	vector of MutList structs, of length numLoci
 * genealogies   	vector of GeneGenealogy structs, of length numLoci (ignored if printGenealogies == 0)
 * numLoci       	number of loci to output
 * printGenealogies indicates whether to print genealogies on line before "segsites: "
 * filename      	output filename
 *
 * (later, add command? positions (which are Uniform(0,1)?))
 * ... positions may not make sense with Mendelian loci
 */
void print_ms_format(MutList * muts, GeneGenealogy * genealogies, int32_t numLoci, int32_t printGenealogies, FILE * fout)
{
	int32_t i;
	NewickTree nt;
	if(printGenealogies)
		NewickTree_init(&nt);
	for(i = 0; i < numLoci-1; i++)
	{
		fprintf(fout, "//\n");
		if(printGenealogies)
		{
			GeneGenealogy_print_tree(&(genealogies[i]), &nt, fout, 0);
			NewickTree_reset(&nt);
		}
		fprintf(fout, "segsites: %i\n", muts[i].S);
		if(muts[i].S > 0)
			MutList_print_ms_format(&(muts[i]), fout);
		fprintf(fout, "\n");
	}
	fprintf(fout, "//\n");
	if(printGenealogies)
	{
		GeneGenealogy_print_tree(&(genealogies[i]), &nt, fout, 0);
		NewickTree_free(&nt);
	}
	fprintf(fout, "segsites: %i\n", muts[numLoci-1].S);
	MutList_print_ms_format(&(muts[numLoci-1]), fout);
	fclose(fout);
	return;
}

void get_init_indivs(Options * opt, int32_t * initIndivs)
{
	int32_t i, j;
	int32_t * demeIndivs = (int32_t *)malloc(sizeof(int32_t) * (size_t)opt->sampleSize);
	CHECKPOINTER(demeIndivs);
	int32_t * demeIndivsBucket = (int32_t *)malloc(sizeof(int32_t) * (size_t)opt->demeSize);
	CHECKPOINTER(demeIndivsBucket);
	for(i = 0; i < opt->numDemes; i++)
	{
		equal_sample_noreplace(opt->sampleSize, opt->demeSize, demeIndivs, demeIndivsBucket);
		for(j = 0; j < opt->sampleSize; j++)
			initIndivs[j+opt->sampleSize*i] = demeIndivs[j] + opt->demeSize * i;
	}
	free(demeIndivsBucket);
	free(demeIndivs);
	return;
}

int32_t main(int argc, char * argv[])
{
	CoalPedigree ped;
	GeneGenealogy tree;
	Options opt;

	int32_t i, j, numIndivs;
	double mu;

    // delete if not needed
	//char optString[] = "D:N:n:M:t:k:r:w:b:saSg:c";

	randseed();

	Options_parse_options(argc, argv, &opt);

	numIndivs = opt.sampleSize * opt.numDemes;

	// get initial individuals
	int32_t * initIndivs = (int32_t *)malloc(sizeof(int32_t) * (size_t)numIndivs);
	CHECKPOINTER(initIndivs);

	get_init_indivs(&opt, initIndivs);

    opt.indivs = initIndivs;

	int32_t * indivs = (int32_t *)malloc(sizeof(int32_t) * (size_t)numIndivs);
	CHECKPOINTER(indivs);


	CoalPedigree_init(&ped, &opt);
	CoalPedigree_shuffle(&ped, &opt);

	if(opt.exportSet)
		CoalPedigree_export_pedigree(&ped, opt.exportFilename);
	// rearrangement probabilities
	if(opt.rearrSet)
	{
		FILE * rearrOut = fopen(opt.rearrFilename, "w");
		if(!rearrOut)
			PERROR("Could not open file to write rearrangement probability.");
        if(opt.numDemes == 1 || opt.numDemes == 2)
            Rearrangement_calculate_reconfiguration_probs_diploid(&ped, initIndivs, numIndivs, opt.rearrGen, rearrOut);
        else
            PERROR("Rearrangement output is supported for only 1 or 2 demes. (-D 1 or -D 2)");
		fclose(rearrOut);
	}
	
    // print ms-data
    MutList * muts = malloc(sizeof(MutList) * opt.numLoci);
    CHECKPOINTER(muts);
    GeneGenealogy * genealogies = malloc(sizeof(GeneGenealogy) * opt.numLoci);
    CHECKPOINTER(genealogies);
    for(i = 0; i < opt.numLoci; i++)
        MutList_init(&(muts[i]), 2*numIndivs);

    mu = opt.theta/(4.0*(double)opt.demeSize);

    for(i = 0; i < opt.numLoci; i++)
    {
        for(j = 0; j < numIndivs; j++)
            indivs[j] = initIndivs[j];
        GeneGenealogy_sim_diploid_gene_genealogy(&ped, &(genealogies[i]), indivs, numIndivs, 0);
        GeneGenealogy_simulate_mutations(&(genealogies[i]), &(muts[i]), mu);
    }
    print_ms_format(muts, genealogies, opt.numLoci, opt.genealogiesSet, stdout);
    for(i = 0; i < opt.numLoci; i++)
    {
        GeneGenealogy_free(&(genealogies[i]));
        MutList_free(&(muts[i]));
    }
    free(muts);
    free(genealogies);


	CoalPedigree_free(&ped);
	free(initIndivs);
	free(indivs);
	return 0;
}

