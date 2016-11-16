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

#ifndef OPTIONS_HEADER
#define OPTIONS_HEADER

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include "definitions.h"

typedef struct options_
{
	double theta;     // 4*N*mu, where N is the deme (subpopulation) size and mu is the per-generation mutation probability
	double migRate;     // 4*N*m, where N is the deme (subpopulation) size and m is the per-generation migration probability
	int32_t numDemes;   // number of demes or subpopulations
	int32_t demeSize;   // number of diploid individuals per deme
	int32_t sampleSize; // number of individuals sampled *per deme*
	int32_t monoecious; // indicator of whether the population is monoecious
	int32_t cyclical;  // indicator of whether the population is cyclical (repeat same relationships every generation)
	int32_t multiplePopnSizes;  // flag indicating whether there are multiple population sizes
	char popnSizeFile[200];    // filename for the multiple population sizes, if any
	int32_t numLoci;    // number of loci to simulate
	int32_t genealogiesSet;   // flag indicating whether genealogies are to be added to output
	int32_t msSet;
	int32_t genepopSet;
	int32_t infAllelesSet;
	char geneFilename[200];
	int32_t sweepSet;
	int32_t s_set;
	double sweep_s;
	int32_t h_set;
	double sweep_h;
	int32_t ts_set;
	int32_t sweep_ts;
	int32_t requireFixation;
	int32_t importSet;
	char importFilename[200];
	int32_t exportSet;
	int32_t numGens;
	char exportFilename[200];
	int32_t loopGen;
	char rearrFilename[200];
	int32_t rearrSet;
	int32_t rearrGen;
	int32_t rearrReps;
    int32_t fixedIBDset;
    int32_t fixedIBDgen;
    int32_t fixedIBDindiv;
    int32_t fixedMigrantSet;
    int32_t fixedMigrantGen;
    int32_t fixedMigrantIndiv;
    int32_t * indivs;
    int32_t sharedAncSet;
    int32_t sharedAncGen;
    int32_t sharedAncIndivs[2];
    // just for pedsim
    int32_t exactSet;
    int32_t numCoalTimes;
    int32_t withinSet;
    int32_t fstSet;
    int32_t numPairs;
    int32_t multiplePairs;
    int32_t singleLineageSet; 
    int32_t numSingleLineages;
    int32_t pedsimNumProbGens;
} Options;

void Options_parse_options(int argc, char ** argv, Options * opt);

#endif

