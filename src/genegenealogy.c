#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "node.h"
#include "mutlist.h"
#include "intvector.h"
#include "newicktree.h"
#include "definitions.h"
#include "random.h"
#include "coalpedigree.h"
#include "genegenealogy.h"

void GeneGenealogy_find_terminal_nodes(Node * curNode, Node ** nodeIDs)
{
	if(curNode->desc1 == NULL)
	{
		//nodeIDs[(*curNodeIdx)++] = curNode;
		nodeIDs[curNode->idx] = curNode;
		return;
	}
	GeneGenealogy_find_terminal_nodes(curNode->desc1, nodeIDs);
	GeneGenealogy_find_terminal_nodes(curNode->desc2, nodeIDs);
	return;
}

// mutIdx0 and mutIdx1 are inclusive
void GeneGenealogy_add_mutation(Node * curNode, MutList * ml, int32_t mutIdx0, int32_t mutIdx1)
{
	int32_t i;
	if(curNode->desc1 == NULL)
	{
		for(i = mutIdx0; i <= mutIdx1; i++)
			MutList_add(ml, curNode, i);
		return;
	}
	GeneGenealogy_add_mutation(curNode->desc1, ml, mutIdx0, mutIdx1);
	GeneGenealogy_add_mutation(curNode->desc2, ml, mutIdx0, mutIdx1);
	
	return;
}

void GeneGenealogy_mutate_branch(MutList * ml, Node * curNode, double mu, int32_t * pcurMutIdx)
{
	int32_t n = (curNode->anc)->time - curNode->time;
	unsigned int numMuts = gsl_ran_binomial(n,mu);
	int32_t i;
	int32_t oldIdx;
	ml->S += numMuts;		// Number of segregating sites is incremented.
	curNode->numMuts = numMuts;
	oldIdx = *pcurMutIdx;
	*pcurMutIdx += (int32_t)numMuts;
	if(numMuts > 0)
		GeneGenealogy_add_mutation(curNode, ml, oldIdx, *pcurMutIdx - 1);
	if(curNode->desc1 == NULL)
		return;
	GeneGenealogy_mutate_branch(ml, curNode->desc1, mu, pcurMutIdx);
	GeneGenealogy_mutate_branch(ml, curNode->desc2, mu, pcurMutIdx);
	return;
}

void GeneGenealogy_mutate_branch_conditional(MutList * ml, Node * curNode, double mu, int32_t * pcurMutIdx, int32_t * pmutIndicator)
{
	int32_t n = (curNode->anc)->time - curNode->time;
	unsigned int numMuts = gsl_ran_binomial(n,mu);
	int32_t i;
	int32_t oldIdx;
	ml->S += numMuts;		// Number of segregating sites is incremented.
	curNode->numMuts = numMuts;
	oldIdx = *pcurMutIdx;
	*pcurMutIdx += (int32_t)numMuts;
	if(numMuts > 0)
	{
		GeneGenealogy_add_mutation(curNode, ml, oldIdx, *pcurMutIdx - 1);
		*pmutIndicator = 1;
	}
	if(curNode->desc1 == NULL)
		return;
	GeneGenealogy_mutate_branch_conditional(ml, curNode->desc1, mu, pcurMutIdx, pmutIndicator);
	GeneGenealogy_mutate_branch_conditional(ml, curNode->desc2, mu, pcurMutIdx, pmutIndicator);
	return;
}

// Assumes ml is already initialized.
void GeneGenealogy_simulate_mutations(GeneGenealogy * gt, MutList * ml, double mu)
{
	int32_t right = 0;
	// Assign nodeIDs in ml.
	// Going to have to do it recursively.
	GeneGenealogy_find_terminal_nodes(gt->MRCA, ml->nodeIDs);
	int32_t i;
	int32_t curMutIdx = 0;
	GeneGenealogy_mutate_branch(ml, (gt->MRCA)->desc1, mu, &curMutIdx);
	GeneGenealogy_mutate_branch(ml, (gt->MRCA)->desc2, mu, &curMutIdx);
	return;
}

int32_t GeneGenealogy_simulate_mutations_conditional(GeneGenealogy * gt, MutList * ml, double mu)
{
	int32_t mutIndicator = 0;
	int32_t right = 0;
	// Assign nodeIDs in ml.
	// Going to have to do it recursively.
	GeneGenealogy_find_terminal_nodes(gt->MRCA, ml->nodeIDs);
	int32_t i;
	int32_t curMutIdx = 0;
	GeneGenealogy_mutate_branch_conditional(ml, (gt->MRCA)->desc1, mu, &curMutIdx, &mutIndicator);
	GeneGenealogy_mutate_branch_conditional(ml, (gt->MRCA)->desc2, mu, &curMutIdx, &mutIndicator);
	return mutIndicator;
}

// using numMuts as the mutation Idx.
void GeneGenealogy_mutate_branch_inf_alleles(Node * curNode, double mu, int32_t * pcurMutIdx)
{
	int32_t n = (curNode->anc)->time - curNode->time;
	//fprintf(stderr,"n = %i, mu = %f, pow(1-mu,n) = %f\n", n, mu, pow(1-mu,n));
	if(rbern(pow(1-mu,n)))
		curNode->numMuts = -1;
	else
		curNode->numMuts = (*pcurMutIdx)++;
	if(curNode->desc1 == NULL)
		return;
	GeneGenealogy_mutate_branch_inf_alleles(curNode->desc1, mu, pcurMutIdx);
	GeneGenealogy_mutate_branch_inf_alleles(curNode->desc2, mu, pcurMutIdx);
	return;
}

int32_t GeneGenealogy_get_inf_allele_type(Node * curNode)
{
	if(curNode->anc == NULL)
		return 1;		// this used to be 0, which is missing data in most cases. corresponding change was made to GeneGenealogy_simulate_inf_alleles_mutations.
	if(curNode->numMuts == -1)
		return GeneGenealogy_get_inf_allele_type(curNode->anc);
	else
		return curNode->numMuts;
}

void GeneGenealogy_simulate_inf_alleles_mutations(GeneGenealogy * gt, Node ** terminalNodes, int32_t * alleleTypes, double mu)
{
	int32_t i, curMutIdx = 2;
	// BUT OF COURSE THE TERMINAL NODES DON'T COME IN ORDER!!! this is the bug.
	GeneGenealogy_find_terminal_nodes(gt->MRCA, terminalNodes);
	// going to use numMuts as the mutation Idx.
	GeneGenealogy_mutate_branch_inf_alleles((gt->MRCA)->desc1, mu, &curMutIdx);
	GeneGenealogy_mutate_branch_inf_alleles((gt->MRCA)->desc2, mu, &curMutIdx);
	// Now, need to fill out the alleleTypes data.
	for(i = 0; i < gt->sampleSize; i++)
		alleleTypes[terminalNodes[i]->idx] = GeneGenealogy_get_inf_allele_type(terminalNodes[i]);
	return;
}

void GeneGenealogy_init(GeneGenealogy * tree, int32_t numIndiv)
{
	//int32_t sz = 2*numIndiv - 1;
	int32_t sz = 2*numIndiv;
	tree->nodeList = (Node *)malloc(sizeof(Node) * sz);
	CHECKPOINTER(tree->nodeList);
	tree->numNodes = sz;
	tree->currentNode = 0;
	tree->sampleSize = numIndiv;
	tree->tMRCA = -1;
	tree->MRCA = NULL;
	tree->nodeListContainer = NULL;
	return;
}

void GeneGenealogy_reset(GeneGenealogy * tree)
{
	tree->numNodes = 0;
	tree->currentNode = 0;
	tree->sampleSize = 0;
	tree->tMRCA = -1;
	tree->MRCA = NULL;
	return;
}

void GeneGenealogy_free(GeneGenealogy * tree)
{
	free(tree->nodeList);
	NodeList_free(tree->nodeListContainer);
	free(tree->nodeListContainer);
	tree->numNodes = 0;
	tree->currentNode = 0;
	tree->sampleSize = 0;
	tree->tMRCA = -1;
	tree->MRCA = NULL;
	return;
}

Node * GeneGenealogy_get_Node(GeneGenealogy * tree)
{
	Node * node = NULL;
	if(tree->currentNode + 1 != tree->numNodes)
	{
		node = &(tree->nodeList[tree->currentNode]);
		(tree->currentNode)++;
		node->anc = NULL;
	}
	else
	{
		PERROR("Accessing non-existent Node address in GeneGenealogy_get_Node().");
	}
	return(node);
}

void GeneGenealogy_sim_gene_genealogy(CoalPedigree * ped, GeneGenealogy * tree, int32_t * indiv, int32_t numIndiv, int32_t loopTime)
{
	NodeList * list = (NodeList *)malloc(sizeof(NodeList));
	CHECKPOINTER(list);
	GeneGenealogy_init(tree, numIndiv);
	NodeList_init(list, tree, indiv, numIndiv);
	tree->nodeListContainer = list;
	Node * cnode, * onode;
	int32_t i,j;
	int32_t t = 0;
	for(t = 0; t < ped->numGenerations; t++)
	{
		// First, advance them.
		for(i = 0; i < list->numNodes; i++)
		{
			cnode = list->nodeList[i];
			cnode->indiv = ped->relationships[t][cnode->indiv][cnode->parent];
			cnode->parent = rbern(0.5);
		}
		// Then, check for coalescence.
		for(i = 0; i < list->numNodes-1; i++)
		{
			cnode = list->nodeList[i];
			for(j = i+1; j < list->numNodes; j++)
			{
				onode = list->nodeList[j];
				// if they're the same chromosome...
				if(cnode->indiv == onode->indiv && cnode->parent == onode->parent)
				{
					list->nodeList[i] = cnode = NodeList_coalesce_Nodes(cnode, onode, t, list, tree);
					j--;
				}
			}
		}
		if(list->numNodes == 1)
		{
			if(t == 0)
				fprintf(stderr,"Coalescence at t=0!\n");
			tree->MRCA = list->nodeList[0];
			return;
		}
	}
	int32_t tp;
	while(1)
	{
		for(tp = loopTime; tp < ped->numGenerations; tp++)
		{
			// First, advance them.
			for(i = 0; i < list->numNodes; i++)
			{
				cnode = list->nodeList[i];
				cnode->indiv = ped->relationships[tp][cnode->indiv][cnode->parent];
				cnode->parent = rbern(0.5);
			}
			// Then, check for coalescence.
			for(i = 0; i < list->numNodes-1; i++)
			{
				cnode = list->nodeList[i];
				for(j = i+1; j < list->numNodes; j++)
				{
					onode = list->nodeList[j];
					// if they're the same chromosome...
					if(cnode->indiv == onode->indiv && cnode->parent == onode->parent)
					{
						list->nodeList[i] = cnode = NodeList_coalesce_Nodes(cnode, onode, t, list, tree);
						// ^ Yes, t, not tp.
						j--;
					}
				}
			}
			if(list->numNodes == 1)
			{
				tree->tMRCA = t;
				tree->MRCA = list->nodeList[0];
				return;
			}
			t++;
		}
	}
	GeneGenealogy_free(tree);
	NodeList_free(list);
	return;
}

void GeneGenealogy_sim_gene_genealogy_head(CoalPedigree * ped, GeneGenealogy * tree, int32_t * indiv, int32_t numIndiv, int32_t finalGen)
{
	NodeList * list = (NodeList *)malloc(sizeof(NodeList));
	CHECKPOINTER(list);
	GeneGenealogy_init(tree, numIndiv);
	NodeList_init(list, tree, indiv, numIndiv);
	Node * cnode, * onode;
	tree->nodeListContainer = list;
	int32_t i,j;
	int32_t t = 0;
	//REPORTI(list->numNodes);
	for(t = 0; t < finalGen; t++)
	{
		// First, advance them.
		for(i = 0; i < list->numNodes; i++)
		{
			cnode = list->nodeList[i];
			cnode->indiv = ped->relationships[t][cnode->indiv][cnode->parent];
			cnode->parent = rbern(0.5);
		}
		// Then, check for coalescence.
		for(i = 0; i < list->numNodes-1; i++)
		{
			cnode = list->nodeList[i];
			for(j = i+1; j < list->numNodes; j++)
			{
				onode = list->nodeList[j];
				// if they're the same chromosome...
				if(cnode->indiv == onode->indiv && cnode->parent == onode->parent)
				{
					list->nodeList[i] = cnode = NodeList_coalesce_Nodes(cnode, onode, t, list, tree);
					j--;
				}
			}
		}
		if(list->numNodes == 1)
		{
			tree->MRCA = list->nodeList[0];
			//REPORTI(list->numNodes);
			return;
		}
	}
	//REPORTI(list->numNodes);
	return;
}

void GeneGenealogy_sim_diploid_gene_genealogy(CoalPedigree * ped, GeneGenealogy * tree, int32_t * indiv, int32_t numIndiv, int32_t loopTime)
{
	int32_t i, j, t = 0;
	NodeList * list = (NodeList *)malloc(sizeof(NodeList));
	CHECKPOINTER(list);
	GeneGenealogy_init(tree, 2*numIndiv);
	NodeList_init_diploid(list, tree, ped, indiv, numIndiv);
	tree->nodeListContainer = list;
	Node * cnode, * onode;
	for(t = 1; t < ped->numGenerations; t++)		// note t = 1; we already advanced them once in the diploid initialization.
	{
		// First, check for coalescence.
		for(i = 0; i < list->numNodes-1; i++)
		{
			cnode = list->nodeList[i];
			for(j = i+1; j < list->numNodes; j++)
			{
				onode = list->nodeList[j];
				// if they're the same chromosome...
				if(cnode->indiv == onode->indiv && cnode->parent == onode->parent)
				{
					list->nodeList[i] = cnode = NodeList_coalesce_Nodes(cnode, onode, t, list, tree);
					j--;
				}
			}
		}
		if(list->numNodes == 1)
		{
			if(t == 0)
				PERROR("Coalescence at t=0!");
			tree->MRCA = list->nodeList[0];
			return;
		}
		// Then, advance them.
		for(i = 0; i < list->numNodes; i++)
		{
			cnode = list->nodeList[i];
			cnode->indiv = ped->relationships[t][cnode->indiv][cnode->parent];
			cnode->parent = rbern(0.5);
		}
	}
	int32_t tp;
	while(1)
	{
		for(tp = loopTime; tp < ped->numGenerations; tp++)
		{
			// First, check for coalescence.
			for(i = 0; i < list->numNodes-1; i++)
			{
				cnode = list->nodeList[i];
				for(j = i+1; j < list->numNodes; j++)
				{
					onode = list->nodeList[j];
					// if they're the same chromosome...
					if(cnode->indiv == onode->indiv && cnode->parent == onode->parent)
					{
						list->nodeList[i] = cnode = NodeList_coalesce_Nodes(cnode, onode, t, list, tree);
						// ^ Yes, t, not tp.
						j--;
					}
				}
			}
			if(list->numNodes == 1)
			{
				tree->tMRCA = t;
				tree->MRCA = list->nodeList[0];
				return;
			}
			t++;
			// Then, advance them.
			for(i = 0; i < list->numNodes; i++)
			{
				cnode = list->nodeList[i];
				cnode->indiv = ped->relationships[tp][cnode->indiv][cnode->parent];
				cnode->parent = rbern(0.5);
			}
		}
	}
}

void GeneGenealogy_sim_diploid_gene_genealogy_debug(CoalPedigree * ped, GeneGenealogy * tree, int32_t * indiv, int32_t numIndiv, int32_t loopTime)
{
	int32_t i, j, t = 0;
	NodeList * list = (NodeList *)malloc(sizeof(NodeList));
	CHECKPOINTER(list);
	GeneGenealogy_init(tree, 2*numIndiv);
	NodeList_init_diploid(list, tree, ped, indiv, numIndiv);
	tree->nodeListContainer = list;
	Node * cnode, * onode;
	for(t = 1; t < ped->numGenerations; t++)		// note t = 1; we already advanced them once in the diploid initialization.
	{
		// First, check for coalescence.
		for(i = 0; i < list->numNodes-1; i++)
		{
			cnode = list->nodeList[i];
			for(j = i+1; j < list->numNodes; j++)
			{
				onode = list->nodeList[j];
				// if they're the same chromosome...
				if(cnode->indiv == onode->indiv && cnode->parent == onode->parent)
				{
					printf("coal %p - %p in gen %i\n", cnode, onode, t);
					list->nodeList[i] = cnode = NodeList_coalesce_Nodes(cnode, onode, t, list, tree);
					j--;
				}
			}
		}
		if(list->numNodes == 1)
		{
			if(t == 0)
				PERROR("Coalescence at t=0!");
			tree->MRCA = list->nodeList[0];
			return;
		}
		// Then, advance them.
		for(i = 0; i < list->numNodes; i++)
		{
			cnode = list->nodeList[i];
			cnode->indiv = ped->relationships[t][cnode->indiv][cnode->parent];
			cnode->parent = rbern(0.5);
		}
	}
	int32_t tp;
	while(1)
	{
		for(tp = loopTime; tp < ped->numGenerations; tp++)
		{
			// First, check for coalescence.
			for(i = 0; i < list->numNodes-1; i++)
			{
				cnode = list->nodeList[i];
				for(j = i+1; j < list->numNodes; j++)
				{
					onode = list->nodeList[j];
					// if they're the same chromosome...
					if(cnode->indiv == onode->indiv && cnode->parent == onode->parent)
					{
						list->nodeList[i] = cnode = NodeList_coalesce_Nodes(cnode, onode, t, list, tree);
						// ^ Yes, t, not tp.
						j--;
					}
				}
			}
			if(list->numNodes == 1)
			{
				tree->tMRCA = t;
				tree->MRCA = list->nodeList[0];
				return;
			}
			t++;
			// Then, advance them.
			for(i = 0; i < list->numNodes; i++)
			{
				cnode = list->nodeList[i];
				cnode->indiv = ped->relationships[tp][cnode->indiv][cnode->parent];
				cnode->parent = rbern(0.5);
			}
		}
	}
}

void GeneGenealogy_sim_gene_genealogy_diploid_head(CoalPedigree * ped, GeneGenealogy * tree, int32_t * indiv, int32_t numIndiv, int32_t finalGen)
{
	NodeList * list = (NodeList *)malloc(sizeof(NodeList));
	CHECKPOINTER(list);
	GeneGenealogy_init(tree, 2*numIndiv);
	NodeList_init_diploid(list, tree, ped, indiv, numIndiv);
	Node * cnode, * onode;
	tree->nodeListContainer = list;
	int32_t i,j;
	int32_t t = 0;
	//REPORTI(list->numNodes);
	for(t = 1; t < finalGen; t++)		// note t = 1; we already advanced them once in the diploid initialization.
	{
		// First, check for coalescence.
		for(i = 0; i < list->numNodes-1; i++)
		{
			cnode = list->nodeList[i];
			for(j = i+1; j < list->numNodes; j++)
			{
				onode = list->nodeList[j];
				// if they're the same chromosome...
				if(cnode->indiv == onode->indiv && cnode->parent == onode->parent)
				{
					list->nodeList[i] = cnode = NodeList_coalesce_Nodes(cnode, onode, t, list, tree);
					j--;
				}
			}
		}
        // BUGFIX previously, would check here to see whether numNodes == 1 and
        // return if so.  now allowing single lineage to traverse pedigree to
        // get final individual, allow correct deme-checking
		for(i = 0; i < list->numNodes; i++)
		{
			cnode = list->nodeList[i];
			cnode->indiv = ped->relationships[t][cnode->indiv][cnode->parent];
			cnode->parent = rbern(0.5);
		}
	}
	//REPORTI(list->numNodes);
	return;
}

inline int32_t GeneGenealogy_check_done(GeneGenealogy * tree)
{
	// tree->numNodes is 2*n-1, the total number of nodes needed.
	// You're done iff 2*n-1 nodes have been accessed.
	if(tree->currentNode == tree->numNodes)
		return(1);
	else
		return(0);
}

// Assumes tree has been constructed already 
int32_t GeneGenealogy_calc_sackin_index(GeneGenealogy * tree)
{
	Node ** leaves = (Node **)malloc(sizeof(Node *) * (size_t)tree->sampleSize);
	CHECKPOINTER(leaves);
	int32_t i, curDist, sackin = 0;
	Node * curNode;
	GeneGenealogy_find_terminal_nodes(tree->MRCA, leaves);
	for(i = 0; i < tree->sampleSize; i++)
	{
		curNode = leaves[i];
		while(curNode->anc != NULL)
		{
			curNode = curNode->anc;
			sackin++;
		}
	}
	free(leaves);
	return(sackin);
}

int32_t GeneGenealogy_calc_colless_index(GeneGenealogy * tree)
{
	int32_t colless = 0;
	Node_calc_colless_index_recursive(tree->MRCA, &colless);
	return colless;
}

Node * GeneGenealogy_get_MRCA(GeneGenealogy * tree)
{
	Node * cnode = &(tree->nodeList[0]);
	while(cnode->anc != NULL)
		cnode = cnode->anc;
	return(cnode);
}

void GeneGenealogy_check_resolved_recursive_(Node * curNode, int32_t * pResolved)
{
	if(curNode->desc1 == NULL || *pResolved == 0)
		return;
	// if there are mutations on the branch above the node...
	if(curNode->numMuts > 0)
	{	
		GeneGenealogy_check_resolved_recursive_(curNode->desc1, pResolved);
		GeneGenealogy_check_resolved_recursive_(curNode->desc2, pResolved);
	}
	else
	{
		//printf("no mutations...\n");
		*pResolved = 0;
	}
	return;
}

int32_t GeneGenealogy_check_resolved(GeneGenealogy * tree)
{
	Node * mrca = GeneGenealogy_get_MRCA(tree);
	int32_t resolved = 1;
	GeneGenealogy_check_resolved_recursive_(mrca->desc1, &resolved);
	GeneGenealogy_check_resolved_recursive_(mrca->desc2, &resolved);
	return resolved;
}

void GeneGenealogy_find_next_node(Node * curNode, IntVector * iv)
{
	// Go down the left side looking for new nodes, then the right side.
	if(curNode->desc1 == NULL)
		return;
	if((curNode->desc1)->numMuts > 0 && (curNode->desc1)->desc1 != NULL)
		GeneGenealogy_add_node_to_distn(curNode->desc1, iv);
	else
		GeneGenealogy_find_next_node(curNode->desc1, iv);
	if((curNode->desc2)->numMuts > 0 && (curNode->desc2)->desc1 != NULL)
		GeneGenealogy_add_node_to_distn(curNode->desc2, iv);
	else
		GeneGenealogy_find_next_node(curNode->desc2, iv);
	return;
}


void GeneGenealogy_add_node_to_distn(Node * curNode, IntVector * iv)
{
	int32_t curSize = 0;
	Node_add_to_size(curNode, &curSize);
	IntVector_add(iv, curSize);
	GeneGenealogy_find_next_node(curNode, iv);
	return;
}

void GeneGenealogy_get_node_size_distn(GeneGenealogy * tree, IntVector * iv)
{
	int32_t curSize;
	Node * mrca = GeneGenealogy_get_MRCA(tree);
	GeneGenealogy_add_node_to_distn(mrca, iv);
	return;
}

// note curIdx is just a pointer pointing to the current value of the index

// Moves up a tree and returns the address of the MRCA.
int32_t GeneGenealogy_calc_new_index(GeneGenealogy * tree)
{
	Node ** leaves = (Node **)malloc(sizeof(Node *) * (size_t)tree->sampleSize);
	CHECKPOINTER(leaves);
	// distances will hold the distance of each node from the root
	int32_t * distances = (int32_t *)malloc(sizeof(int32_t) * (size_t)tree->sampleSize);
	int32_t i, j, curDist;
	Node * curNode;
	GeneGenealogy_find_terminal_nodes(tree->MRCA, leaves);
	for(i = 0; i < tree->sampleSize; i++)
	{
		curNode = leaves[i];
		distances[i] = 0;
		while(curNode->anc != NULL)
		{
			curNode = curNode->anc;
			distances[i]++;
		}
	}
	twoints curDists;
	Node * curNodes[2];
	int32_t currentlyGreater;
	Node * ultAnc = GeneGenealogy_get_MRCA(tree);
	int32_t depth = 0;
	int32_t newIdx = 0;
	for(i = 0; i < tree->sampleSize-1; i++)
	{
		for(j = i + 1; j < tree->sampleSize; j++)
		{
			curNodes[0] = leaves[i];
			curDists[0] = distances[i];
			curNodes[1] = leaves[j];
			curDists[1] = distances[j];
			//while(!(curNodes[0] == ultAnc && curNodes[1] == ultAnc))
			while(1)
			{
				if(curDists[0] == curDists[1] && curNodes[0] == curNodes[1])
				{
					depth = 0;
					while(curNodes[0]->anc != NULL)
					{
						curNodes[0] = curNodes[0]->anc;
						depth++;
					}
					break;
				}
				currentlyGreater = (curDists[0] <= curDists[1]);
				curNodes[currentlyGreater] = curNodes[currentlyGreater]->anc;
				curDists[currentlyGreater]--;
			}
			newIdx += depth;
		}
	}
	free(distances);
	free(leaves);
	return(newIdx);
}


void GeneGenealogy_print_tree_recursive_(GeneGenealogy * tree, NewickTree * nt, Node * curAnc, int32_t mutations)
{
	char staging[50];
	int32_t stagingLen;
	if(curAnc->desc1 == NULL && curAnc->desc2 == NULL)		// End of recursion...
	{
		//sprintf(staging, ":%i", tree->tMRCA - curAnc->time);
		if(mutations)
			sprintf(staging, ":%i (m=%i)",(curAnc->anc)->time - curAnc->time, curAnc->numMuts);
		else
			sprintf(staging, ":%i",(curAnc->anc)->time - curAnc->time);
		stagingLen = index(staging,'\0') - staging;
		NewickTree_add(nt, staging, stagingLen);
		return;
	}
	NewickTree_add(nt,"(", 1);
	GeneGenealogy_print_tree_recursive_(tree, nt, curAnc->desc1, mutations);
	NewickTree_add(nt,",", 1);
	GeneGenealogy_print_tree_recursive_(tree, nt, curAnc->desc2, mutations);
	NewickTree_add(nt,")", 1);
	//sprintf(staging, ":%i", tree->tMRCA - curAnc->time);
	if(curAnc->anc != NULL)
	{
		if(mutations)
			sprintf(staging, ":%i (m = %i)",(curAnc->anc)->time - curAnc->time, curAnc->numMuts);
		else
			sprintf(staging, ":%i",(curAnc->anc)->time - curAnc->time);
		stagingLen = index(staging,'\0') - staging;
		NewickTree_add(nt, staging, stagingLen);
	}
	return;
}

void GeneGenealogy_print_tree(GeneGenealogy * tree, NewickTree * nt, FILE * fout, int32_t mutations)
{
	GeneGenealogy_print_tree_recursive_(tree, nt, GeneGenealogy_get_MRCA(tree), mutations);
	NewickTree_add(nt, ";", 1);
	fprintf(fout, "%s\n", nt->line);
	return;
}

