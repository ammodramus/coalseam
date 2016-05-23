#ifndef GENEGENEALOGY_HEADER
#define GENEGENEALOGY_HEADER

#include <stdint.h>
#include "node.h"
#include "mutlist.h"
#include "intvector.h"
#include "newicktree.h"
#include "definitions.h"
#include "coalpedigree.h"

typedef struct genegenealogy 
{
	Node * nodeList;
	NodeList * nodeListContainer;
	int32_t numNodes;
	int32_t currentNode;
	int32_t sampleSize;
	int32_t tMRCA;
	Node * MRCA;
} GeneGenealogy;

void GeneGenealogy_find_terminal_nodes(Node * curNode, Node ** nodeIDs);
void GeneGenealogy_add_mutation(Node * curNode, MutList * ml, int32_t mutIdx0, int32_t mutIdx1);
void GeneGenealogy_mutate_branch(MutList * ml, Node * curNode, double mu, int32_t * pcurMutIdx);
void GeneGenealogy_mutate_branch_conditional(MutList * ml, Node * curNode, double mu, int32_t * pcurMutIdx, int32_t * pmutIndicator);
void GeneGenealogy_simulate_mutations(GeneGenealogy * gt, MutList * ml, double mu);
int32_t GeneGenealogy_simulate_mutations_conditional(GeneGenealogy * gt, MutList * ml, double mu);
void GeneGenealogy_mutate_branch_inf_alleles(Node * curNode, double mu, int32_t * pcurMutIdx);
int32_t GeneGenealogy_get_inf_allele_type(Node * curNode);
void GeneGenealogy_simulate_inf_alleles_mutations(GeneGenealogy * gt, Node ** terminalNodes, int32_t * alleleTypes, double mu);
void GeneGenealogy_init(GeneGenealogy * tree, int32_t numIndiv);
void GeneGenealogy_reset(GeneGenealogy * tree);
void GeneGenealogy_free(GeneGenealogy * tree);
Node * GeneGenealogy_get_Node(GeneGenealogy * tree);
inline int32_t GeneGenealogy_check_done(GeneGenealogy * tree);
int32_t GeneGenealogy_calc_sackin_index(GeneGenealogy * tree);
int32_t GeneGenealogy_calc_colless_index(GeneGenealogy * tree);
Node * GeneGenealogy_get_MRCA(GeneGenealogy * tree);
void GeneGenealogy_check_resolved_recursive_(Node * curNode, int32_t * pResolved);
int32_t GeneGenealogy_check_resolved(GeneGenealogy * tree);
void GeneGenealogy_find_next_node(Node * curNode, IntVector * iv);
void GeneGenealogy_add_node_to_distn(Node * curNode, IntVector * iv);
void GeneGenealogy_get_node_size_distn(GeneGenealogy * tree, IntVector * iv);
int32_t GeneGenealogy_calc_new_index(GeneGenealogy * tree);
void GeneGenealogy_print_tree_recursive_(GeneGenealogy * tree, NewickTree * nt, Node * curAnc, int32_t mutations);
void GeneGenealogy_print_tree(GeneGenealogy * tree, NewickTree * nt, FILE * fout, int32_t mutations);
void GeneGenealogy_sim_gene_genealogy(CoalPedigree * ped, GeneGenealogy * tree, int32_t * indiv, int32_t numIndiv, int32_t loopTime);
void GeneGenealogy_sim_gene_genealogy_head(CoalPedigree * ped, GeneGenealogy * tree, int32_t * indiv, int32_t numIndiv, int32_t finalGen);
void GeneGenealogy_sim_diploid_gene_genealogy(CoalPedigree * ped, GeneGenealogy * tree, int32_t * indiv, int32_t numIndiv, int32_t loopTime);

#endif
