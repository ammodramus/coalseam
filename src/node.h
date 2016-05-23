#ifndef NODE_HEADER
#define NODE_HEADER

#include <stdint.h>

typedef struct genegenealogy GeneGenealogy;
typedef struct coalpedigree CoalPedigree;

typedef struct node
{
	int32_t time;
	int32_t indiv;
	int32_t parent;
	int32_t numMuts;
	int32_t idx;
	struct node * desc1;
	struct node * desc2;
	struct node * anc;
} Node;

typedef struct nodelist
{
	Node ** nodeList;
	int32_t numNodes;
	int32_t sizeNowNodeList;
} NodeList;

void NodeList_init(NodeList * list, GeneGenealogy * tree, int32_t * indiv, int32_t numIndiv);
void NodeList_init_diploid(NodeList * list, GeneGenealogy * tree, CoalPedigree * ped, int32_t * indiv, int32_t numIndiv);
inline void NodeList_free(NodeList * nl);
void NodeList_remove_Node(Node * node, NodeList * list);
Node * NodeList_coalesce_Nodes(Node * node1, Node * node2, int32_t time, NodeList * list, GeneGenealogy * tree);
void Node_calc_num_descendents_recursive(Node * curNode, int32_t * curCount);
int32_t Node_calc_num_descendents(Node * curNode);
void Node_calc_colless_index_recursive(Node * curNode, int32_t * curIdx);
void Node_add_to_size(Node * curNode, int32_t * curSize);

#endif
