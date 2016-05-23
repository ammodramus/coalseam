#ifndef NEWICKTREE_HEADER
#define NEWICKTREE_HEADER

#include <stdint.h>

typedef struct newicktree
{
	char * line;
	int32_t curIdx;
	int32_t curSize;
} NewickTree;

void NewickTree_init(NewickTree * nt);
void NewickTree_add(NewickTree * nt, char * toadd, int32_t toaddLength);
void NewickTree_reset(NewickTree * nt);
void NewickTree_free(NewickTree * nt);

#endif
