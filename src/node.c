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
#include <stdint.h>
#include "genegenealogy.h"
#include "coalpedigree.h"
#include "random.h"
#include "node.h"
#include "definitions.h"

void NodeList_init(NodeList * list, GeneGenealogy * tree, int32_t * indiv, int32_t numIndiv)
{
	list->nodeList = (Node **)malloc(sizeof(Node *) * (size_t)numIndiv);
	CHECKPOINTER(list->nodeList);
	int i;
	for(i = 0; i < numIndiv; i++)
	{
		list->nodeList[i] = GeneGenealogy_get_Node(tree);
		CHECKPOINTER(list->nodeList[i]);
		list->nodeList[i]->desc1 = NULL;
		list->nodeList[i]->desc2 = NULL;
		list->nodeList[i]->anc = NULL;
		list->nodeList[i]->time = 0;
		list->nodeList[i]->indiv = indiv[i];
		list->nodeList[i]->parent = rbern(0.5);
		list->nodeList[i]->numMuts = 0;
		list->nodeList[i]->idx = i;
	}
	list->numNodes = numIndiv;
	return;
}

	//NodeList_init_diploid(&list, tree, ped, indiv, numIndiv);
void NodeList_init_diploid(NodeList * list, GeneGenealogy * tree, CoalPedigree * ped, int32_t * indiv, int32_t numIndiv)
{
	list->nodeList = (Node **)malloc(sizeof(Node *) * 2*(size_t)numIndiv);
	CHECKPOINTER(list->nodeList);
	int i,j;
	for(i = 0; i < numIndiv; i++)
	{
		list->nodeList[2*i] = GeneGenealogy_get_Node(tree);
		CHECKPOINTER(list->nodeList[2*i]);
		list->nodeList[2*i]->desc1 = NULL;
		list->nodeList[2*i]->desc2 = NULL;
		list->nodeList[2*i]->anc = NULL;
		list->nodeList[2*i]->time = 0;
		list->nodeList[2*i]->indiv = ped->relationships[0][indiv[i]][0];
		list->nodeList[2*i]->parent = rbern(0.5);
		list->nodeList[2*i]->numMuts = -2;
		list->nodeList[2*i]->idx = 2*i;
		list->nodeList[2*i+1] = GeneGenealogy_get_Node(tree);
		CHECKPOINTER(list->nodeList[2*i+1]);
		list->nodeList[2*i+1]->desc1 = NULL;
		list->nodeList[2*i+1]->desc2 = NULL;
		list->nodeList[2*i+1]->anc = NULL;
		list->nodeList[2*i+1]->time = 0;
		list->nodeList[2*i+1]->indiv = ped->relationships[0][indiv[i]][1];
		list->nodeList[2*i+1]->parent = rbern(0.5);
		list->nodeList[2*i+1]->numMuts = -2;
		list->nodeList[2*i+1]->idx = 2*i+1;
	}
	list->numNodes = 2*numIndiv;
	return;
}

inline void NodeList_free(NodeList * nl)
{
	free(nl->nodeList);
	return;
}

void NodeList_remove_Node(Node * node, NodeList * list)
{
	int32_t i,j;
	for(i = 0; i < list->numNodes; i++)
	{
		if(list->nodeList[i] == node)
		{
			for(j = i+1; j < list->numNodes; j++)
				list->nodeList[j-1] = list->nodeList[j];
			break;
		}
	}
	(list->numNodes)--;
	return;
}

Node * NodeList_coalesce_Nodes(Node * node1, Node * node2, int32_t time, NodeList * list, GeneGenealogy * tree)
{
	Node * anc = GeneGenealogy_get_Node(tree);
	CHECKPOINTER(anc);
	node1->anc = node2->anc = anc;
	anc->desc1 = node1;
	anc->desc2 = node2;
	anc->indiv = node1->indiv;		// note node1->indiv == node2->indiv.
	anc->parent = node1->parent;
	anc->time = time;
	// Always remove the second node, since it's the onode being compared against cnode.
	NodeList_remove_Node(node2, list);		// Be careful with changing the length of lists while you're looping through them!
	return(anc);
}

void Node_calc_num_descendents_recursive(Node * curNode, int32_t * curCount)
{
	int32_t numDescendents;
	if(curNode->desc1 == NULL)
	{
		(*curCount)++;
		return;
	}
	Node_calc_num_descendents_recursive(curNode->desc1,curCount);
	Node_calc_num_descendents_recursive(curNode->desc2,curCount);
	return;
}

int32_t Node_calc_num_descendents(Node * curNode)
{
	int32_t numDescendents = 0;
	Node_calc_num_descendents_recursive(curNode, &numDescendents);
	return numDescendents;
}

void Node_calc_colless_index_recursive(Node * curNode, int32_t * curIdx)
{

	if(curNode->desc1 == NULL)
		return;
	int32_t numDesc1, numDesc2, diff;
	numDesc1 = Node_calc_num_descendents(curNode->desc1);
	numDesc2 = Node_calc_num_descendents(curNode->desc2);
	diff = numDesc1 - numDesc2;
	*curIdx += (diff > 0 ? diff : -1*diff);
	Node_calc_colless_index_recursive(curNode->desc1, curIdx);
	Node_calc_colless_index_recursive(curNode->desc2, curIdx);
	return;
}

// TODO: Make this recursive with a private_ function and update all
// calls.
void Node_add_to_size(Node * curNode, int32_t * curSize)
{
	if(curNode->desc1 == NULL)
	{
		fprintf(stderr, "terminal node called in Node_add_to_size\n");
		return;
	}
	if((curNode->desc1)->numMuts == 0 && (curNode->desc1)->desc1 != NULL)
	{
		//(*curSize)++;
		Node_add_to_size(curNode->desc1, curSize);
	}
	else
		(*curSize)++;
	if((curNode->desc2)->numMuts == 0 && (curNode->desc2)->desc2 != NULL)
	{
		//(*curSize)++;
		Node_add_to_size(curNode->desc2, curSize);
	}
	else
		(*curSize)++;
	return;
}
