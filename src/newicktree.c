#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "newicktree.h"
#include "definitions.h"

void NewickTree_init(NewickTree * nt)
{
	int32_t i;
	nt->line = (char *)malloc(sizeof(char) * (size_t)100);
	CHECKPOINTER(nt->line);
	nt->line[0] = '\0';
	nt->curIdx = 0;
	nt->curSize = 100;
	return;
}

void NewickTree_add(NewickTree * nt, char * toadd, int32_t toaddLength)
{
	size_t newSize;
	// -1 because of the null character!
	if(nt->curIdx + toaddLength > nt->curSize - 1)
	{
		newSize = nt->curSize + 100;
		nt->line = (char *)realloc(nt->line, sizeof(char) * newSize);
		CHECKPOINTER(nt->line);
		nt->curSize += 100;
	}
	strncat(nt->line, toadd, toaddLength);
	nt->curIdx = index(nt->line,'\0') - nt->line;
	return;
}

void NewickTree_reset(NewickTree * nt)
{
	nt->curIdx = 0;
	nt->line[0] = '\0';
	nt->curSize = 100;
	return;
}

void NewickTree_free(NewickTree * nt)
{
	free(nt->line);
	nt->curIdx = 0;
	nt->curSize = 0;
	return;
}
