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
