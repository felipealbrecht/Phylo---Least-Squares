/***************************************************************************
* result_pair.h
* (C) 2006 Felipe Albrecht (felipe.albrecht@gmail.com)
 
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 
# result_pair.h:
# Header of save the result of a pair and do the conversion
 
****************************************************************/


#ifndef _RESULT_PAIR_H_
#define _RESULT_PAIR_H_

#include "pih_lib.h"
#include "scop2dist.h"

typedef struct __result_pair *result_pair_t;
struct __result_pair
{
	char *sequence_1;
	char *sequence_2;
	double value;
};

result_pair_t result_pair_create(char *seq1, char *seq2, double value);
int result_pair_destroy(void *data);
void result_pair_print(char *id, void *data);
void* result_pair_clone(void *data);
DISTMAT *convert_to_distmat(hash_table_t result_table);

#endif
