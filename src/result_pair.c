/***************************************************************************
* result_pair.c
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
 
# result_pair.c:
# Save the result of a pair and do the conversion
 
****************************************************************/


#include <assert.h>

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "result_pair.h"
#include "pih_lib.h"
#include "scop2dist.h"

// teste
#ifdef _TEST_RP_
int main()
{
	char *seq1 = "seq_1";
	char *seq0 = "seq_0";
	char *seq2 = "seq_2";
	char *seq3 = "seq_3";
	char *seq4 = "seq_4";
	char *seq5 = "seq_5";
	char *seq6 = "seq_6";
	char *seq7 = "seq_7";
	char *seq8 = "seq_8";
	char *seq9 = "seq_9";

	result_pair_t rp_00 = result_pair_create(seq0, seq0, 0.0);
	result_pair_t rp_01 = result_pair_create(seq0, seq1, 0.1);
	result_pair_t rp_02 = result_pair_create(seq0, seq2, 0.2);
	result_pair_t rp_03 = result_pair_create(seq0, seq3, 0.3);
	result_pair_t rp_04 = result_pair_create(seq0, seq4, 0.4);
	result_pair_t rp_05 = result_pair_create(seq0, seq5, 0.5);
	result_pair_t rp_06 = result_pair_create(seq0, seq6, 0.6);
	result_pair_t rp_07 = result_pair_create(seq0, seq7, 0.7);
	result_pair_t rp_08 = result_pair_create(seq0, seq8, 0.8);
	result_pair_t rp_09 = result_pair_create(seq0, seq9, 0.9);

	result_pair_t rp_11 = result_pair_create(seq1, seq1, 1.1);
	result_pair_t rp_12 = result_pair_create(seq1, seq2, 1.2);
	result_pair_t rp_13 = result_pair_create(seq1, seq3, 1.3);
	result_pair_t rp_14 = result_pair_create(seq1, seq4, 1.4);
	result_pair_t rp_15 = result_pair_create(seq1, seq5, 1.5);
	result_pair_t rp_16 = result_pair_create(seq1, seq6, 1.6);
	result_pair_t rp_17 = result_pair_create(seq1, seq7, 1.7);
	result_pair_t rp_18 = result_pair_create(seq1, seq8, 1.8);
	result_pair_t rp_19 = result_pair_create(seq1, seq9, 1.9);

	result_pair_t rp_22 = result_pair_create(seq2, seq2, 2.2);
	result_pair_t rp_23 = result_pair_create(seq2, seq3, 2.3);
	result_pair_t rp_24 = result_pair_create(seq2, seq4, 2.4);
	result_pair_t rp_25 = result_pair_create(seq2, seq5, 2.5);
	result_pair_t rp_26 = result_pair_create(seq2, seq6, 2.6);
	result_pair_t rp_27 = result_pair_create(seq2, seq7, 2.7);
	result_pair_t rp_28 = result_pair_create(seq2, seq8, 2.8);
	result_pair_t rp_29 = result_pair_create(seq2, seq9, 2.9);

	result_pair_t rp_33 = result_pair_create(seq3, seq3, 3.3);
	result_pair_t rp_34 = result_pair_create(seq3, seq4, 3.4);
	result_pair_t rp_35 = result_pair_create(seq3, seq5, 3.5);
	result_pair_t rp_36 = result_pair_create(seq3, seq6, 3.6);
	result_pair_t rp_37 = result_pair_create(seq3, seq7, 3.7);
	result_pair_t rp_38 = result_pair_create(seq3, seq8, 3.8);
	result_pair_t rp_39 = result_pair_create(seq3, seq9, 3.9);

	result_pair_t rp_44 = result_pair_create(seq4, seq4, 4.4);
	result_pair_t rp_45 = result_pair_create(seq4, seq5, 4.5);
	result_pair_t rp_46 = result_pair_create(seq4, seq6, 4.6);
	result_pair_t rp_47 = result_pair_create(seq4, seq7, 4.7);
	result_pair_t rp_48 = result_pair_create(seq4, seq8, 4.8);
	result_pair_t rp_49 = result_pair_create(seq4, seq9, 4.9);

	result_pair_t rp_55 = result_pair_create(seq5, seq5, 5.5);
	result_pair_t rp_56 = result_pair_create(seq5, seq6, 5.6);
	result_pair_t rp_57 = result_pair_create(seq5, seq7, 5.7);
	result_pair_t rp_58 = result_pair_create(seq5, seq8, 5.8);
	result_pair_t rp_59 = result_pair_create(seq5, seq9, 5.9);

	result_pair_t rp_66 = result_pair_create(seq6, seq6, 6.6);
	result_pair_t rp_67 = result_pair_create(seq6, seq7, 6.7);
	result_pair_t rp_68 = result_pair_create(seq6, seq8, 6.8);
	result_pair_t rp_69 = result_pair_create(seq6, seq9, 6.9);

	result_pair_t rp_77 = result_pair_create(seq7, seq7, 7.7);
	result_pair_t rp_78 = result_pair_create(seq7, seq8, 7.8);
	result_pair_t rp_79 = result_pair_create(seq7, seq9, 7.9);

	result_pair_t rp_88 = result_pair_create(seq8, seq8, 8.8);
	result_pair_t rp_89 = result_pair_create(seq8, seq9, 8.9);

	result_pair_t rp_99 = result_pair_create(seq9, seq9, 9.9);


	hash_table_t result_table = hash_table_create();

	/* Populando 0 */
	hash_table_t this_table = hash_table_create();

	hash_table_add(this_table, rp_00->sequence_2, rp_00);
	hash_table_add(this_table, rp_01->sequence_2, rp_01);
	hash_table_add(this_table, rp_02->sequence_2, rp_02);
	hash_table_add(this_table, rp_03->sequence_2, rp_03);
	hash_table_add(this_table, rp_04->sequence_2, rp_04);
	hash_table_add(this_table, rp_05->sequence_2, rp_05);
	hash_table_add(this_table, rp_06->sequence_2, rp_06);
	hash_table_add(this_table, rp_07->sequence_2, rp_07);
	hash_table_add(this_table, rp_08->sequence_2, rp_08);
	hash_table_add(this_table, rp_09->sequence_2, rp_09);

	hash_table_add(result_table, rp_00->sequence_1, this_table);

	/* Populando 1 */
	this_table = hash_table_create();

	hash_table_add(this_table, rp_11->sequence_2, rp_11);
	hash_table_add(this_table, rp_12->sequence_2, rp_12);
	hash_table_add(this_table, rp_13->sequence_2, rp_13);
	hash_table_add(this_table, rp_14->sequence_2, rp_14);
	hash_table_add(this_table, rp_15->sequence_2, rp_15);
	hash_table_add(this_table, rp_16->sequence_2, rp_16);
	hash_table_add(this_table, rp_17->sequence_2, rp_17);
	hash_table_add(this_table, rp_18->sequence_2, rp_18);
	hash_table_add(this_table, rp_19->sequence_2, rp_19);

	hash_table_add(result_table, rp_11->sequence_1, this_table);

	/* Populando 2 */
	this_table = hash_table_create();

	hash_table_add(this_table, rp_22->sequence_2, rp_22);
	hash_table_add(this_table, rp_23->sequence_2, rp_23);
	hash_table_add(this_table, rp_24->sequence_2, rp_24);
	hash_table_add(this_table, rp_25->sequence_2, rp_25);
	hash_table_add(this_table, rp_26->sequence_2, rp_26);
	hash_table_add(this_table, rp_27->sequence_2, rp_27);
	hash_table_add(this_table, rp_28->sequence_2, rp_28);
	hash_table_add(this_table, rp_29->sequence_2, rp_29);

	hash_table_add(result_table, rp_22->sequence_1, this_table);

	/* Populando 3 */
	this_table = hash_table_create();

	hash_table_add(this_table, rp_33->sequence_2, rp_33);
	hash_table_add(this_table, rp_34->sequence_2, rp_34);
	hash_table_add(this_table, rp_35->sequence_2, rp_35);
	hash_table_add(this_table, rp_36->sequence_2, rp_36);
	hash_table_add(this_table, rp_37->sequence_2, rp_37);
	hash_table_add(this_table, rp_38->sequence_2, rp_38);
	hash_table_add(this_table, rp_39->sequence_2, rp_39);

	hash_table_add(result_table, rp_33->sequence_1, this_table);

	/* Populando 4 */
	this_table = hash_table_create();

	hash_table_add(this_table, rp_44->sequence_2, rp_44);
	hash_table_add(this_table, rp_45->sequence_2, rp_45);
	hash_table_add(this_table, rp_46->sequence_2, rp_46);
	hash_table_add(this_table, rp_47->sequence_2, rp_47);
	hash_table_add(this_table, rp_48->sequence_2, rp_48);
	hash_table_add(this_table, rp_49->sequence_2, rp_49);

	hash_table_add(result_table, rp_44->sequence_1, this_table);

	/* Populando 5 */
	this_table = hash_table_create();

	hash_table_add(this_table, rp_55->sequence_2, rp_55);
	hash_table_add(this_table, rp_56->sequence_2, rp_56);
	hash_table_add(this_table, rp_57->sequence_2, rp_57);
	hash_table_add(this_table, rp_58->sequence_2, rp_58);
	hash_table_add(this_table, rp_59->sequence_2, rp_59);

	hash_table_add(result_table, rp_55->sequence_1, this_table);

	/* Populando 6 */
	this_table = hash_table_create();

	hash_table_add(this_table, rp_66->sequence_2, rp_66);
	hash_table_add(this_table, rp_67->sequence_2, rp_67);
	hash_table_add(this_table, rp_68->sequence_2, rp_68);
	hash_table_add(this_table, rp_69->sequence_2, rp_69);

	hash_table_add(result_table, rp_66->sequence_1, this_table);

	/* Populando 7 */
	this_table = hash_table_create();

	hash_table_add(this_table, rp_77->sequence_2, rp_77);
	hash_table_add(this_table, rp_78->sequence_2, rp_78);
	hash_table_add(this_table, rp_79->sequence_2, rp_79);

	hash_table_add(result_table, rp_77->sequence_1, this_table);

	/* Populando 8 */
	this_table = hash_table_create();

	hash_table_add(this_table, rp_88->sequence_2, rp_88);
	hash_table_add(this_table, rp_89->sequence_2, rp_89);

	hash_table_add(result_table, rp_88->sequence_1, this_table);

	/* Populando 9 */
	this_table = hash_table_create();

	hash_table_add(this_table, rp_99->sequence_2, rp_99);

	hash_table_add(result_table, rp_99->sequence_1, this_table);


	convert_to_distmat(result_table);
	fprintf(stderr, "FIM!\n");

	return 0;
}

#endif

result_pair_t result_pair_create(char *seq1, char *seq2, double value)
{
	result_pair_t rp;

	rp = (result_pair_t) malloc(sizeof(struct __result_pair));
	assert(rp != NULL);
	memset(rp, '\0', sizeof(struct __result_pair));

	rp->sequence_1 = seq1;
	rp->sequence_2 = seq2;
	rp->value = value;

	return rp;
}

int result_pair_destroy(void *data)
{
	result_pair_t *rp = (result_pair_t *) data;

	if (*rp == NULL) {
		return 0;
	}

//	fprintf(stderr, "destruindo RP (%p) %s %s \n", *rp, (*rp)->sequence_1, (*rp)->sequence_2);
	free(*rp);
	*rp = NULL;

	return 1;
}

void* result_pair_clone(void *data)
{
	result_pair_t rp = (result_pair_t) data;

	result_pair_t clone = (result_pair_t) malloc(sizeof(struct __result_pair));
	assert(clone != NULL);

	clone->sequence_1 = rp->sequence_1;
	clone->sequence_2 = rp->sequence_2;
	clone->value = rp->value;

	return (void *) clone;
}


void result_pair_print(char *id, void *data)
{

	result_pair_t rp = (result_pair_t) data;
	assert(rp != NULL);

	fprintf(stderr, "%s%s X %s : %10.2e\n", id, rp->sequence_1, rp->sequence_2, rp->value);
}

int get_pos(char *id, hash_table_t table_position, int *actual)
{
	int pos = (int ) hash_table_get(table_position, id);

	if (pos == 0) {
		(*actual)++;
		pos = *actual; // Somar um para nao confundir o 0 com o NULL do hash_table_get
		hash_table_add(table_position, id, (void *) pos);
	}

	return pos -1;
}

/*
 * A entrada deste eh um hash contando varios hashes sendo cada um de uma sequencia vs outras
 */
DISTMAT *convert_to_distmat(hash_table_t result_table)
{
	DISTMAT* distmat = NULL;
	char *seq1 = NULL;
	char *seq2 = NULL;
	int actual_pos = 0;
	cell_t cell = NULL;
	cell_t cell_inner = NULL;
	hash_table_t table_position = hash_table_create();
	hash_table_t table_seq = NULL;

	assert(result_table != NULL);

	distmat = DISTMATalloc(result_table->count);

	iterator_t table_result_iterator = hash_table_iterator(result_table);
	iterator_t table_seq_iterator = NULL;

	result_pair_t rp = NULL;

	int pos, pos_inner;

	while (table_result_iterator->has_next(table_result_iterator)) {
		cell = table_result_iterator->next(table_result_iterator);

		seq1 = cell->id;
		pos = get_pos(seq1, table_position, &actual_pos);

		table_seq = (hash_table_t) cell->data;
		table_seq_iterator = hash_table_iterator(table_seq);

		while (table_seq_iterator->has_next(table_seq_iterator)) {
			cell_inner = table_seq_iterator->next(table_seq_iterator);
			seq2 = cell_inner->id;
			rp = (result_pair_t) cell_inner->data;

			pos_inner = get_pos(seq2, table_position, &actual_pos);

			distmat->dist[pos][pos_inner] = distmat->dist[pos_inner][pos] = rp->value;
		}

	}

	iterator_t table_position_iterator = hash_table_iterator(table_position);

	while (table_position_iterator->has_next(table_position_iterator)) {
		cell = table_position_iterator->next(table_position_iterator);

		distmat->taxa[((int) cell->data) - 1] = cell->id;
	}


	return distmat;
}


