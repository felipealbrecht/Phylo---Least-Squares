#include <assert.h>

#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include "result_pair.h"
#include "pih_lib.h"

#include "readdist.h"

/*
from: http://evolution.genetics.washington.edu/phylip/doc/distance.html
 
 The input format for distance data is straightforward.
 The first line of the input file contains the number of species. 
 There follows species data, starting, as with all other programs, with a species name. 
 The species name is ten characters long, and must be padded out with blanks if shorter. 
 For each species there then follows a set of distances to all the other species 
  (options selected in the programs' menus allow the distance matrix to be upper or lower triangular or square). 
 The distances can continue to a new line after any of them. 
 
 If the matrix is lower-triangular, 
  the diagonal entries (the distances from a species to itself) will not be read by the programs. 
  
  If they are included anyway,   they will be ignored by the programs, 
  except for the case where one of them starts a new line, in which case the program will mistake it for a species name and get very confused.
 
 
Example:
  -----//cut here//----------------
     5
	 Alpha      0.000 1.000 2.000 3.000 3.000
	 Beta       1.000 0.000 2.000 3.000 3.000
	 Gamma      2.000 2.000 0.000 3.000 3.000
	 Delta      3.000 3.000 0.000 0.000 1.000
	 Epsilon    3.000 3.000 3.000 1.000 0.000
  -----//cut here//----------------
 
Example triagonal:
  -----//cut here//----------------
   14
   Mouse     
   Bovine      1.7043
   Lemur       2.0235  1.1901
   Tarsier     2.1378  1.3287  1.2905
   Squir Monk  1.5232  1.2423  1.3199  1.7878
   Jpn Macaq   1.8261  1.2508  1.3887  1.3137  1.0642
   Rhesus Mac  1.9182  1.2536  1.4658  1.3788  1.1124  0.1022
   Crab-E.Mac  2.0039  1.3066  1.4826  1.3826  0.9832  0.2061  0.2681
   BarbMacaq   1.9431  1.2827  1.4502  1.4543  1.0629  0.3895  0.3930  0.3665
   Gibbon      1.9663  1.3296  1.8708  1.6683  0.9228  0.8035  0.7109  0.8132  0.7858
   Orang       2.0593  1.2005  1.5356  1.6606  1.0681  0.7239  0.7290  0.7894  0.7140  0.7095
   Gorilla     1.6664  1.3460  1.4577  1.5935  0.9127  0.7278  0.7412  0.8763  0.7966  0.5959  0.4604
   Chimp       1.7320  1.3757  1.7803  1.7119  1.0635  0.7899  0.8742  0.8868  0.8288  0.6213  0.5065  0.3502
   Human       1.7101  1.3956  1.6661  1.7599  1.0557  0.6933  0.7118  0.7589  0.8542  0.5612  0.4700  0.3097  0.2712
  -----//cut here//----------------
 
  \todo: cuidado com os memory leaks, principalmente depois de usar o get_text, deve-se liberar a memoria!
*/

#ifdef _DEBUG_READDIST_
int main(int argc, char **argv)
{

	values_table_t vt = read_dist_file_from_paup("../data/domains.NX");
	FILE *output = stdout;

	if (argc == 2) {
		output = fopen(argv[1], "w+");
	}

	values_table_print(output, vt);

        values_table_destroy(&vt);

	return 0;
}
#endif


static inline int pass_spaces(FILE *fp)
{
	int count = 0;
	char c;

	while ( isspace(c = getc(fp)) ) {
		count++;
	}

	/* back the last one */
	count--;
	ungetc(c, fp);

	return count;
}

#define STRING_SIZE 256

static inline char* get_text(FILE *fp)
{
	char *tmp = (char *) malloc(sizeof(char) * STRING_SIZE);
	memset(tmp, '\0', STRING_SIZE);

	fscanf(fp, "%s", tmp);

	return tmp;
}

static inline int is_number(char *number)
{
	assert(number != NULL);

	int len = strlen(number);
	int i;
	int dots = 0; /* To make the use of dots, Valid: 2.0 . Invalid: 2.6.2 */

	for (i = 0; i < len; i++) {

		if (number[i] == '.') {
			if (dots == 1) {
				return 0;
			} else {
				dots = 1;
			}

		} else {

			if (!isdigit(number[i])) {
				return 0;
			}
		}
	}

	return 1;
}


value_pair_t value_pair_create(int taxon_1, int taxon_2, double value)
{
	value_pair_t vp = (value_pair_t) malloc(sizeof(struct __value_pair));

	assert(vp != NULL);

	vp->taxon_1 = taxon_1;
	vp->taxon_2 = taxon_2;
	vp->value   = value;

	return vp;
}

int value_pair_destroy(void *v)
{
	assert(v != NULL);
	value_pair_t *value_pair = (value_pair_t *) v;
	free(*value_pair);
	value_pair = NULL;

	return 1;
}

char *int_to_string(values_table_t values_table, size_t number)
{
	assert(number < values_table->size);

	int chars = 2;
        size_t num = number;
	char *n = NULL;


	if (values_table->number_cache[number] != NULL) {
		return values_table->number_cache[number];
	}

	while (num / 10 > 0) {
		chars++;
		num = num / 10;
	}

	n = calloc(sizeof(char), chars);
	sprintf(n, "%ld", number);

	values_table->number_cache[number] = n;

	return n;
}

values_table_t values_table_create(size_t size)
{
	values_table_t values_table = (values_table_t) malloc(sizeof(struct __values_table));

	assert(values_table != NULL);

	values_table->size = size;

	values_table->number_cache = (char **) malloc(sizeof(char *) * size);
	memset(values_table->number_cache, '\0', sizeof(char *) * size);
	assert(values_table->number_cache != NULL);

	values_table->names = hash_table_create();
	assert(values_table->names != NULL);

	values_table->values = hash_table_create();
	assert(values_table->values != NULL);

	return values_table;
}

int destroy_values(void *v)
{
	assert(v != NULL);
	hash_table_t* value_table = (hash_table_t *) v;
	hash_table_destroy_all(value_table, value_pair_destroy);
	*value_table = NULL;

	return 1;
}

void values_table_destroy(values_table_t *values_table)
{
	size_t i;

	if (*values_table == NULL) {
		return;
	}

	hash_table_destroy_all(&(*values_table)->names, destroy_string);
	hash_table_destroy_all(&(*values_table)->values, destroy_values);

	for (i = 0; i < (*values_table)->size; i++) {
		free((*values_table)->number_cache[i]);
	}
	free((*values_table)->number_cache);

	free(*values_table);
	values_table = NULL;

}

unsigned int values_table_get_size(values_table_t values_table)
{
	assert(values_table != NULL);

	return values_table->names->keys->size;
}

/**
 * Return  0: Okay
 * Return -1: The pos alread in use
 */
unsigned int values_table_add_name(values_table_t values_table, unsigned int pos, char *name)
{
	assert(values_table != NULL);
	assert(name != NULL);

	char *pos_to_string = int_to_string(values_table, pos);
	void* value = hash_table_get(values_table->names, pos_to_string);

	if (value != NULL) {
		return -1;
	} else {
		hash_table_add(values_table->names, pos_to_string, (void *) name);
		return 0;
	}
}

iterator_t values_table_name_iterator(values_table_t values_table)
{
	iterator_t iterator;

	iterator = (iterator_t) malloc(sizeof(struct __iterator));
	assert(iterator != NULL);

	iterator->list = values_table->names->keys;
	iterator->actual = NULL;
	iterator->has_next = __list_has_next;
	iterator->next = __list_next;

	return iterator;
}

char *values_table_get_name(values_table_t values_table, unsigned int pos)
{
	assert(values_table != NULL);

	char *pos_to_string = int_to_string(values_table, pos);
	char *name = hash_table_get(values_table->names, pos_to_string);

	return name;
}

void* values_table_add_value(values_table_t values_table, unsigned int pos_1, unsigned int pos_2, double value)
{
	assert(values_table != NULL);

	char *pos_bigger_to_string = NULL;
	char *pos_smaller_to_string = NULL;
	unsigned int max;
	unsigned int min;

	hash_table_t values = NULL;

	if (pos_1 >= pos_2) {
		pos_bigger_to_string = int_to_string(values_table, pos_1);
		pos_smaller_to_string = int_to_string(values_table, pos_2);
		max = pos_1;
		min = pos_2;
	} else {
		pos_bigger_to_string = int_to_string(values_table, pos_2);
		pos_smaller_to_string = int_to_string(values_table, pos_1);
		max = pos_2;
		min = pos_1;
	}

	values = hash_table_get(values_table->values, pos_bigger_to_string);
	if (values == NULL) {
		values = hash_table_create();
		hash_table_add(values_table->values, pos_bigger_to_string, values);
	}

	return hash_table_add(values, pos_smaller_to_string, value_pair_create(max, min, value));
}


/**TODO: Fazer uma versao que ja receba strings como parematro! */
double values_table_get_value(values_table_t values_table, unsigned int pos_1, unsigned int pos_2)
{
	assert(values_table != NULL);

	char *pos_bigger_to_string = NULL;
	char *pos_smaller_to_string = NULL;

	value_pair_t vp;

	hash_table_t values = NULL;

	if (pos_1 >= pos_2) {
		pos_bigger_to_string = int_to_string(values_table, pos_1);
		pos_smaller_to_string = int_to_string(values_table, pos_2);
	} else {
		pos_bigger_to_string = int_to_string(values_table, pos_2);
		pos_smaller_to_string = int_to_string(values_table, pos_1);
	}

	values = hash_table_get(values_table->values, pos_bigger_to_string);
	assert (values != NULL);

	vp = hash_table_get(values, pos_smaller_to_string);

	return vp->value;
}

values_table_t read_dist_file_from_phylip(char *file_name)
{
	assert(file_name != NULL);

	FILE* fp = NULL;
	char *text;
	char *name = NULL;
	char *aux = NULL;
	unsigned int j, i, qtd_taxons;

	fp = fopen(file_name, "r");
	if (fp == NULL) {
		return NULL;
	}

	pass_spaces(fp);

	text = get_text(fp);

	if (!is_number(text)) {
		/* SHIT! File error! */
		/* \todo: do a erro message */
		assert(1 == 2);
	}

	qtd_taxons = atoi(text);
	free(text);

	values_table_t values_table = values_table_create(qtd_taxons);

	/* ler o nome dos taxons */
	for (i = 0; i < qtd_taxons; i++) {
		name = get_text(fp);

		values_table_add_name(values_table, i, name);

		for (j = 0; j < qtd_taxons; j++) {
			pass_spaces(fp);
			aux = get_text(fp);
			if (!is_number(aux)) {
				fprintf(stderr, "\nO text '%s' deve ser um valor numerico\n", aux);
				fclose(fp);
				return NULL;
			}
			if (i <= j) {
				values_table_add_value(values_table, i, j, atof(aux));
			}
			free(aux);
		}
	}

	fclose(fp);
	return values_table;
}

int check_text(char* text, char *expected, char *error_message)
{
	assert(text != NULL);
	assert(expected != NULL);

	if (!(strcmp(text, expected) == 0)) {
		fprintf(stderr, error_message);
		return 0;
	}

	return 1;
}


char* get_from_expression(char *expr)
{
	assert(expr != NULL);

	char *number;
	int i, len, initial = -1, final = -1;

	len = strlen(expr);

	for (i = 0; i < len; i++) {
		if (!isdigit(expr[i])) {
			if (initial == -1) {
				continue;
			} else {
				final = i;
				break;
			}

		} else {
			if (initial == -1) {
				initial = i;
			}
		}
	}

	int length = final - initial + 1;

	number = malloc(sizeof(char) * length);
	memset(number, '\0', sizeof(char) * length);

	for (i = initial; i < final; i++) {
		number[i - initial] = expr[i];
	}

	return number;
}

/*
 * A really STUPID parser for PAUP* distance files
 */
values_table_t read_dist_file_from_paup(char *file_name)
{
	assert(file_name != NULL);

	FILE* fp = fopen(file_name, "r");
	char *text;

	assert(fp != NULL);

	text = get_text(fp);
	if (!check_text(text, "#NEXUS", "'#NEXUS' nao foi encontrado\n")) {
		return NULL;
	}
	free(text);


	text = get_text(fp);
	if (!check_text(text, "begin", "'begin' nao foi encontrado\n")) {
		return NULL;
	}
	free(text);

	text = get_text(fp);
	if (!check_text(text, "taxa;", "'taxa;' nao foi encontrado\n")) {
		return NULL;
	}
	free(text);

	text = get_text(fp);
	if (!check_text(text, "dimensions", "'dimensions' nao foi encontrado\n")) {
		return NULL;
	}
	free(text);

	text = get_text(fp);
	char *number = get_from_expression(text);
	if (number == NULL) {
		return NULL;
	}
	free(text);

	int qtd = atoi(number);
	free(number);

	text = get_text(fp);
	if (!check_text(text, "taxlabels", "'taxlabels' nao foi encontrado\n")) {
		return NULL;
	}
	free(text);

	values_table_t values_table = values_table_create(qtd);

	int i = 0;
	for (i = 0; i < qtd; i++) {
		text = get_text(fp);
		values_table_add_name(values_table, i, text);
	}

	text = get_text(fp);
	if (!check_text(text, ";", "';' nao foi encontrado\n")) {
		return NULL;
	}
	free(text);


	text = get_text(fp);
	if (!check_text(text, "end;", "'end;' nao foi encontrado\n")) {
		return NULL;
	}
	free(text);

	text = get_text(fp);
	if (!check_text(text, "begin", "'begin' nao foi encontrado\n")) {
		return NULL;
	}
	free(text);

	text = get_text(fp);
	if (!check_text(text, "distances;", "'distances;' nao foi encontrado\n")) {
		return NULL;
	}
	free(text);

	text = get_text(fp);
	if (!check_text(text, "format", "'format' nao foi encontrado\n")) {
		return NULL;
	}
	free(text);

	text = get_text(fp);
	if (!check_text(text, "triangle=lower;", "'triangle=lower;' nao foi encontrado\n")) {
		return NULL;
	}
	free(text);

	text = get_text(fp);
	if (!check_text(text, "matrix", "'matrix' nao foi encontrado\n")) {
		return NULL;
	}
	free(text);

	text = get_text(fp);
	if (!check_text(text, "[", "'[' nao foi encontrado\n")) {
		return NULL;
	}
	free(text);

	char *value = NULL;
	// Uses the content for to do a check...
	for (i = 0; i < qtd; i++) {
		text = get_text(fp);
		value = values_table_get_name(values_table, i);
		assert(text != NULL);
		assert(value != NULL);
		assert(strcmp(text, value) == 0);
		free(text);
	}

	text = get_text(fp);
	if (!check_text(text, "]", "']' nao foi encontrado\n")) {
		return NULL;
	}
	free(text);

	int j;
	for (i = 0; i < qtd; i++) {
		text = get_text(fp);
		free(text);
		for(j = 0; j <= i; j++) {
			text = get_text(fp);
			if (!is_number(text)) {
				fprintf(stderr, "\nO text '%s' deve ser um valor numerico\n", text);
				return NULL;
			}
			values_table_add_value(values_table, i, j, atof(text));
			free(text);
		}
	}

	fclose(fp);
	return values_table;
}


void values_table_print(FILE *output, values_table_t values_table)
{
	unsigned int i, j;
	unsigned int qtd_taxons = values_table_get_size(values_table);
	char *name = NULL;

	fprintf(output, "Taxons: %d\n", qtd_taxons);

	for (i = 0; i < qtd_taxons; i++) {
		name = values_table_get_name(values_table, i);
		fprintf(output, "%-10s", (char *) name);

		for (j = 0; j < qtd_taxons; j++) {
			fprintf(output, "%.3f ", values_table_get_value(values_table, i, j));
		}
		fprintf(output, "\n");
	}
}


double** vector_d_create(values_table_t v, size_t *size)
{
	size_t i;
	size_t s = values_table_get_size(v);

	assert(s > 1);
	size_t total = --s;
	while (--s) {
		total += s;
	}

	double** dd  = (double **) malloc(sizeof(double *) * total);
	assert(dd != NULL);

	for (i = 0; i < total; i++) {
		dd[i] = (double *) malloc(sizeof(double) * 1);
		assert(dd[i] != NULL);
		dd[i][0] = 0.0;
	}

	*size = total;

	return dd;
}


double ** vector_d_popule(values_table_t v, double **d, size_t d_size)
{
	list_t keys = list_sort_by_data(v->names->keys, compare_strings);

	iterator_t it = list_iterator(keys);
	iterator_t it_2 = NULL;
	cell_t cell, cell_2;
	size_t total = 0;

	while (it->has_next(it)) {
		cell = it->next(it);
		it_2 = list_iterator_continues(it);
		
		while (it_2->has_next(it_2)) {
			cell_2 = it_2->next(it_2);
			//fprintf(stderr, "Distancia entre '%s'(%s)  e '%s'(%s) eh %f\n", 
			//	(char *) cell->data, cell->id, (char *) cell_2->data, cell_2->id,
			//	values_table_get_value(v, atoi(cell->id), atoi(cell_2->id))
			//	);
			d[total][0] = values_table_get_value(v, atoi(cell->id), atoi(cell_2->id));
			total++;
			assert(total <= d_size);
		}
		list_iterator_destroy(&it_2);
	}
	return d;
}
