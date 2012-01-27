#ifndef _READDIST_H_
#define _READDIST_H_ 1

#include <stdio.h>

#include "pih_lib.h"

typedef struct __value_pair   *value_pair_t;
typedef struct __values_table *values_table_t;

values_table_t read_dist_file_from_phylip(char *file_name);
values_table_t read_dist_file_from_paup(char *file_name);
void write_dist_file(char *file_name, hash_table_t dist);

values_table_t values_table_create(size_t size);
void values_table_destroy(values_table_t *values_table);

size_t values_table_get_size(values_table_t values_table);
size_t values_table_add_name(values_table_t values_table, size_t, char *name);
char *values_table_get_name(values_table_t values_table, size_t pos);
void* values_table_add_value(values_table_t values_table, size_t pos_1, size_t pos_2, double value);
double values_table_get_value(values_table_t values_table, size_t pos_1, size_t pos_2);
void values_table_print(FILE* output, values_table_t values_table);
iterator_t values_table_name_iterator(values_table_t values_table);

int check_text(char* text, char *expected, char *error_menssage);

double** vector_d_create(values_table_t v, size_t *size);
double** vector_d_popule(values_table_t v, double **d, size_t d_size);


char *int_to_string(values_table_t values_table, size_t number);

struct __value_pair
{
	size_t taxon_1;
	size_t taxon_2;
	double value;
};

struct __values_table
{
	size_t 		 size;
	char** 		 number_cache;
	hash_table_t names;
	hash_table_t values;
};

#endif
