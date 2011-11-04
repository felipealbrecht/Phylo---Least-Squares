#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pih_lib.h"

#include "tree.h"

#define CHECK(v) if(v==NULL){fprintf(stderr, "Erro de alocacao em %s:%d", __FILE__, __LINE__); return NULL;}

static tree_item_t __item_get_down(char *taxon, tree_item_t item, list_t path);
static tree_item_t __internal_node_get_down(char *taxon, internal_node_t internal_node, list_t path);
static tree_item_t __tree_get_down(char *taxon, tree_t item, list_t path);

#ifdef _DEBUG_
static char* create_id(char* id_one, char *id_two, char *sep);
#endif

#define GET_PARENT(i) i->parent
#define TREE_ROOT "tree-root"

list_t tree_leaf_walk(taxon_node_t leaf_from, tree_t tree);
list_t leaf_walk_to(taxon_node_t leaf_from, taxon_node_t leaf_to);


double** tree_create_matrix(tree_t tree, size_t *lines, size_t *columns)
{
	assert(tree != NULL);

	list_t path;
	list_t taxons = list_sort(tree->taxon_nodes->keys, compare_strings);
	iterator_t it, it_2, in_iterator;
	cell_t cell, cell_2, in_cell;
	taxon_node_t leaf_from, leaf_to;

	size_t i,j;
	size_t c = combination(taxons->size, 2);
	double **matrix = (double **) malloc(sizeof(double) * c);

	CHECK(matrix);
	for (i = 0; i < c; i++) {
		matrix[i] = (double *) malloc(sizeof(double) * tree->internal_nodes->keys->size);
		CHECK(matrix[i]);
	}

	it = list_iterator(taxons);
	i = 0;
	while (it->has_next(it)) {
		cell = (cell_t) it->next(it);
		leaf_from = (taxon_node_t) cell->data;
		it_2 = list_iterator_continues(it);
		while (it_2->has_next(it_2)) {
			cell_2 = (cell_t) it_2->next(it_2);
			leaf_to = (taxon_node_t) cell_2->data;

			path = leaf_walk_to(leaf_from, leaf_to);
                        
			j = 0;
			in_iterator = list_iterator(tree->internal_nodes->keys);
			while (in_iterator->has_next(in_iterator)) {
				in_cell = (cell_t) in_iterator->next(in_iterator);
				if (list_search(path, in_cell->id) == NULL) {
					matrix[i][j] = 0.0;
				} else {
					matrix[i][j] = 1.0;
				}
				j++;
			}

			list_destroy(&path);
			list_iterator_destroy(&in_iterator);
			i++;
		}
                list_iterator_destroy(&it_2);
	}
        list_iterator_destroy(&it);

        free(taxons);

	*lines = c;
	*columns = tree->internal_nodes->keys->size;

	return matrix;
}

void tree_matrix_destroy(double ***matrix, size_t lines)
{
	size_t i;

        for (i = 0; i < lines; i++) {
		free((*matrix)[i]);
	}
	free(*matrix);
	*matrix = NULL;
}


list_t leaf_walk_to(taxon_node_t leaf_from, taxon_node_t leaf_to)
{
	assert(leaf_from != NULL);
	assert(leaf_to != NULL);
	assert(leaf_from->self != NULL);
	assert(leaf_to->self != NULL);
	assert(leaf_from->self->parent != NULL);
	assert(leaf_to->self->parent != NULL);

	list_t path = list_create();

	list_add(path, leaf_from->self->parent_branch->id,
	         leaf_from->self->parent_branch);

	if (leaf_from->self->parent == leaf_to->self->parent) {
		list_add(path, leaf_to->self->parent_branch->id,
		         leaf_from->self->parent_branch);
		return path;
	}

	// UP to the root
	tree_item_t item = leaf_from->self->parent;
	internal_node_t internal_node;
	while (item->type == INTERNAL_NODE) {
		internal_node = item->item.internal_node;
		if (list_search(internal_node->content, leaf_to->taxon) != NULL) {
			break;
		}
		list_add(path, item->parent_branch->id, item->parent_branch);
		item = item->parent;
	}

	__item_get_down(leaf_to->taxon, item, path);

	return path;
}

static tree_item_t __item_get_down(char *taxon, tree_item_t item, list_t path)
{
	assert(taxon != NULL);
	assert(item != NULL);

	if (item->type == INTERNAL_NODE) {
		if (list_search(item->item.internal_node->content, taxon) != NULL) {
			return __internal_node_get_down(taxon, item->item.internal_node, path);
		} else {
			return NULL;
		}


	} else if (item->type == LEAF)  {
		if (item->item.taxon_node->taxon == taxon) {
			return item;
		} else {
			return NULL;
		}

	} else if (item->type == TREE) {
		return __tree_get_down(taxon, item->item.tree, path);
	}


	fprintf(stderr, "Nao encontrado em nenhum sub item! Um Bug!!!!\n");
	assert(1 == 2);
}

static tree_item_t __internal_node_get_down(char *taxon, internal_node_t internal_node, list_t path)
{
	assert(taxon != NULL);
	assert(internal_node != NULL);

	if ((internal_node->item_a->type == INTERNAL_NODE) &&
	        (list_search(internal_node->item_a->item.internal_node->content, taxon) != NULL)) {
		list_add(path, internal_node->branch_a->id, internal_node->branch_a);
		return __internal_node_get_down(taxon, internal_node->item_a->item.internal_node, path);

	} else if ((internal_node->item_b->type == INTERNAL_NODE) &&
	           (list_search(internal_node->item_b->item.internal_node->content, taxon) != NULL)) {
		list_add(path, internal_node->branch_b->id, internal_node->branch_b);
		return __internal_node_get_down(taxon, internal_node->item_b->item.internal_node, path);

	} else if ((internal_node->item_a->type == LEAF) &&
	           (internal_node->item_a->item.taxon_node->taxon == taxon)) {
		list_add(path, internal_node->branch_a->id, internal_node->branch_a);
		return internal_node->item_a;

	} else if ((internal_node->item_b->type == LEAF) &&
	           (internal_node->item_b->item.taxon_node->taxon == taxon)) {
		list_add(path, internal_node->branch_b->id, internal_node->branch_b);
		return internal_node->item_b;

	} else {
		fprintf(stderr, "Nao encontrado em nenhum sub ramo. Um bug!!\n");
		assert(1 == 2);
	}
}



static tree_item_t __tree_get_down(char *taxon, tree_t tree, list_t path)
{
	assert(taxon != NULL);
	assert(tree != NULL);
	tree_item_t child;

	child = __item_get_down(taxon, tree->item_a, path);
	if (child != NULL) {
		list_add(path, tree->branch_a->id, tree->branch_a);
		return child;
	}

	child = __item_get_down(taxon, tree->item_b, path);
	if (child != NULL) {
		list_add(path, tree->branch_b->id, tree->branch_b);
		return child;
	}

	child = __item_get_down(taxon, tree->item_c, path);
	if (child != NULL) {
		list_add(path, tree->branch_c->id, tree->branch_c);
		return child;
	}

	if (child == NULL) {
		fprintf(stderr, "Nao encontrado em nenhum ramo da arvore. Um bug *muito* provavelmente!!\n");
		assert(1 ==2);
	}

	return child;
}

#ifdef _DEBUG_
static char* create_id(char* id_one, char *id_two, char *sep)
{
	assert(id_one != NULL);
	assert(id_two != NULL);

	size_t len_id_one = strlen(id_one);
	size_t len_id_two = strlen(id_two);
	size_t len_sep    = strlen(sep);

	size_t s = len_id_one + len_id_two + len_sep + 1;

	char *path_id = (char *) malloc(sizeof(char) * s);
	memset(path_id, '\0', s);

	memcpy(path_id, id_one, len_id_one);
	memcpy(path_id + len_id_one, sep, len_sep);
	memcpy(path_id + len_id_one + len_sep, id_two, len_id_two);

	return path_id;
}
#endif
