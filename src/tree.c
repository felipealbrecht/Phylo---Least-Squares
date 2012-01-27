#include <assert.h>

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "pih_lib.h"

#include "result_pair.h"
#include "tree_matrix.h"
#include "readdist.h"

#include "tree.h"

#define MAX_NAME 256

#define A_TO_INTERNAL 0
#define B_TO_INTERNAL 1
#define C_TO_INTERNAL 2

#define TREE_ROOT "tree-root"

#define CHECK(v) if(v==NULL){fprintf(stderr, "Erro de alocacao em %s:%d", __FILE__, __LINE__); return NULL;}

static int global_iteration;
static int global_iteration_tree_count;

static char* create_id(char* id_one, char *id_two, char *sep);
static char* get_id(tree_item_t item);

static char *tree_item_get_id(tree_item_t item);

static void tree_item_destroy(tree_item_t *tree_item);

static tree_branch_t tree_branch_create(tree_item_t item_a, tree_item_t item_b, tree_t tree);
static void* tree_branch_clone(void *data);
static size_t tree_branch_destroy_fast(void *data);
static size_t tree_branch_destroy(tree_branch_t *tree_branch, tree_t tree);

static void print_repr_tree_item(values_table_t values_table, tree_item_t tree_item, string_buffer_t sb);
static void print_repr_internal_node(values_table_t values_table, internal_node_t internal_node, string_buffer_t sb);


static void parent_add_child(tree_item_t parent, tree_item_t child, tree_t tree);
static void internal_node_add_child(internal_node_t internal_node, tree_item_t child, tree_t tree);
static void tree_add_child(tree_t tree, tree_item_t child);


static tree_item_t parent_remove_child(tree_item_t parent, tree_item_t child, tree_t tree);
static tree_item_t internal_node_remove_child(internal_node_t internal_node, tree_item_t child, tree_t tree);
static tree_item_t tree_remove_child(tree_t tree, tree_item_t child);

static int __tree_update_internal_nodes_distances(tree_item_t tree_item, internal_node_t internal_node_added, tree_item_t actual_node, tree_item_t previous_node, double parent_distance);

static void internal_node_add_content(internal_node_t internal_node, taxon_node_t taxon_node);
static hash_table_t tree_item_get_distances(tree_item_t item);
static tree_item_t tree_item_create(tree_t tree, tree_item_t parent, node_type_t type, void *item, double parent_distance);

static size_t triple_comparer(void *one, void *two);
static size_t triple_destroy(void *data);

/**
 * Calcula a distancia entre os tres taxons e retorna uma matriz com 3 elementos contendo as distancias
 * dos taxons para o no central.
 * Primeiro distancia eh do 'a', segunda do 'b' e tercira do 'c'.
 *
 **/
double* calcule_distances(double d_ab, double d_ac, double d_bc)
{
	double *values = (double *) malloc(sizeof(double) * 3);
	assert(values != NULL);

	values[A_TO_INTERNAL] = (d_ab + d_ac - d_bc) / 2;
	values[B_TO_INTERNAL] = (d_ab + d_bc - d_ac) / 2;
	values[C_TO_INTERNAL] = (d_ac + d_bc - d_ab) / 2;

	return values;
}


static char *tree_item_get_id(tree_item_t item)
{
	if (item == NULL) {
		return "null";
	}

	switch (item->type) {
		case LEAF: return item->item.taxon_node->taxon;
		case INTERNAL_NODE: return item->item.internal_node->id;
		case TREE: return item->item.tree->id;
		default: return NULL;
	}
}

static hash_table_t tree_item_get_distances(tree_item_t item)
{
	assert(item != NULL);

	switch (item->type) {
		case INTERNAL_NODE: return item->item.internal_node->distances;
		case TREE: return item->item.tree->distances;
		default: fprintf(stderr, "O item deve ser um internal node ou tree e eh um %d!!\n", item->type); assert(1 == 2); return NULL;
	}
}

/** 
 * Cria um item da arvore com os dados informados e 
 * faz a associacao birirecional entre o tree_item e o item
 **/
static tree_item_t tree_item_create(tree_t tree, tree_item_t parent, node_type_t type, void *item, double parent_distance)
{
	assert(item != NULL);

	char *id = NULL;
	tree_item_t tree_item = (tree_item_t) malloc(sizeof( struct __tree_item));
	assert(tree_item != NULL);
	memset(tree_item, '\0', sizeof(struct __tree_item));

	tree_item->parent = parent;
	tree_item->type = type;
	tree_item->parent_distance = parent_distance;
	

	switch (type) {
		case LEAF: 	{ 
			tree_item->item.taxon_node = item; 
			((taxon_node_t) item )->self = tree_item; 
			id = ((taxon_node_t) item )->taxon;
			break;
		}
		case INTERNAL_NODE:	{
			tree_item->item.internal_node = item; 
			((internal_node_t) item )->self = tree_item;
			id = ((internal_node_t) item )->id;
			break;
		}
		case TREE: {
			tree_item->item.tree = item; 
			((tree_t) item )->self = tree_item;
			id = ((tree_t) item )->id;
			break;
		}

		default: {
			fprintf(stderr, "O item eh %d\n", type);
			assert(1 == 2);
		}
	}
	
	if (parent != NULL) {
		parent_add_child(parent, tree_item, tree);
	}

	return tree_item;
}

static void parent_add_child(tree_item_t parent, tree_item_t child, tree_t tree)
{
	switch(parent->type) {
		case INTERNAL_NODE: {
			internal_node_add_child(parent->item.internal_node, child, tree);
			if (child->type == LEAF) {
				internal_node_add_content(parent->item.internal_node, child->item.taxon_node);
			}
			break;
		}
		
		case TREE: {
			tree_add_child(parent->item.tree, child);
			break;
		}

		default: {
			fprintf(stderr, "O item do tipo %d nao pode ser parent\n", parent->type);
			assert(1 == 2);
		}
	}
}

static tree_item_t parent_remove_child(tree_item_t parent, tree_item_t child, tree_t tree)
{
	switch(parent->type) {
		case INTERNAL_NODE: {
			return internal_node_remove_child(parent->item.internal_node, child, tree);
		}
		
		case TREE: {
			return tree_remove_child(parent->item.tree, child);
		}

		default: {
			fprintf(stderr, "O item do tipo %d nao pode ser parent\n", parent->type);
			assert(1 == 2);
		}
	}
}


static void internal_node_add_child(internal_node_t internal_node, tree_item_t child, tree_t tree)
{
	assert(internal_node != NULL);
	assert(child != NULL);

	assert(child != internal_node->item_a);
	assert(child != internal_node->item_b);

	if (internal_node->item_a == NULL) {
		internal_node->item_a = child;
		internal_node->branch_a = tree_branch_create(internal_node->self, child, tree);
		child->parent_branch = internal_node->branch_a;
	} else if (internal_node->item_b == NULL) {
		internal_node->item_b = child;
		internal_node->branch_b = tree_branch_create(internal_node->self, child, tree);
		child->parent_branch = internal_node->branch_b;
	} else {
		fprintf(stderr, "O internal_node %s nao tem mais espacos para filhos!\n", internal_node->id);
		assert(1 == 2);
	}
}


static tree_item_t internal_node_remove_child(internal_node_t internal_node, tree_item_t child, tree_t tree)
{
	assert(internal_node != NULL);
	assert(child != NULL);

	tree_item_t tree_item = NULL;

	if (internal_node->item_a == child) {
		internal_node->item_a = NULL;
		tree_branch_destroy(&(internal_node->branch_a), tree);
		tree_item = child;
	} else if (internal_node->item_b == child) {
		internal_node->item_b = NULL;
		tree_branch_destroy(&(internal_node->branch_b), tree);
		tree_item = child;
	} 

	return tree_item;
}


static void tree_add_child(tree_t tree, tree_item_t child)
{
	assert(tree != NULL);
	assert(child != NULL);

	assert(child != tree->item_a);
	assert(child != tree->item_b);
	assert(child != tree->item_c);

	if (tree->item_a == NULL) {
		tree->item_a = child;
		tree->branch_a = tree_branch_create(tree->self, child, tree);
		child->parent_branch = tree->branch_a;
	} else if (tree->item_b == NULL) {
		tree->item_b = child;
		tree->branch_b = tree_branch_create(tree->self, child, tree);
		child->parent_branch = tree->branch_b;
	} else if (tree->item_c == NULL) {
		tree->item_c = child;
		tree->branch_c = tree_branch_create(tree->self, child, tree);
		child->parent_branch = tree->branch_c;
	}  else {
//		fprintf(stderr, "[%p e %p]\n", tree->item_b, tree->item_b);
		fprintf(stderr, "O tree %s nao possui mais espaco (%s)\n", tree->id, tree_item_get_id(child));
		assert(1 == 2);
	}
}


static tree_item_t tree_remove_child(tree_t tree, tree_item_t child)
{
	assert(tree != NULL);
	assert(child != NULL);

	tree_item_t item = NULL;

	if (tree->item_a == child) {
		item = child;
		tree->item_a = NULL;
		tree_branch_destroy(&(tree->branch_a), tree);
	} else if (tree->item_b == child) {
		item = child;
		tree->item_b = NULL;
		tree_branch_destroy(&(tree->branch_b), tree);
	} else if (tree->item_c == child) {
		item = child;
		tree->item_c = NULL;
		tree_branch_destroy(&(tree->branch_c), tree);
	} 

	return item;
}

tree_item_t tree_item_clone(tree_item_t item, tree_item_t parent, tree_t tree)
{
	assert(item != NULL);

	if (tree != NULL) {
		assert(item->type != TREE || parent != NULL);
		assert(item->type != TREE || tree != NULL);
	}

	tree_item_t item_clone = (tree_item_t) malloc(sizeof(struct __tree_item));
	assert(item_clone != NULL);

	item_clone->type       		= item->type;
	item_clone->parent  	    = parent;
	item_clone->parent_distance = item->parent_distance;

	switch (item->type) {
		case LEAF: {
			taxon_node_clone(item->item.taxon_node, tree, item_clone); 
			break;
		}
		
		case INTERNAL_NODE: { 
			internal_node_clone(item->item.internal_node, tree, item_clone); 
			break;
		}

		case TREE: {
			tree_clone(item->item.tree, item_clone); 
			break;
		}

		default: {
			fprintf(stderr, "Item do tipo %d invalido!\n", item->type);
			assert(1 == 2);
		}
	}

	if (parent != NULL) {
		parent_add_child(parent, item_clone, tree);
	}



	return item_clone;
}

tree_t tree_create(values_table_t values_table, taxon_triple_t triple)
{
	assert(triple != NULL);

	double *values;
	tree_t tree = (tree_t) malloc(sizeof(struct __tree));
	assert(tree != NULL);

	tree_item_create(tree, NULL, TREE, tree, 0.0);
	tree->id = create_tree_identifier(0, global_iteration_tree_count);
	tree->iteration  = global_iteration;
	tree->tree_count = global_iteration_tree_count;
	tree->internal_nodes = hash_table_create();
	tree->branchs = hash_table_create();
	tree->taxon_nodes = hash_table_create();
	tree->value = 0.0;
	tree->last_identifier = 0;
	tree->distances = hash_table_create();
	tree->taxon_remains = list_clone(values_table->names->keys);

        tree->matrix = NULL;
        tree->matrix_lines = 0;
        tree->matrix_columns = 0;

        tree->item_a = NULL;
        tree->item_b = NULL;
        tree->item_c = NULL;

	values = calcule_distances(triple->d_ab, triple->d_ac, triple->d_bc);

	taxon_node_create(tree, tree->self, triple->taxon_a, values[A_TO_INTERNAL]);
	tree->branch_a = tree_branch_create(tree->self, tree->item_a, tree);
	
	taxon_node_create(tree, tree->self, triple->taxon_b, values[B_TO_INTERNAL]);
	tree->branch_b = tree_branch_create(tree->self, tree->item_b, tree);

	taxon_node_create(tree, tree->self, triple->taxon_c, values[C_TO_INTERNAL]);
	tree->branch_c = tree_branch_create(tree->self, tree->item_c, tree);

        free(values);
	
	assert(tree->taxon_nodes->keys->size == 3);

	fprintf(stderr, "Arvore %s criada\n", tree->id);

	return tree;
}

char *create_tree_identifier(int iteration, int tree_count)
{
	char *id = (char *) malloc(sizeof(char) * MAX_NAME);
	memset(id, '\0', sizeof(char) * MAX_NAME);

	sprintf(id, "Tree:%d:%d", iteration, tree_count);

	return id;
}

tree_t tree_clone(tree_t tree, tree_item_t self_clone)
{
	assert(tree != NULL);
	assert(self_clone != NULL);

	tree_t tree_clone = (tree_t) malloc(sizeof(struct __tree));
	assert(tree_clone != NULL);

	self_clone->item.tree = tree_clone;

	tree_clone->self  			= self_clone;
	tree_clone->iteration       = global_iteration;
	tree_clone->tree_count      = global_iteration_tree_count;
	tree_clone->previous_tree   = tree;
	tree_clone->id 				= create_tree_identifier(tree_clone->iteration, tree_clone->tree_count);

	tree_clone->distances       = hash_table_clone_all(tree->distances, result_pair_clone);
	tree_clone->branchs         = hash_table_clone_all(tree->branchs, tree_branch_clone);

	tree_clone->internal_nodes  = hash_table_create();
	tree_clone->taxon_nodes     = hash_table_create();
	tree_clone->taxon_remains   = list_clone(tree->taxon_remains);
	tree_clone->last_identifier = tree->last_identifier;
	tree_clone->value           = tree->value;

	tree_clone->item_a = NULL;
	tree_clone->item_b = NULL;
	tree_clone->item_c = NULL;

	tree_clone->item_a = tree_item_clone(tree->item_a, tree_clone->self, tree_clone);
	tree_clone->branch_a = tree_branch_create(tree_clone->self, tree_clone->item_a, tree_clone);

	tree_clone->item_b = tree_item_clone(tree->item_b, tree_clone->self, tree_clone);
	tree_clone->branch_b = tree_branch_create(tree_clone->self, tree_clone->item_b, tree_clone);

	tree_clone->item_c = tree_item_clone(tree->item_c, tree_clone->self, tree_clone);
	tree_clone->branch_c = tree_branch_create(tree_clone->self, tree_clone->item_c, tree_clone);
        
	global_iteration_tree_count++;

	return tree_clone;
}

int tree_destroy(tree_t *tree)
{
	assert(*tree != NULL);

	free((*tree)->id);
	hash_table_destroy_all(&(*tree)->internal_nodes, internal_node_destroy);
	hash_table_destroy_all(&(*tree)->taxon_nodes, taxon_node_destroy);
	hash_table_destroy_all(&(*tree)->branchs, tree_branch_destroy_fast);
	hash_table_destroy_all(&(*tree)->distances, result_pair_destroy);
	list_destroy(&(*tree)->taxon_remains);

	tree_item_destroy(&(*tree)->self);

	if ((*tree)->matrix != NULL) {
		tree_matrix_destroy(&((*tree)->matrix), (*tree)->matrix_lines);
	}
	assert((*tree)->matrix == NULL);
	free(*tree);	*tree = NULL;
	return 1;
}

size_t internal_node_destroy(void *data)
{
	internal_node_t *internal_node = (internal_node_t *) data;
	
	// Nao pode remover as distancias porque elas sao reaproveitadas entre diversas arvores
	hash_table_destroy_all( &(*internal_node)->distances, result_pair_destroy);
	assert((*internal_node)->distances == NULL);
	list_destroy(&(*internal_node)->content);
	tree_item_destroy(&(*internal_node)->self);

	free( (*internal_node)->id );
	free(*internal_node);
	*internal_node = NULL;

	return 1;
}

size_t taxon_node_destroy(void *data)
{
	taxon_node_t *taxon_node = (taxon_node_t *) data;
        //TODO: duplicate at allocation moment
    //    free((*taxon_node)->taxon);
	tree_item_destroy(&(*taxon_node)->self);
	free( *taxon_node );
	*taxon_node = NULL;

	return 1;
}

/**
 * Create a tree_item with in your inner a taxon node
 **/
tree_item_t taxon_node_create(tree_t tree, tree_item_t parent, char *taxon, double parent_distance)
{
	assert(tree != NULL);
	assert(parent != NULL);

	taxon_node_t taxon_tree_node = (taxon_node_t) malloc(sizeof(struct __taxon_node));
	assert(taxon_tree_node != NULL);

	taxon_tree_node->taxon = strdup(taxon);

	tree_item_t tree_item = tree_item_create(tree, parent, LEAF, taxon_tree_node, parent_distance);
	assert(tree_item != NULL);

	hash_table_add(tree->taxon_nodes, taxon_tree_node->taxon, taxon_tree_node);
	list_remove_by_id(tree->taxon_remains, taxon);
	
	assert(tree_item->type == LEAF);
	return tree_item;
}

taxon_node_t taxon_node_clone(taxon_node_t taxon_node, tree_t tree, tree_item_t self_clone)
{
	assert(taxon_node != NULL);
	assert(tree != NULL);
	assert(self_clone != NULL);

	taxon_node_t taxon_clone = (taxon_node_t) malloc(sizeof(struct __taxon_node));
	memset(taxon_clone, '\0', sizeof(struct __taxon_node));
	assert(taxon_clone != NULL);

	self_clone->item.taxon_node = taxon_clone;

	taxon_clone->self  = self_clone;
	taxon_clone->taxon = taxon_node->taxon;
	hash_table_add(tree->taxon_nodes, taxon_clone->taxon, taxon_clone);

	return taxon_clone;
}

/**
 * Create a tree_item with in your inner a internal taxon node
 **/
tree_item_t internal_node_create(tree_t tree, tree_item_t parent, double parent_distance)
{
	assert(tree != NULL);
	assert(parent != NULL);

	internal_node_t internal_tree_node = (internal_node_t) malloc(sizeof(struct __internal_node));
	assert(internal_tree_node != NULL);
	memset(internal_tree_node, '\0', sizeof(struct __internal_node));

	internal_tree_node->id        = create_internal_node_identifier(tree);
	internal_tree_node->distances = hash_table_create();
	internal_tree_node->content   = list_create();

	hash_table_add(tree->internal_nodes, internal_tree_node->id, internal_tree_node);

	tree_item_t tree_item = tree_item_create(tree, parent, INTERNAL_NODE, (void *) internal_tree_node, parent_distance);
	assert(tree_item != NULL);

	assert(parent->type == TREE || parent->type == INTERNAL_NODE);
	internal_tree_node->self->parent_distance = parent_distance;
	tree_update_internal_nodes_distances(internal_tree_node, parent, NULL,  parent_distance);

	assert(tree_item->type == INTERNAL_NODE);

	return tree_item;
}

/**
 * Funcao recursiva que clona os internal nodes da tree e adiciona no tree_clone
 **/
internal_node_t internal_node_clone(internal_node_t internal_node, tree_t tree_clone, tree_item_t self_clone)
{
	assert(internal_node != NULL);
	assert(tree_clone != NULL);
	assert(self_clone != NULL);
	
	internal_node_t internal_clone = (internal_node_t) malloc(sizeof(struct __internal_node));
	assert(internal_clone != NULL);
	memset(internal_clone, '\0', sizeof(struct __internal_node));

	self_clone->item.internal_node = internal_clone;

	internal_clone->self      = self_clone;
	internal_clone->id        = strdup(internal_node->id);
	internal_clone->distances = hash_table_clone_all(internal_node->distances, result_pair_clone);
	internal_clone->content   = list_clone(internal_node->content);

	hash_table_add(tree_clone->internal_nodes, internal_clone->id, internal_clone);

	tree_item_clone(internal_node->item_a, internal_clone->self, tree_clone);
	tree_item_clone(internal_node->item_b, internal_clone->self, tree_clone);

	return internal_clone;
}

static void tree_item_destroy(tree_item_t *tree_item)
{
	assert(*tree_item != NULL);

	if (*tree_item == NULL) {
		return;
	}

	free(*tree_item);
	*tree_item = NULL;
}


static void internal_node_add_content(internal_node_t internal_node, taxon_node_t taxon_node)
{
	assert(internal_node != NULL);
	assert(taxon_node != NULL);

	list_add(internal_node->content, taxon_node->taxon, taxon_node);
	if (internal_node->self->parent->type == INTERNAL_NODE) {
		internal_node_add_content(internal_node->self->parent->item.internal_node, taxon_node);
	}
}

int tree_update_internal_nodes_distances(internal_node_t internal_node_added, tree_item_t actual_node, tree_item_t previous_node, double parent_distance)
{
	assert(internal_node_added != NULL);
	assert(actual_node != NULL);

	int total = 0;

	hash_table_add(internal_node_added->distances, tree_item_get_id(actual_node), 
		result_pair_create(internal_node_added->id, tree_item_get_id(actual_node), parent_distance));

	hash_table_add(tree_item_get_distances(actual_node), internal_node_added->id, 
		result_pair_create(tree_item_get_id(actual_node), internal_node_added->id, parent_distance));


	if ((actual_node->parent != NULL) && (actual_node->parent != previous_node)) {
		total += tree_update_internal_nodes_distances(internal_node_added, actual_node->parent, actual_node, parent_distance + actual_node->parent_distance);
	}

	if (actual_node->type == TREE) {
		total += __tree_update_internal_nodes_distances(actual_node->item.tree->item_a, internal_node_added, actual_node, previous_node, parent_distance);
		total += __tree_update_internal_nodes_distances(actual_node->item.tree->item_b, internal_node_added, actual_node, previous_node, parent_distance);
		total += __tree_update_internal_nodes_distances(actual_node->item.tree->item_c, internal_node_added, actual_node, previous_node, parent_distance);
	} else if (actual_node->type == INTERNAL_NODE) {
		total += __tree_update_internal_nodes_distances(actual_node->item.internal_node->item_a, 
			internal_node_added, actual_node, previous_node, parent_distance);
		total += __tree_update_internal_nodes_distances(actual_node->item.internal_node->item_b, 
			internal_node_added, actual_node, previous_node, parent_distance);
	}

	return total;
}

static int __tree_update_internal_nodes_distances(tree_item_t tree_item, internal_node_t internal_node_added, tree_item_t actual_node, tree_item_t previous_node, double parent_distance)
{
	int total = 0;

	if ((tree_item->type == INTERNAL_NODE) && (tree_item->item.internal_node != internal_node_added) && (tree_item != previous_node)) {
		total++;
		total += tree_update_internal_nodes_distances(internal_node_added, tree_item, actual_node, parent_distance + tree_item->parent_distance);
	}
	
	return total;
}


static tree_branch_t tree_branch_create(tree_item_t item_a, tree_item_t item_b, tree_t tree)
{
	assert(item_a != NULL);
	assert(item_b != NULL);

	tree_branch_t branch = (tree_branch_t) malloc(sizeof(struct __tree_branch));
	CHECK(branch);
	
	branch->id = create_id(get_id(item_a), get_id(item_b), "_<->_");
	branch->item_a = item_a;
	branch->item_b = item_b;
	hash_table_add(tree->branchs, branch->id, branch);
	
	return branch;
}


static void* tree_branch_clone(void *data)
{
	tree_branch_t branch = (tree_branch_t) data;

	// Aqui ta sinistro, porque o item_a e item_b nao sao os mesmos da arvore
	// porem como eh usado apenas dos "pais", eles nao serao destruidos...
	tree_branch_t branch_clone = (tree_branch_t) malloc(sizeof(struct __tree_branch));
	CHECK(branch_clone);

	branch_clone->id = strdup(branch->id);
	branch_clone->item_a = branch->item_a;
	branch_clone->item_b = branch->item_b;

	return branch_clone;
}

static size_t tree_branch_destroy_fast(void *data)
{
	tree_branch_t *tree_branch = (tree_branch_t *) data;

	if(*tree_branch == NULL) {
		return 0;
	}

	free((*tree_branch)->id);
	free(*tree_branch);
	*tree_branch = NULL;

	return 1;
}

static size_t tree_branch_destroy(tree_branch_t *tree_branch, tree_t tree) 
{
	if (*tree_branch == NULL) {
		return 0;
	}

	hash_table_remove_by_id(tree->branchs, (*tree_branch)->id);
	tree_branch_destroy_fast(tree_branch);

	return 1;
}


inline char* create_internal_node_identifier(tree_t tree)
{
	char *id = (char *) malloc(sizeof(char) * MAX_NAME);

	memset(id, '\0', sizeof(char) * MAX_NAME);

	sprintf(id, "Internal Node %d da tree %s(%p)", tree->last_identifier, tree->id, tree);
	tree->last_identifier++;

	return id;
}

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

static char* get_id(tree_item_t item)
{
    assert(item != NULL);

	switch (item->type) {
		case INTERNAL_NODE:
			return item->item.internal_node->id;
			
		case LEAF:
			return item->item.taxon_node->taxon;
			
		case TREE:
			return TREE_ROOT;
			
		default:
			return NULL;
	}
}



/**
 * Compare the values os two trees
 * If are the same, return 0, if tree_one are bigger, return 1,
 * otherwise return -1
 **/
int tree_comparer(const void *one, const void *two)
{
	assert(one != NULL);
	assert(two != NULL);

	if (one == two) {
		return 0;
	}

	tree_t tree_one, tree_two;

	tree_one = (tree_t) one;
	tree_two = (tree_t) two;

	if (tree_one->value == tree_two->value) {
		return 0;
	}

	else if (tree_one->value > tree_two->value) {
		return 1;
	}

	else {
		return -1;
	}

}

static size_t triple_comparer(void *one, void *two)
{
	assert(one != NULL);
	assert(two != NULL);

	if (one == two) {
		return 0;
	}

	taxon_triple_t triple_one, triple_two;

	triple_one = (taxon_triple_t) one;
	triple_two = (taxon_triple_t) two;

	if (triple_one->total_distance == triple_two->total_distance) {
		return 0;
	} 

	else if (triple_one->total_distance > triple_two->total_distance) {
		return 1;
	}

	else {
		return -1;
	}

}


/**
 * Retorna uma lista contendo todas os trios de taxons 
 *  possiveis atraves do values_table passado
 *
 * Os trios sao compostos por inteiros na forma de strings que representam 
 * a posicao do taxon na values_table.
 */
list_t create_triples(values_table_t values_table)
{
	size_t i, j, k, size;
	char *num_1, *num_2, *num_3;
	taxon_triple_t triple = NULL;

	list_t triples_list = list_create();

	size = values_table_get_size(values_table);

	double distance_ab;
	double distance_ac;
	double distance_bc;

	for (i = 0; i < size; i++) {
		num_1 = int_to_string(values_table, i);
		for (j = i + 1; j < size; j++) {
			num_2 = int_to_string(values_table, j);
			distance_ab = values_table_get_value(values_table, i, j);

			for (k = j + 1; k < size; k++) {
				num_3 = int_to_string(values_table, k);

				distance_ac = values_table_get_value(values_table, i, k);
				distance_bc = values_table_get_value(values_table, j, k);

				// Nao se pode criar uma tripla com dois ou mais valores iguais
				if (i == j) {
					fprintf(stderr, "%lu %lu\n", i, j);
					assert(i != j);
				}
				
				if (i == k) {
					fprintf(stderr, "%lu %lu\n", i, k);
					assert(i != k);
				}

				if (j == k) {
					fprintf(stderr, "%lu %lu\n", j, k);
					assert(j != k);
				}

				triple = malloc(sizeof(struct __taxon_triple));

				triple->taxon_a = strdup(num_1);
				triple->taxon_b = strdup(num_2);
				triple->taxon_c = strdup(num_3);

				triple->d_ab = distance_ab;
				triple->d_ac = distance_ac;
				triple->d_bc = distance_bc;

				triple->total_distance = distance_ab + distance_ac + distance_bc;

				list_push(triples_list, num_1, triple);
			}
		}
	}

	return triples_list;
}

static size_t triple_destroy(void *data)
{
	taxon_triple_t* triple = (taxon_triple_t*) data;

        free((*triple)->taxon_a);
        free((*triple)->taxon_b);
        free((*triple)->taxon_c);

	free(*triple);

	return 1;
}

void print_triple(char *id, void* data)
{
	taxon_triple_t triple = (taxon_triple_t) data;
	printf("%s: [%s, %s, %s](%lf , %lf , %lf) total: %lf\n", id, triple->taxon_a, triple->taxon_b, triple->taxon_c, 
		triple->d_ab, triple->d_ac, triple->d_bc, triple->total_distance);
}

void print_tree(char *id, void* data)
{
	tree_t tree = (tree_t) data;
	printf("Tree %s least square %lf\n", id, tree->value);
}

/*Filtra os trios criados, removendo os abaixo */
list_t filter_triples(list_t* list, unsigned int saved)
{
	assert(list != NULL);
	assert(saved > 0);

	unsigned int i;
	cell_t cell;
	taxon_triple_t triple;
	list_t remains_triples = list_create();
	iterator_t iterator = list_iterator(*list);

	for(i = 0; i < saved && iterator->has_next(iterator); i++) {
		cell = iterator->next(iterator);
		triple = cell->data;
		list_push(remains_triples, triple->taxon_a, triple);
	}

        while(iterator->has_next(iterator)) {
            cell = iterator->next(iterator);
            triple = cell->data;
            triple_destroy(&triple);
        }
        
        list_destroy(list);
	list_iterator_destroy(&iterator);

	return remains_triples;
}

/**
 * To create the trees seeds.
 * Return a list with all posibles trees seeds
 * @deprecada
 **/
list_t trees_create(values_table_t values_table, list_t triples)
{
	assert(values_table != NULL);

	iterator_t iterator = list_iterator(triples);

	cell_t cell = NULL;
	taxon_triple_t taxon_triple = NULL;
	list_t trees = list_create();
	tree_t tree = NULL;

	int total = 0;

	/**
	 * Criar todas as arvores iniciais
	 **/
	while (iterator->has_next(iterator)) {
		cell = iterator->next(iterator);

		taxon_triple = cell->data;
		tree = tree_create(values_table, taxon_triple);

		list_push(trees, tree->id, tree);
		total++;
	}
	list_iterator_destroy(&iterator);

	fprintf(stderr, "foram criadas %d arvores\n", total);

	return trees;
}


/**
 * Calcula a distancia de todos entre todos os taxons.
 **/
double calculate_total_distance_tree(tree_t tree)
{
	iterator_t iterator_one, iterator_two;
	cell_t cell_one, cell_two;
	taxon_node_t taxon_node_one, taxon_node_two;
	double total = 0;

	iterator_one = list_iterator(tree->taxon_nodes->keys);
	while (iterator_one->has_next(iterator_one)) {
		cell_one = iterator_one->next(iterator_one);
		taxon_node_one = (taxon_node_t) cell_one->data;

		iterator_two = list_iterator_continues(iterator_one);
		while (iterator_two->has_next(iterator_two)) {
			cell_two = iterator_two->next(iterator_two);
			taxon_node_two = (taxon_node_t) cell_two->data;

			if (taxon_node_one != taxon_node_two) {
				total += calcule_taxons_distance(tree, 
					(taxon_node_t) taxon_node_one, (taxon_node_t) taxon_node_two);
			}
		}
		list_iterator_destroy(&iterator_two);
	}
	list_iterator_destroy(&iterator_one);

	return total;
}

list_t add_taxon_trees(values_table_t values_table, list_t tree_seeds)
{
	cell_t cell_tree, cell_remain, cell_pos;
	iterator_t trees_iterator, positions_iterator;
	taxon_node_t neighbor;
	tree_t tree, new_tree;
	list_t new_trees;
	
	new_trees = list_create();
	trees_iterator = list_iterator(tree_seeds);
	global_iteration_tree_count = 0;
	while(trees_iterator->has_next(trees_iterator)) {
		cell_tree = trees_iterator->next(trees_iterator);
		tree = (tree_t) cell_tree->data;
			
		cell_remain = tree->taxon_remains->first;	
		positions_iterator = list_iterator(tree->taxon_nodes->keys);
		while (positions_iterator->has_next(positions_iterator)) {
			cell_pos = positions_iterator->next(positions_iterator);
			neighbor = (taxon_node_t) cell_pos->data;

			assert(tree != NULL);
			assert(tree->self != NULL);

			new_tree  = tree_item_clone(tree->self, NULL, NULL)->item.tree;
			assert(new_tree->taxon_nodes->keys->size == tree->taxon_nodes->keys->size);

			// pegar o no correto, pois o neighbor pego eh da tree que foi clonada
			neighbor = hash_table_get(new_tree->taxon_nodes, neighbor->taxon);
			assert(neighbor != NULL);
			new_tree = tree_add_taxon(values_table, new_tree, cell_remain->id, neighbor);

			size_t lines, columns;
			new_tree->matrix = tree_create_matrix(new_tree, &lines, &columns);
                        new_tree->matrix_lines = lines;
                        new_tree->matrix_columns = columns;

			new_tree->value = calcule_tree_least_square(values_table, new_tree);
			list_add(new_trees, new_tree->id, new_tree);
			new_tree = NULL;
		}

		list_iterator_destroy(&positions_iterator);
	}

	list_iterator_destroy(&trees_iterator);
	trees_iterator = list_iterator(tree_seeds);
	while (trees_iterator->has_next(trees_iterator)) {
		cell_tree = trees_iterator->next(trees_iterator);
		tree = (tree_t) cell_tree->data;
		tree_destroy(&tree);
	}
	list_destroy(&tree_seeds);
	list_iterator_destroy(&trees_iterator);

	trees_iterator = list_iterator(new_trees);
	list_t rtrees = list_create();
	int total = 0;
	while (trees_iterator->has_next(trees_iterator) && total < 1 ) {
		cell_tree = trees_iterator->next(trees_iterator);
		list_add(rtrees, ((tree_t) cell_tree->data)->id, cell_tree->data);
		total++;
	}

	tree_t tree_to_destroy;
	while (trees_iterator->has_next(trees_iterator)) {
		cell_tree = trees_iterator->next(trees_iterator);
		tree_to_destroy = (tree_t) cell_tree->data;
		tree_destroy(&tree_to_destroy);
	}
	list_iterator_destroy(&trees_iterator);

	return rtrees;
}


/****
             [ Taxon X ]
			     |           /-- [Taxon id]
				 |          /
				 o-----------[ Taxon Node ]  
				 |
				 |
 			 [ Taxon Y ]
*/
		
tree_t tree_add_taxon(values_table_t values_table, tree_t tree, char *taxon_id, taxon_node_t neighbor_node)
{
	assert(values_table != NULL);
	assert(tree != NULL);
	assert(taxon_id != NULL);
	assert(neighbor_node != NULL);

	tree_item_t tree_item = parent_remove_child(neighbor_node->self->parent, neighbor_node->self, tree);
	assert(tree_item != NULL);

	// Posicionar um novo parent para o neighbor
	tree_item_t internal_node = internal_node_create(tree, neighbor_node->self->parent, 0.60);	
	
	// Adicionar o antigo
	neighbor_node->self->parent = internal_node;

	assert(internal_node->item.internal_node->item_a != neighbor_node->self);
	assert(internal_node->item.internal_node->item_b != neighbor_node->self);

	internal_node->item.internal_node->item_a = neighbor_node->self;
	internal_node->item.internal_node->branch_a = tree_branch_create(internal_node, neighbor_node->self, tree);
	neighbor_node->self->parent_branch = internal_node->item.internal_node->branch_a;
	list_add(internal_node->item.internal_node->content, neighbor_node->taxon, neighbor_node);
	
	// Adiciona o nodo taxon_id
	taxon_node_create(tree, internal_node, taxon_id, 1.1);

	return tree;
}

double calcule_tree_least_square(values_table_t values_table, tree_t tree)
{
	assert(values_table != NULL);
	assert(tree != NULL);

	// Populate distances between the internal nodes
	iterator_t taxon_iterator_one = NULL;
	iterator_t taxon_iterator_two = NULL;

	cell_t cell_one;
	cell_t cell_two;

	double taxons_distance;
	double matrix_distance;

	double diference = 0.0;
	double least_square = 0.0;

	taxon_iterator_one = list_iterator(tree->taxon_nodes->keys);
	while (taxon_iterator_one->has_next(taxon_iterator_one)) {
		cell_one = taxon_iterator_one->next(taxon_iterator_one);

		taxon_iterator_two = list_iterator_continues(taxon_iterator_one);
		while (taxon_iterator_two->has_next(taxon_iterator_two)) {
			cell_two = taxon_iterator_two->next(taxon_iterator_two);

			if (cell_one->data != cell_two->data) {
				taxons_distance = calcule_taxons_distance(tree, 
					(taxon_node_t) cell_one->data, (taxon_node_t) cell_two->data);

				matrix_distance = values_table_get_value(values_table, 
						atoi( ((taxon_node_t) cell_one->data)->taxon),
						atoi( ((taxon_node_t) cell_two->data)->taxon)
						);

				diference = taxons_distance - matrix_distance;
				least_square += (diference * diference);
			}
		}
		list_iterator_destroy(&taxon_iterator_two);
	}
	list_iterator_destroy(&taxon_iterator_one);

	return least_square;
}

double calcule_taxons_distance(tree_t tree, taxon_node_t taxon_one, taxon_node_t taxon_two)
{
	assert(tree != NULL);
	assert(taxon_one != NULL);
	assert(taxon_two != NULL);

	tree_item_t parent_one, parent_two;
		
	if (taxon_one->self->parent == taxon_two->self->parent) {
		return taxon_one->self->parent_distance + taxon_two->self->parent_distance;
	} else {
		parent_one = taxon_one->self->parent;
	   	parent_two = taxon_two->self->parent;

		// TODO: algum dia ver porque alguns estao no parent_one e outros no parent_two
		// fprintf(stderr, "%s e %s\n", tree_item_get_id(parent_two), tree_item_get_id(parent_one));

		result_pair_t rp = hash_table_get(tree_item_get_distances(parent_two), tree_item_get_id(parent_one));
		if (rp == NULL) {
			rp = hash_table_get(tree_item_get_distances(parent_one), tree_item_get_id(parent_two));
		}
		assert(rp != NULL);

		return rp->value + taxon_one->self->parent_distance + taxon_two->self->parent_distance;
	}
}

void print_tree_values(list_t trees)
{
	assert(trees != NULL);

	iterator_t iterator;
	cell_t cell;
	tree_t tree;

	iterator = list_iterator(trees);

	while(iterator->has_next(iterator)) {
		cell = iterator->next(iterator);
		tree = cell->data;

		fprintf(stderr, "%s Least Square:%lf\n", tree->id, tree->value);
	}

	list_iterator_destroy(&iterator);
}


static void print_repr_tree_item(values_table_t values_table, tree_item_t tree_item, string_buffer_t sb)
{
	if (tree_item->type == LEAF) {
		string_buffer_add(sb, 
			"%s:%.5lf", 
			values_table_get_name(values_table, atoi(tree_item->item.taxon_node->taxon)), 
			tree_item->parent_distance);

	} else if (tree_item->type == INTERNAL_NODE) {
		print_repr_internal_node(values_table, tree_item->item.internal_node, sb);
		string_buffer_add(sb, ":%.5lf", tree_item->item.internal_node->self->parent_distance);

	} else {
		fprintf(stderr, "tree_item_t (%p) possui um tipo invalido (%d)\n", tree_item, tree_item->type);
	}
}

static void print_repr_internal_node(values_table_t values_table, internal_node_t internal_node, string_buffer_t sb)
{
	assert(internal_node != NULL);

	string_buffer_add(sb, "(");
	print_repr_tree_item(values_table, internal_node->item_a, sb);
	string_buffer_add(sb, ",");
	print_repr_tree_item(values_table, internal_node->item_b, sb);
	string_buffer_add(sb, ")");
}


char *print_reprtree(values_table_t values_table, tree_t tree)
{
	assert(tree != NULL);
	string_buffer_t sb = string_buffer_create();
	
	string_buffer_add(sb, "(");

	print_repr_tree_item(values_table, tree->item_a, sb);		
	string_buffer_add(sb, ",");
	print_repr_tree_item(values_table, tree->item_b, sb);		
	string_buffer_add(sb, ",");
	print_repr_tree_item(values_table, tree->item_c, sb);		

	string_buffer_add(sb, ")");
	string_buffer_add(sb, ";");

	return string_buffer_to_string(sb);
}


list_t trees_construct(values_table_t values_table, list_t trees)
{
	assert(values_table != NULL);
	assert(trees != NULL);

	size_t remain_taxons = values_table_get_size(values_table) - 3;

	while (remain_taxons > 0) {
		fprintf(stderr, "Iteracao %d\n", global_iteration);
		trees = add_taxon_trees(values_table, trees);
		remain_taxons--;
		global_iteration++;
	}

	return trees;
}

#ifdef _LS_RUN_
int main() 
{
        values_table_t values_table = read_dist_file_from_phylip("../data/really_simple_matrix");
	//values_table_t values_table = read_dist_file_from_phylip("../data/tcc_matrix");
	//values_table_t values_table = read_dist_file_from_phylip("../data/simple_matrix");
	//values_table_t values_table = read_dist_file_from_phylip("../data/not_so_simple_matrix");
	//values_table_t values_table = read_dist_file_from_phylip("../data/more_complex_matrix");
	//values_table_t values_table = read_dist_file_from_phylip("../data/49_taxons");
	//values_table_t values_table = read_dist_file_from_paup("../data/domains.NX");
	
        values_table_print(stderr, values_table);
	
        list_t triples = create_triples(values_table);
	fprintf(stderr, "There was created %ld triples... ", triples->size);
	list_t continuing_triples = filter_triples(&triples, 1);
	fprintf(stderr, "keep running %ld triples.\n", continuing_triples->size);

	list_t trees = trees_create(values_table, continuing_triples);
	
	list_destroy_all(&continuing_triples, triple_destroy);

	trees = trees_construct(values_table, trees);

	print_tree_values(trees);

        tree_t tree;
	cell_t cell_tree;
	iterator_t trees_iterator = list_iterator(trees);
	while (trees_iterator->has_next(trees_iterator)) {
		cell_tree = trees_iterator->next(trees_iterator);
		tree = (tree_t) cell_tree->data;
	        char *result = print_reprtree(values_table, tree);
	        fprintf(stderr, "%s", result);
                free(result);
		tree_destroy(&tree);
	}

	list_destroy(&trees);
	list_iterator_destroy(&trees_iterator);
        values_table_destroy(&values_table);
	

	return 0;
}
#endif

