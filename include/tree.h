#ifndef _TREE_H_
#define _TREE_H_ 1

#include "readdist.h"

typedef struct __taxon_triple *taxon_triple_t;

typedef struct __tree_item     *tree_item_t;
typedef struct __internal_node *internal_node_t;
typedef struct __taxon_node    *taxon_node_t;
typedef struct __tree_branch   *tree_branch_t;
typedef struct __tree          *tree_t;


struct __taxon_triple {
	char *taxon_a;
	char *taxon_b;
	char *taxon_c;

	double d_ab;
	double d_ac;
	double d_bc;

	double total_distance;
};

typedef enum {
	LEAF,
	INTERNAL_NODE,
	TREE
} node_type_t;

struct __tree_item {
	tree_item_t   parent;
	tree_branch_t parent_branch;
	node_type_t   type;
	double        parent_distance;

	union {
		taxon_node_t taxon_node;
		internal_node_t internal_node;
		tree_t tree;
	} item;
};

struct __taxon_node {
	tree_item_t 	self;
	char  	        *taxon;
};

struct __internal_node {
	char  	        *id;
	tree_item_t 	self;

	tree_item_t 	item_a;
	tree_item_t 	item_b;
	tree_branch_t 	branch_a;	
	tree_branch_t 	branch_b;	

	hash_table_t 	distances; // a cache between the internal nodes and roots distances!

	list_t 			content;   // nodes are in this internal node
};

struct __tree {
	tree_item_t 	self;
	char 		*id;
	int		iteration;        // iteration
	int		tree_count;       // tree pos in your iteration
	unsigned int    last_identifier;  // used to generate internal nodes identifier
	double          value;            // least square value

	tree_t          previous_tree;    // tree "parent"
	list_t 		taxon_remains;    // taxons that arent in the tree

	tree_item_t 	item_a;
        tree_item_t 	item_b;
        tree_item_t 	item_c;

        tree_branch_t 	branch_a;
        tree_branch_t 	branch_b;
        tree_branch_t 	branch_c;

	hash_table_t    internal_nodes;   // hash table with all internal nodes and the distances between them
	hash_table_t    taxon_nodes;      // hash table with all taxon nodes and the distances between them
	hash_table_t    branchs;    

	hash_table_t 	distances;

	size_t 		matrix_lines;
	size_t 		matrix_columns;
	double 		**matrix; 		// matrix 
};


struct __tree_branch {
	char *id;
	tree_item_t item_a;
	tree_item_t item_b;
};

char *print_reprtree(values_table_t values_table, tree_t tree);

list_t trees_create(values_table_t values_table, list_t triples);
tree_t tree_create(values_table_t values_table, taxon_triple_t triple);
char *create_tree_identifier(int actual_iteration, int tree_count);
int tree_destroy(tree_t *tree);

tree_item_t internal_node_create(tree_t tree, tree_item_t parent, double parent_distance);
int internal_node_destroy(void *data);

tree_item_t taxon_node_create(tree_t tree, tree_item_t parent, char *identifier, double parent_distance);
int taxon_node_destroy(void *data);


list_t add_taxon_trees(values_table_t values_table, list_t tree_seeds);

tree_item_t tree_item_clone(tree_item_t item, tree_item_t parent, tree_t tree);
taxon_node_t taxon_node_clone(taxon_node_t taxon_node, tree_t tree, tree_item_t item_clone);
internal_node_t internal_node_clone(internal_node_t internal_node, tree_t tree_clone, tree_item_t tree_item_clone);
tree_t tree_clone(tree_t tree, tree_item_t self_clone);
void tree_item_print(tree_item_t tree_item);

list_t create_triples(values_table_t values_table);
list_t filter_triples(list_t *list, unsigned int saved);

tree_item_t tree_node_create();
int tree_comparer(const void *one, const void *two);

list_t trees_evaluate(list_t trees, int max_trees, int min_trees);
double remove_max_and_calcule_new_max(list_t list, int max);

void print_triple(char *id, void* data);

inline char* create_internal_node_identifier(tree_t tree);

double*	calcule_distances(double d_ab, double d_ac, double d_bc);

double calcule_taxons_distance(tree_t tree, taxon_node_t taxon_one, taxon_node_t taxon_two);

double tree_add_taxons(values_table_t values_table, tree_t tree);
int tree_update_internal_nodes_distances(internal_node_t internal_node_added, tree_item_t actual_node, tree_item_t previous_node, double parent_distance);

tree_item_t* get_internal_node_other_children(internal_node_t parent_node, char *neighbor_identifier);
double calculate_total_distance_tree(tree_t tree);

tree_t tree_add_taxon(values_table_t values_table, tree_t tree, char *taxon_id, taxon_node_t neighbor_node);
double node_minimize_distance(values_table_t values_table, tree_t tree, char* neighbor, internal_node_t parent, char *taxon);
double calcule_tree_least_square(values_table_t values_table, tree_t tree);
double visit_neighbor(tree_t tree, internal_node_t neighbor, double acumuled_distance, list_t visiteds);

tree_t rearrange_distances(values_table_t values_table, tree_t tree);

void print_tree_values(list_t trees);

#endif
