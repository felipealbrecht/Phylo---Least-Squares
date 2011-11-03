#include "pih_lib.h"

#include "tree.h"

double** tree_create_matrix(tree_t tree, size_t *lines, size_t *columns);
void tree_matrix_destroy(double ***matrix, size_t lines);
list_t leaf_walk_to(taxon_node_t leaf_from, taxon_node_t leaf_to);
