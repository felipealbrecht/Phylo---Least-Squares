void
DISTMATdestroy(DISTMAT *distmat);

DISTMAT
*get_scores(DISTMAT *distmat, char *listfile_name, int blast);

void
print_NX_distmat(DISTMAT *distmat, char *NXfile_name);

void
print_NX_distmat_file(DISTMAT *distmat, FILE* NXfile_ptr);

#endif
