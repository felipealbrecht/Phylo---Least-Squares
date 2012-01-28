/***************************************************************************
* scop2dist.c
* (C) 2005 Douglas L. Theobald
 
* Problems in this file report to Felipe Albrecht (felipe.albrecht@gmail.com)

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
 
# scop2dist.c:
# Read compass or BLAST output and create a matrix in PAUP* format
 
****************************************************************/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#ifdef __linux__
#include <getopt.h>
#endif
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <errno.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include "scop2dist.h"

static long skipspace(FILE *afp);

static long skipspace(FILE *afp)
{
	char ch = ' ';
	while (!feof(afp) && isspace(ch)) {
		ch = getc(afp);
	}
	return(ftell(afp));
}


size_t get_ntax(char *listfile_name)
{
	FILE          *listfile = NULL;
	char           ch = 'A';
	int            ntax = 0;

	listfile = fopen(listfile_name, "r");
	if (listfile == NULL) {
		fprintf(stderr, "\n file \"%s\" does not exist. \n\n", listfile_name);
		return (EXIT_FAILURE);
	}

	while(!feof(listfile)) {
		if (isspace(ch)) {
			skipspace(listfile);
			++ntax;
		}
		ch = getc(listfile);
	}

	fclose(listfile);
	return (ntax);
}


DISTMAT *DISTMATalloc(size_t ntax)
{
	int             i;
	static DISTMAT *distmat = NULL;

	distmat = (DISTMAT *) malloc(sizeof(DISTMAT));

	distmat->taxa = (char **)   malloc(ntax * sizeof(char *));
	distmat->dist = (double **) malloc(ntax * sizeof(double *));
	distmat->flag = (int **)    malloc(ntax * sizeof(int *));

	for (i = 0; i < ntax; ++i) {
		distmat->dist[i] = (double *) calloc(ntax,       sizeof(double));
		distmat->flag[i] = (int *)    calloc(ntax,       sizeof(int));
	}

	distmat->ntax = ntax;
	return(distmat);
}


void
DISTMATdestroy(DISTMAT *distmat)
{
	int             i;

	for (i = 0; i < distmat->ntax; ++i) {
		free(distmat->dist[i]);
		free(distmat->flag[i]);
		free(distmat->taxa[i]);
	}

	free(distmat->taxa);
	free(distmat->dist);
	free(distmat->flag);
	free(distmat);
}


void
print_NX_distmat(DISTMAT *distmat, char *NXfile_name)
{
	FILE          *NXfile_ptr = NULL;

        NXfile_ptr  = fopen(NXfile_name, "w");
        print_NX_distmat_file(distmat, NXfile_ptr);
	fclose(NXfile_ptr);
}


void
print_NX_distmat_file(DISTMAT *distmat, FILE* NXfile_ptr)
{
	int            i, j;

	fprintf(NXfile_ptr, "#NEXUS\n\n");
	fprintf(NXfile_ptr, "begin taxa;\n");
	fprintf(NXfile_ptr, "  dimensions ntax=%lu;\n", distmat->ntax);
	fprintf(NXfile_ptr, "  taxlabels\n");

	for (i = 0; i < distmat->ntax; ++i)
		fprintf(NXfile_ptr, "    %-s\n", distmat->taxa[i]);

	fprintf(NXfile_ptr, "  ;\nend;\n\n");

	fprintf(NXfile_ptr, "begin distances;\n");
	fprintf(NXfile_ptr, "  format\n");
	fprintf(NXfile_ptr, "    triangle=lower;\n");
	/*fprintf(NXfile_ptr, "    missing=?;\n");*/
	fprintf(NXfile_ptr, "  matrix\n");

	fprintf(NXfile_ptr, "               [");
	for (i = 0; i < distmat->ntax; ++i)
		fprintf(NXfile_ptr, " %8.8s", distmat->taxa[i]);
	fprintf(NXfile_ptr, "]\n");

	for (i = 0; i < distmat->ntax; ++i) /* for all taxa (down the row) */
	{
		fprintf(NXfile_ptr, "    %-12.12s", distmat->taxa[i]);
		for (j = 0; j <= i; ++j) /* for each column from 0 up to the diagonal */
		{
			fprintf(NXfile_ptr, " %8.6f", distmat->dist[i][j]);
		}

		fprintf(NXfile_ptr, "\n"); /* end the row */
	}

	fprintf(NXfile_ptr, "               [");
	for (i = 0; i < distmat->ntax; ++i)
		fprintf(NXfile_ptr, " %8.8s", distmat->taxa[i]);
	fprintf(NXfile_ptr, "]\n");

	fprintf(NXfile_ptr, ";\nend;\n\n");
}
