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

long
skipspace(FILE *afp)
{
	char ch = ' ';
	while (!feof(afp) && isspace(ch)) {
		ch = getc(afp);
	}
	return(ftell(afp));
}

char
*itoa(int value, char *string, int radix)
{
	char            tmp[33];
	char           *tp = tmp;
	int             i;
	unsigned        v;
	int             sign;
	char           *sp;

	if (radix > 36 || radix <= 1) {
		errno = EDOM;
		return 0;
	}
	sign = (radix == 10 && value < 0);
	if (sign)
		v = -value;
	else
		v = (unsigned) value;
	while (v || tp == tmp) {
		i = v % radix;
		v = v / radix;
		if (i < 10)
			*tp++ = i + '0';
		else
			*tp++ = i + 'a' - 10;
	}

	if (string == 0)
		string = (char *) malloc((tp - tmp) + sign + 1);
	sp = string;

	if (sign)
		*sp++ = '-';
	while (tp > tmp)
		*sp++ = *--tp;
	*sp = 0;
	return string;
}


void
Usage(char *program_name)
{
	fprintf(stderr, "\n  Usage:                                                                      \n");
	fprintf(stderr, "  %s [options] <score file list>                                                \n", program_name);
	fprintf(stderr, "  scop2dist takes a list of exhaustive pairwise COMPASS or BLAST score files    \n");
	fprintf(stderr, "  and constructs a NEXUS distance matrix file from them                         \n");
	fprintf(stderr, "                                                                                \n");
	fprintf(stderr, "  Options:                                                                      \n");
	fprintf(stderr, "    -h   help                                                                   \n");
	fprintf(stderr, "    -b   parse BLAST output instead of COMPASS output                           \n");
	fprintf(stderr, "    -e   max E-value to consider as significant-- all scores truncated here     \n");
	fprintf(stderr, "    -tN  transform scores according to method #N (currently 1-4)                \n\n");
	exit(EXIT_SUCCESS);
}


char
*getroot(char *filename)
{
	char           *p = NULL;
	char           *rootname = NULL;

	rootname = malloc((strlen(filename) + 2) * sizeof(char));

	strcpy(rootname, filename);
	p = strrchr(rootname, '.'); /* find where file extension is */
	if (p == NULL) {
		return (rootname);
	} else {
		*p = '\0';
		return (rootname);
	}
}


DISTMAT
*get_distmat(char *listfile_name, int blast)
{
	int           ntax;
	DISTMAT      *distmat = NULL;

	ntax = get_ntax(listfile_name);
	distmat = DISTMATalloc(ntax);
	get_scores(distmat, listfile_name, blast);

	return(distmat);
}


int
get_ntax(char *listfile_name)
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


DISTMAT
*DISTMATalloc(int ntax)
{
	int             i;
	int             name_size = 256;
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
	}

	free(distmat->taxa);
	free(distmat->dist);
	free(distmat->flag);
	free(distmat);
}


DISTMAT
*get_scores(DISTMAT *distmat, char *listfile_name, int blast)
{
	char                 scorefile_appendix[256];
	char                 buff[200];
	int                  i, j;
	double               Evalue;
	double               log_Evalue;
	FILE                *scorefile_ptr = NULL;
	FILE                *listfile_ptr = NULL;
	char                 scorefile_name[256];
	char                *p = NULL;
	char               **nothing = NULL;

	if (blast == 0)
		strcpy(scorefile_appendix, "_compass.out");
	else if (blast == 1)
		strcpy(scorefile_appendix, "_blast.out");

	listfile_ptr = fopen(listfile_name, "r");

	for (i = 0; i < distmat->ntax; ++i) /* for all taxa (down the row) */
	{
		fgets(buff, 256, listfile_ptr);
		sscanf(buff, "%16s", distmat->taxa[i]);
	}

	for (i = 0; i < distmat->ntax; ++i) /* for all taxa (down the row) */
	{
		for (j = 0; j <= i; ++j) /* for each column */
		{
			strcpy(scorefile_name, distmat->taxa[i]);
			strcat(scorefile_name, "_2_");
			strcat(scorefile_name, distmat->taxa[j]);
			strcat(scorefile_name, scorefile_appendix);
			scorefile_ptr = fopen(scorefile_name, "r");

			if (scorefile_ptr == NULL)
			{
				strcpy(scorefile_name, distmat->taxa[j]);
				strcat(scorefile_name, "_2_");
				strcat(scorefile_name, distmat->taxa[i]);
				strcat(scorefile_name, scorefile_appendix);
				scorefile_ptr = fopen(scorefile_name, "r");

				if (scorefile_ptr == NULL) {
					fprintf(stderr,
					        "\n ERROR: scorefile for taxa \"%s\" and \"%s\" not found. \n\n",
					        distmat->taxa[i], distmat->taxa[j]);
					distmat->dist[i][j] = distmat->dist[j][i] = 0.0;
					continue;
				}
			}

			while (fgets(buff, sizeof(buff), scorefile_ptr))
			{
				if (blast == 0)
					p = strstr(buff, "Evalue = ");
				else if (blast == 1)
					p = strstr(buff, "Expect = ");

				/*p = strstr(buff, "Smith-Waterman score = ");*/
				if (p != NULL) {
					Evalue = strtod(p + 9, nothing);
					printf("EValue: %f\n", Evalue);
					log_Evalue = -log10(Evalue);
					/*                     fprintf(stdout, */
					/*                               "%-40.40s  %10.2le  %10.2lf\n", */
					/*                               scorefile_name, Evalue, log_Evalue); */

					distmat->dist[i][j] = distmat->dist[j][i] = Evalue;
					break;
				}
			}

			if (p == NULL)
			{
				fprintf(stderr,
				        "\n ERROR: could not find E-value for file \"%s\". ",
				        scorefile_name);
				Evalue = 100.0;
			}

			fclose(scorefile_ptr);
		}
	}

	fprintf(stdout, "\nDone getting E-values. ");

	fclose(listfile_ptr);
	return (distmat);
}


void
transform_scores(DISTMAT *distmat, int transform, double evalue, double bonferroni)
{
	int            i, j;
	double         selfscore = 0.0;
	/* double         smallest = DBL_MAX; */
	double         largest = -DBL_MAX;

	if (bonferroni == 0.0)
		bonferroni = (double) distmat->ntax;

	switch (transform) {
	case 1: /* use raw E-values */
		break;

	case 2: /* convert E-values to P-values */
	{
		for (i = 0; i < distmat->ntax; ++i)
			for (j = 0; j < distmat->ntax; ++j)
				if (distmat->dist[i][j] > 1e-10)
					distmat->dist[i][j] = 1.0 - exp(-distmat->dist[i][j]);
		break;
	}

	case 3: /* Use raw E-values, but truncate those > 1 to 1 */
	{
		for (i = 0; i < distmat->ntax; ++i)
			for (j = 0; j < distmat->ntax; ++j)
				if (distmat->dist[i][j] > -log(evalue))
					distmat->dist[i][j] = 1.0;
		break;
	}

	case 4: /* -log(-log(E-value)) transform,
				                                   dist=0 estimated from average of transformed self-scores */
	{ /* that means some scores, esp. final self-scores, may be negative */
		for (i = 0; i < distmat->ntax; ++i)
		{
			for (j = 0; j <= i; ++j) {
				/* To correct for database size */
				if (distmat->dist[i][j] > (evalue / bonferroni))
					distmat->flag[i][j] = 1;
				else if (distmat->dist[i][j] <= 0.0) /* shouldn't be negative, but ... */
					distmat->dist[i][j] = -log(-log(DBL_MIN));
				else /* do the infamous transform -- 0.57722 is the digamma(1), expected average score for EVD */
					distmat->dist[i][j] = -log(-log(distmat->dist[i][j]) - 0.57722);
			}
		}

		/*             for (i = 0; i < distmat->ntax; ++i) */
		/*                 for (j = 0; j <= i; ++j) */
		/*                     printf("\n%16.5e", distmat->dist[i][j]); */

		/* find the average of self-scores */
		selfscore = 0.0;
		for (i = 0; i < distmat->ntax; ++i)
			selfscore += distmat->dist[i][i];
		selfscore /= (double) distmat->ntax;

		fprintf(stdout, "\nselfscore = %10.3lf", selfscore);

		/* set all self-scores to zero */
		for (i = 0; i < distmat->ntax; ++i)
			distmat->dist[i][i] = 0.0;

		/* subtract the self-score average constant from each distance */
		for (i = 0; i < distmat->ntax; ++i)
			for (j = 0; j < i; ++j)
				distmat->dist[i][j] -= selfscore;

		/* set all negative distances to zero */
		for (i = 0; i < distmat->ntax; ++i)
			for (j = 0; j < i; ++j)
				if (distmat->dist[i][j] < 0.0)
					distmat->dist[i][j] = 0.0;

		/* find largest distance */
		for (i = 0; i < distmat->ntax; ++i)
			for (j = 0; j < i; ++j)
				if (distmat->dist[i][j] > largest)
					largest = distmat->dist[i][j];

		fprintf(stdout, "\nlargest = %10.3lf\n", largest);

		for (i = 0; i < distmat->ntax; ++i)
			for (j = 0; j < i; ++j)
				if (distmat->flag[i][j] == 1)
					distmat->dist[i][j] = largest;

		/* rescale all values to be between 0 and 1 --
		   this is algorithmically superfluous,
		   but it makes the distance matrix look purty */
		for (i = 0; i < distmat->ntax; ++i)
			for (j = 0; j < i; ++j)
				distmat->dist[i][j] /= largest;

		/* copy to other triangle of matrix */
		for (i = 0; i < distmat->ntax; ++i)
			for (j = 0; j < i; ++j)
				distmat->dist[j][i] = distmat->dist[i][j];

		/*             for (i = 0; i < distmat->ntax; ++i) */
		/*                 for (j = 0; j < distmat->ntax; ++j) */
		/*                     printf("\n%16.5e", distmat->dist[i][j]); */

		break;
	}

	default: {
		fprintf(stderr,
		        "\nERROR5: E-value transform designator '%c' unrecognized \n\n",
		        transform);
		exit(EXIT_FAILURE);
	}
	}
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
	fprintf(NXfile_ptr, "  dimensions ntax=%d;\n", distmat->ntax);
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

#ifdef _SCOP2DIST_
int
main(int argc, char *argv[])
{
	char           *listfile_name = NULL;
	char           *NXfile_name = NULL;
	char           *program_name = NULL;
	int             i, j, option, narguments, nosig, flag;
	int             transform = 6;
	int             blast = 0;
	double          evalue = 0.01, bonferroni = 0.0;
	DISTMAT        *distmat = NULL;

	listfile_name = malloc(256 * sizeof(char));
	NXfile_name   = malloc(256 * sizeof(char));

	program_name = argv[0]; /* save program name for later use */

	while ((option = getopt(argc, argv, "B:be:ht:")) != -1) {
		switch (option) {
		case 'B': /* Bonferroni constant, if different from number of taxa */
			bonferroni = (double) strtod(optarg, NULL);
			break;
		case 'b': /* parse BLAST output, not COMPASS */
			blast = 1;
			break;
		case 'e': /* maximum E-value to call significant, all else greater truncated to max dist */
			evalue = (double) strtod(optarg, NULL);
			break;
		case 'h':
			Usage(program_name);
			exit(EXIT_FAILURE);
			break;
		case 't':
			transform = (int) strtol(optarg, NULL, 10);
			break;
		default: {
			fprintf(stderr,
			        "\nBad option '-%c'\n\n", optopt);
			Usage(program_name);
		}
		}
	}
	narguments = argc - optind; /* number of nonoption args */
	argv += optind; /* now argv is set with first arg = argv[0] */

	if (narguments == 0)
		Usage(program_name);

	strcpy(listfile_name, argv[0]);
	NXfile_name = getroot(listfile_name);
	strcat(NXfile_name, ".NX");

	distmat = get_distmat(listfile_name, blast);
	transform_scores(distmat, transform, evalue, bonferroni);

	/* print out the taxa that have no significant hits to anything */
	flag = 0;
	for (i = 0; i < distmat->ntax; ++i) {
		nosig = 0;
		for (j = 0; j < distmat->ntax; ++j) {
			if (distmat->dist[i][j] < 0.99 && i != j)
				nosig += 1;
		}

		if (nosig == 0) {
			if (flag == 0)
				printf("\n\nThese taxa have no significant hits to anything: ");
			printf("\n%s %3d", distmat->taxa[i], i+1);
			flag = 1;
		}
	}

	/* print out the taxa pairs with zero distances (may be redundant) */
	flag = 0;
	for (i = 0; i < distmat->ntax; ++i) {
		for (j = 0; j < i; ++j) {
			if (distmat->dist[i][j] < 1e-2) {
				if (flag == 0)
					printf("\n\nThese taxa pairs have zero distances (possibly redundant): ");
				printf("\n%s %3d and %s %3d", distmat->taxa[i], i+1, distmat->taxa[j], j+1);
				flag = 1;
			}
		}
	}

	print_NX_distmat(distmat, NXfile_name);

	DISTMATdestroy(distmat);

	return (EXIT_SUCCESS);
}

#endif
