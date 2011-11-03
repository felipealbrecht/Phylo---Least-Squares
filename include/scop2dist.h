/***************************************************************************
* scop2dist.h
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
 
# scop2dist.h:
# Header of read compass or BLAST output and create a matrix in PAUP* format
 
****************************************************************/

#ifndef _SCOP2DIST_H_
#define _SCOP2DIST_H_ 1

typedef struct
{
	int             ntax; /* number of taxa */
	char          **taxa; /* pointer to array of taxa names */
	double        **dist; /* pointer to array of distances */
	int           **flag; /* pointer to array of flags for each pairwise dist */
}
DISTMAT;

void
Usage(char *program_name);

char
*getroot(char *filename);

DISTMAT
*get_distmat(char *listfile_name, int blast);

int
get_ntax(char *listfile_name);

DISTMAT
*DISTMATalloc(int ntax);

void
DISTMATdestroy(DISTMAT *distmat);

DISTMAT
*get_scores(DISTMAT *distmat, char *listfile_name, int blast);

void
transform_scores(DISTMAT *distmat, int transform, double evalue, double bonferroni);

void
print_NX_distmat(DISTMAT *distmat, char *NXfile_name);

#endif
