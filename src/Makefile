# Least square phylogenies inference
#
# Felipe Albrecht (04/11/2011)

.SUFFIXES:
.SUFFIXES: .c .o .do

CC = gcc

PIH_LIB_INCLUDE=../../Pih-Lib/include
PIH_LIB_C=../../Pih-Lib/src/pih_lib.c
PIH_LIB_O=../../Pih-Lib/src/pih_lib.o
PIH_LIB_O=../../Pih-Lib/src/pih_lib.do

INCLUDE = -I../include -I$(PIH_LIB_INCLUDE)

CFLAGS = -Wall -Wextra $(INCLUDE) -D_LS_RUN_
DEBUGFLAGS = -Wall -Wextra -g -O0 $(INCLUDE) -D_LS_RUN_
LIBS = -lm

NAME = -o
EXECNAME = leastsquares
EXECNAMEDEBUG = leastsquares_debug

STANDALONEFILE=tree.c tree_matrix.c readdist.c  result_pair.c scop2dist.c $(PIH_LIB_C)
STANDALONEOBJS=tree.o tree_matrix.o readdist.o  result_pair.o scop2dist.o $(PIH_LIB_O)
DEBUGOBJS=tree.do tree_matrix.do readdist.do  result_pair.do scop2dist.do $(PIH_LIB_O)

.c.o:;  $(CC) -c -o $@ $(CFLAGS) $< 
.c.do:; $(CC) -c -o $@ $(DEBUGFLAGS) -c $<

CC = gcc


all: $(STANDALONEOBJS)
	$(CC) $(CFLAGS) $(NAME) $(EXECNAME) $(STANDALONEOBJS) $(LIBS)

debug: $(DEBUGOBJS)
	$(CC) $(CFLAGS) $(NAME) $(EXECNAMEDEBUG) $(DEBUGOBJS) $(LIBS)

readdist:
	gcc readdist.c -D_DEBUG_READDIST_ -I../include -I../../Pih-Lib/include/ ../../Pih-Lib/src/pih_lib.c -g -o readdist

result_pair:
	gcc result_pair.c scop2dist.c -D_TEST_RP_ -I../include -I../../Pih-Lib/include/ ../../Pih-Lib/src/pih_lib.c -g -o result_pair -lm

clean: 
	rm -rf $(STANDALONEOBJS) $(DEBUGOBJS)
	rm -rf $(EXECNAME) $(EXECNAMEDEBUG)
	rm -rf cscope.files  cscope.in.out  cscope.out  cscope.po.out
	rm -rf readdist result_pair

