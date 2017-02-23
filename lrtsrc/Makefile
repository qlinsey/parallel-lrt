#
#  Makefile
#  LRT
#
#  Copyright (c) 2013 The University of Sydney. All rights reserved.
#
#
# generate dependencies automatically
#

# command line tools for building tool
CC = gcc
CCX = g++
LDLIBS = -lm -lpthread -lrt

ifdef FAST
CCFLAGS = -Wall -std=c99 -O4 -msse4 -DFAST
else
ifdef DEBUG
CCFLAGS = -Wall -std=c99 -g -DDEBUG_LEVEL=2
else
CCFLAGS = -Wall -std=c99 -D_POSIX_C_SOURCE=199309
endif
endif

LINK.o = $(CC) -g 

ifdef CUDA
CUDA_INSTALL_PATH = /usr/local/cuda4.0
CUDA_INCLUDE = -I. -I$(CUDA_INSTALL_PATH)/include
CUDA_LIB = -L$(CUDA_INSTALL_PATH)/lib64 -lcudart
NVCCFLAGS = $(CUDA_INCLUDE) -ccbin /usr/local/gcc-4.4.6_x86_64/bin/gcc -O4 -arch=sm_13
NVCC = $(CUDA_INSTALL_PATH)/bin/nvcc
CCFLAGS += -fPIC  -DCUDA
ifdef FAST
NVCCFLAGS += -DFAST
endif
LDLIBS += $(CUDA_LIB) -lstdc++
endif

ifdef MPI
MPICC = mpicc
endif

LRT_OBJ=lrt.o io.o naive_lrt.o fast_lrt.o multicore_lrt.o misc.o prefix.o heap.o gpusim_lrt.o 
ifdef CUDA
LRT_OBJ += gpu_lrt.o
endif




VALGRIND = 

# list all source codes source codes
CSOURCES= $(wildcard *.c)

# tests
TESTS= $(wildcard regression/*.dat)

ifdef CUDA
CUDATESTS = $(TESTS:.dat=.txt.4)
else
CUDATESTS =
endif

# executable 
EXES = lrt mk_grid 

ifdef MPI
EXES += mpi_lrt
MPITESTS = $(TESTS:.dat=.txt.5)
else
MPITESTS = 
endif

#####################################################################################
# phony targets

.PHONY: all clean gentest test

# create all executables
all: $(EXES) 

# clean everything 
clean: 
	rm -f *.o *.d $(EXES)
	rm -f $(TESTS:.dat=.txt.0) $(TESTS:.dat=.txt.1) $(TESTS:.dat=.txt.2) $(TESTS:.dat=.txt.3) $(CUDATESTS) $(MPITESTS)

# generate reference files for testcases in regression 
# using the naive lrt version
gentest: $(TESTS:.dat=.txt) 

# run for various methods (suffix is method number)
test: $(TESTS:.dat=.txt.0) $(TESTS:.dat=.txt.1) $(TESTS:.dat=.txt.2)  $(TESTS:.dat=.txt.3) $(CUDATESTS) $(MPITESTS)

#####################################################################################
# Rules
#

# generates object files from C files (depend on dependency files that are automatically generated)
%.o: %.c %.d
	$(CC) $(CCFLAGS) -c $< -o $@

ifdef CUDA
# generates object files from cuda files
%.o: %.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@
endif


# generate dependency files from C files 
%.d: %.c
	@echo "create dependence for $<";$(SHELL) -ec '$(CC) -MG -MM $(CCFLAGS) $< \
	| sed '\''s/\($*\)\.o[ :]*/\1.o $@ : /g'\'' > $@; \
	[ -s $@ ] || rm -f $@'

#####################################################################################
# Dependencies 

# include all dependency files 
include $(CSOURCES:.c=.d)

# dependencies of executables
lrt: $(LRT_OBJ)

mk_grid: mk_grid.o misc.o io.o

mpi_lrt: misc.o prefix.o heap.o io.o mpi_lrt.d
	$(MPICC) -o mpi_lrt mpi_lrt.c misc.o prefix.o heap.o io.o -lm 

#####################################################################################
# Tests
 
# generate text (result) files from grid data files (used for generation) 
%.txt: %.dat lrt
	@echo "Generate reference file for test $<"
	@./lrt -m 0 -o $@ $<

# run regression test for naive method
%.txt.0: %.dat lrt
	@echo "Test $< using naive lrt"
	@./lrt -m 0 -o $@ $<
	@cmp $@ $(subst .dat,.txt,$<)

# run regression test for fast method
%.txt.1: %.dat lrt
	@echo "Test $< using fast lrt"
	@$(VALGRIND) ./lrt -m 1 -o $@ $<
	@cmp $@ $(subst .dat,.txt,$<)

# run regression test for multicore method
%.txt.2: %.dat lrt
	@echo "Test $< using multicore lrt"
	@$(VALGRIND) ./lrt -m 2 -o $@ $<
	@cmp $@ $(subst .dat,.txt,$<)

# run regression test for gpu simulator method
%.txt.3: %.dat lrt
	@echo "Test $< using gpusim lrt"
	@$(VALGRIND) ./lrt -m 3 -o $@ $<
	@cmp $@ $(subst .dat,.txt,$<)


# run regression test for gpu method
%.txt.4: %.dat lrt
	@echo "Test $< using gpu lrt"
	@$(VALGRIND) ./lrt -m 4 -o $@ $<
	@cut -f1 -d' ' $@ >1.txt
	@cut -f1 -d' ' $(subst .txt.4,.txt,$@) >2.txt
	@cmp 1.txt 2.txt

# run regression test for mpi method
%.txt.5: %.dat lrt
	@echo "Test $< using mpi lrt"
	@$(VALGRIND) mpirun -np 4 ./mpi_lrt -o $@ $<
	@cmp $@ $(subst .dat,.txt,$<)
