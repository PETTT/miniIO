
# MPI Compilers
CC = mpicc
CXX = mpicxx
FC = mpif90

# Optimization and debugging flags, (un)comment as desired
OPT = -O3
#DBG = -g

# Flags; leave OPT & DBG alone but add whatever else you want
CFLAGS = $(OPT) $(DBG)
CXXFLAGS = $(OPT) $(DBG)
FCFLAGS = $(OPT) $(DBG)
LDFLAGS = $(OPT) $(DBG)
LIBS = 

# You should probably leave the rest below alone

# Common rules for those needing osn
#   just add $(OSNOBJ) to list of dependencies (and objects if not using $^)
#   and add $(OSNINC) for includes

OSNINC = -I../osn
OSNOBJ = ../osn/open-simplex-noise.o

#.PHONY: $(OSNOBJ)

$(OSNOBJ): ../osn/open-simplex-noise.h ../osn/open-simplex-noise.c
	$(MAKE) -C ../osn open-simplex-noise.o

# Don't let OBNOBJ be the default goal
# (a gnu make-ism, so be careful with weird makes)
.DEFAULT_GOAL :=

