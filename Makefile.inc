
# MPI Compiler wrappers preferred; just compilers otherwise
CC = mpicc
CXX = mpicxx
FC = mpif90

# Optimization and debugging flags, (un)comment as desired
OPT = -O3
#DBG = -g

# Add HDF5 library, set to install directory
HDF5_DIR = /scr/brtnfld/trunk/hdf5

# Flags, set as needed; add MPI flags if wrappers weren't available 
CFLAGS = 
CXXFLAGS = 
FCFLAGS =
LDFLAGS =  -lm
LIBS =

# You should probably leave the rest below alone

CFLAGS += $(OPT) $(DBG) -I..
CXXFLAGS += $(OPT) $(DBG) -I..
FCFLAGS += $(OPT) $(DBG) -I..
LDFLAGS += $(OPT) $(DBG) -I..

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

