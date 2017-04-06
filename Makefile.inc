
# MPI Compiler wrappers preferred; just compilers otherwise

HOSTNAME := $(shell hostname)
ifeq ($(findstring .mil,$(HOSTNAME)),.mil)
  ifeq ($(findstring topaz,$(HOSTNAME)),topaz)
   CC = icc
   FC = ifort
   LIBS = -lmpi
  else
   CC = cc
   FC = ftn
  endif
else
  CC ?= mpicc
  CXX ?= mpicxx
  FC ?= mpif90
endif

# Optimization and debugging flags, (un)comment as desired
OPT = -O3
DBG = -g

# Add HDF5 library, set to install directory
# ** set via the HDF5_DIR environment variable **

# Flags, set as needed; add MPI flags if wrappers weren't available 
CFLAGS =
CXXFLAGS = 
FCFLAGS =
LDFLAGS = 
LIBS = -lm

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

### if ADIOS_HOME is set then enable ADIOS ###
ifdef ADIOS_HOME
  CFLAGS += -DHAS_ADIOS -I$(ADIOS_HOME)/include
  LIBS += -L$(ADIOS_HOME)/lib -ladios -lmxml
endif

### if HDF5_DIR is set then enable HDF5 ###
ifdef HDF5_DIR
  LIBS += $(HDF5_DIR)/lib/libhdf5.a -lz -ldl
  INCLUDE += -I$(HDF5_DIR)/include
  CFLAGS += -DHAS_HDF5
endif

