
# MPI Compiler wrappers preferred; just compilers otherwise
ifeq ($(LIBS),)
  LIBS =
endif

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
OPT = -O0
DBG = -ggdb -g3

# Add HDF5 library, set to install directory
# ** set via the HDF5_DIR environment variable **

# Flags, set as needed; add MPI flags if wrappers weren't available
ifeq ($(CFLAGS),)
  CFLAGS =
endif
ifeq ($(CXXFLAGS),)
  CXXFLAGS =
endif
ifeq ($(FCFLAGS),)
  FCFLAGS =
endif
ifeq ($(LDFLAGS),)
  LDFLAGS =
endif
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

### if NC_DIR is set then enable NetCDF-C ###
ifdef NC_DIR
   LIBS += $(NC_DIR)/lib/libnetcdf.a $(NC_DIR)/lib/libhdf5_hl.a $(NC_DIR)/lib/libhdf5.a -lz -ldl
   INCLUDE += -I$(NC_DIR)/include
   CFLAGS += -DHAS_NC
endif

### if HDF5_DIR is set then enable HDF5 ###
ifdef HDF5_DIR
  LIBS += $(HDF5_DIR)/lib/libhdf5.a -lz -ldl
# setenv LD_LIBRARY_PATH $HDF5_DIR/lib
#  LIBS += -L$(HDF5_DIR)/lib/ -lhdf5 -lz -ldl  
  INCLUDE += -I$(HDF5_DIR)/include
  CFLAGS += -DHAS_HDF5
endif
