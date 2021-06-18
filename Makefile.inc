
# MPI Compiler wrappers preferred; just compilers otherwise
ifeq ($(LIBS),)
  LIBS =
endif

CC ?= mpicc

# Optimization and debugging flags, (un)comment as desired
OPT ?= -O2
DBG ?= -g

# Add HDF5 library, set to install directory
# ** set via the HDF5_DIR environment variable **

# Flags, set as needed; add MPI flags if wrappers weren't available
ifeq ($(CFLAGS),)
  CFLAGS =
endif
ifeq ($(LDFLAGS),)
  LDFLAGS =
endif
LIBS += -lm

# You should probably leave the rest below alone

CFLAGS += $(OPT) $(DBG) -I..
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
  CFLAGS += -DHAS_ADIOS $(shell $(ADIOS_HOME)/bin/adios_config -c)
  LIBS += $(shell $(ADIOS_HOME)/bin/adios_config -l)
endif

### if NC_DIR is set then enable NetCDF-C ###
ifdef NC_DIR
   LIBS += $(NC_DIR)/lib/libnetcdf.a $(NC_DIR)/lib/libhdf5_hl.a $(NC_DIR)/lib/libhdf5.a -lz -ldl
   INCLUDE += -I$(NC_DIR)/include
   CFLAGS += -DHAS_NC
endif

### if HDF5_DIR is set then enable HDF5 ###
ifdef HDF5_DIR

  H5PCC = $(HDF5_DIR)/bin/h5pcc

  TESTING = $(shell if [ -f /usr/bin/h5pcc.mpich ]; then echo "yes"; else echo "no"; fi)
  
  ifeq ($(TESTING),yes)
    H5PCC = /usr/bin/h5pcc.mpich
  endif

  ## PROCESS h5pcc to get linking options ##

  LINK_LIBS = $(shell $(H5PCC) -showconfig | grep "AM_LDFLAGS:" | sed -n -e 's/^.*AM_LDFLAGS: //p')
  LINK_LIBS += $(shell $(H5PCC) -showconfig | grep "Extra libraries:" | sed -n -e 's/^.*Extra libraries: //p')

  LIBDIR = $(shell if [ -d $$HDF5_DIR/lib ]; then echo "$$HDF5_DIR/lib"; else echo "$$HDF5_DIR/lib64"; fi)

  LIBS += $(LIBDIR)/libhdf5.a $(LINK_LIBS)
  INCLUDE += -I$(HDF5_DIR)/include
  CFLAGS += -DHAS_HDF5
endif
