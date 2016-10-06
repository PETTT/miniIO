
include ../Makefile.inc

OBJS = amr.o cubes.o
SRCS = amr.c cubes.c

### Add Output Modules Here ###

# VTK Output Module
ENABLE_VTKOUT = 0
ifeq ($(ENABLE_VTKOUT),1)
  OBJS += vtkout.o
  SRCS += vtkout.c
  CFLAGS += -DHAS_VTKOUT
endif

# ADIOS Common
ENABLE_ADIOS = 0
ifeq ($(ENABLE_ADIOS),1)
  ADIOS_HOME ?= $(HOME)/ace4/adios
  MXML_HOME ?= $(HOME)/ace4/mxml
  OBJS += adiosamr.o
  SRCS += adiosamr.c
  CFLAGS += -DHAS_ADIOS -I$(ADIOS_HOME)/include
  LIBS += -L$(ADIOS_HOME)/lib -ladios -L$(MXML_HOME)/lib -lmxml
endif

### if HDF5_DIR set then enable HDF5 ###
ifdef HDF5_DIR
  LIBS += $(HDF5_DIR)/lib/libhdf5.a -lz -ldl
  INCLUDE += -I$(HDF5_DIR)/include
  OBJS += hdf5amr.o
  SRCS += hdf5amr.c
  CFLAGS += -DHAS_HDF5
endif

### End of Output Modules ###

CFLAGS += $(OSNINC)

.PHONY: clean depend

amr: $(OSNOBJ) $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f *.o amr *.h5

depend:
	makedepend -- $(CFLAGS) -- -Y $(SRCS)    # ignore warnings about system headers

.c.o:
	$(CC) $(CFLAGS) -c $< $(INCLUDE)

# DO NOT DELETE
