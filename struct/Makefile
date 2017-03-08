
include ../Makefile.inc

OBJS = struct.o
SRCS = struct.c

# ADIOS Common
#ENABLE_ADIOS = 1
ifeq ($(ENABLE_ADIOS),1)
  ADIOS_HOME ?= $(HOME)/ace4/adios
  CFLAGS += -DHAS_ADIOS -I$(ADIOS_HOME)/include
  LIBS += -L$(ADIOS_HOME)/lib -ladios -lmxml
  OBJS += adiosstruct.o
  SRCS += adiosstruct.c
endif

### if HDF5_DIR is set then enable HDF5 ###
ifdef HDF5_DIR
  LIBS += $(HDF5_DIR)/lib/libhdf5.a -lz -ldl
  INCLUDE += -I$(HDF5_DIR)/include
  OBJS += hdf5struct.o
  SRCS += hdf5struct.c 
  CFLAGS += -DHAS_HDF5
endif

### End of Output Modules ###

CFLAGS += $(OSNINC)

.PHONY: clean depend

struct: $(OSNOBJ) $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f *.o struct

depend:
	makedepend -- $(CFLAGS) -- -Y $(SRCS)    # ignore warnings about system headers

.c.o:
	$(CC) $(CFLAGS) -c $< $(INCLUDE)

# DO NOT DELETE
