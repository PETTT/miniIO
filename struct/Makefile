
include ../Makefile.inc

OBJS = struct.o
SRCS = struct.c

### Add Output Modules Here ###
OBJS += adiosstruct.o
SRCS += adiosstruct.c

# Example
#OBJS += mymod.o
#SRCS += mymod.c
#CFLAGS += -DHAS_MYMOD

# ADIOS Full Output Module
CFLAGS += -DHAS_ADIOS -I$(HOME)/ace4/adios/include
LIBS += -L$(HOME)/ace4/adios/lib -L$(HOME)/ace4/mxml/lib -ladios -lmxml

### if HDF5_DIR set then enable HDF5 ###
ifdef HDF5_DIR
  LIBS += $(HDF5_DIR)/lib/libhdf5.a -lz -ldl
  INCLUDE += -I$(HDF5_DIR)/include
  OBJS += hdf5i.o hdf5p.o
  SRCS += hdf5i.c hdf5p.c
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
	$(CC) $(CFLAGS) -c $<

# DO NOT DELETE
