
include ../Makefile.inc

OBJS = struct.o splitspace.o
SRCS = struct.c splitspace.c

ifdef ADIOS_HOME
  OBJS += adiosstruct.o
  SRCS += adiosstruct.c
endif

ifdef HDF5_DIR
  OBJS += hdf5struct.o
  SRCS += hdf5struct.c
endif

ifdef NETCDF_DIR
	OBJS += netcdfstruct.o
	SRCS += netcdfstruct.c
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
