
include ../Makefile.inc

OBJS = amr.o cubes.o
SRCS = amr.c cubes.c

### Add Output Modules Here ###

# VTK Output Module
ENABLE_VTKOUT = 1
ifeq ($(ENABLE_VTKOUT),1)
  OBJS += vtkout.o
  SRCS += vtkout.c
  CFLAGS += -DHAS_VTKOUT
endif

ifdef ADIOS_HOME
  OBJS += adiosamr.o
  SRCS += adiosamr.c
endif

ifdef HDF5_DIR
  OBJS += hdf5amr.o
  SRCS += hdf5amr.c
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
