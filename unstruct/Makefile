
include ../Makefile.inc

OBJS = unstruct.o
SRCS = unstruct.c

### Add Output Modules Here ###

# PRZM Output Module
OBJS += przm.o
SRCS += przm.c
CFLAGS += -DHAS_PRZM

# ADIOS Output Module
ifdef ADIOS_HOME
  OBJS += adiosunstruct.o
  SRCS += adiosunstruct.c
endif

ifdef HDF5_DIR
  OBJS += hdf5unstruct.o
  SRCS += hdf5unstruct.c
endif

### End of Output Modules ###

CFLAGS += $(OSNINC)

.PHONY: clean depend

unstruct: $(OSNOBJ) $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS) 

clean:
	rm -f *.o unstruct

depend:
	makedepend -- $(CFLAGS) -- -Y $(SRCS)    # ignore warnings about system headers

.c.o:
	$(CC) $(CFLAGS) -c $< $(INCLUDE)

# DO NOT DELETE

unstruct.o: ../osn/open-simplex-noise.h przm.h adiosunstruct.h
przm.o: ../pdirs.h przm.h
adiosunstruct.o: adiosunstruct.h
hdf5unstruct.o: ../pdirs.h

