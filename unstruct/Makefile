
include ../Makefile.inc

OBJS = unstruct.o
SRCS = unstruct.c

### Add Output Modules Here ###

# PRZM Output Module
OBJS += przm.o
SRCS += przm.c
CFLAGS += -DHAS_PRZM

### if HDF5_DIR is set then enable HDF5 ###
ifdef HDF5_DIR
  LIBS = $(HDF5_DIR)/lib/libhdf5.a -lz -ldl
  INCLUDE = -I$(HDF5_DIR)/include
  OBJS += hdf5.o
  SRCS += hdf5.c
  CFLAGS += -DHAS_HDF5
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

unstruct.o: ../osn/open-simplex-noise.h przm.h
przm.o: ../pdirs.h przm.h
