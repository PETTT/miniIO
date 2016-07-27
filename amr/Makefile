
include ../Makefile.inc

OBJS = amr.o
SRCS = amr.c

### Add Output Modules Here ###

# Example
#OBJS += mymod.o
#SRCS += mymod.c
#CFLAGS += -DHAS_MYMOD

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

amr: $(OSNOBJ) $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f *.o amr

depend:
	makedepend -- $(CFLAGS) -- -Y $(SRCS)    # ignore warnings about system headers

.c.o:
	$(CC) $(CFLAGS) -c $< $(INCLUDE)

# DO NOT DELETE
