
include ../Makefile.inc

OBJS = unstruct.o
SRCS = unstruct.c

### Add Output Modules Here ###

# PRZM Output Module
OBJS += przm.o
SRCS += przm.c
CFLAGS += -DHAS_PRZM

# ADIOS Output Module
ENABLE_ADIOS = 1
ifeq ($(ENABLE_ADIOS),1)
  ADIOS_HOME ?= $(HOME)/ace4/adios
  CFLAGS += -DHAS_ADIOS -I$(ADIOS_HOME)/include
  LIBS += -L$(ADIOS_HOME)/lib -ladios -lmxml
  OBJS += adiosunstruct.o
  SRCS += adiosunstruct.c
endif

# if HDF5_DIR is set then enable HDF5 Output Module ###
ifdef HDF5_DIR
  LIBS += $(HDF5_DIR)/lib/libhdf5.a -lz -ldl
  INCLUDE += -I$(HDF5_DIR)/include
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

unstruct.o: ../osn/open-simplex-noise.h przm.h adiosunstruct.h
przm.o: ../pdirs.h przm.h
adiosunstruct.o: adiosunstruct.h
hdf5.o: ../pdirs.h

