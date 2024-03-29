
include ../Makefile.inc

OBJS = iso.o sfc.o cartiso.o
SRCS = iso.c sfc.c cartiso.c

### Add Output Modules Here ###

# PVTI Full Output Module
ENABLE_PVTI = 1
ifeq ($(ENABLE_PVTI),1)
  OBJS += pvti.o
  SRCS += pvti.c
  CFLAGS += -DHAS_PVTI
endif

# PVTP Isosurface Output Module
ENABLE_PVTP = 1
ifeq ($(ENABLE_PVTP),1)
  OBJS += pvtp.o
  SRCS += pvtp.c
  CFLAGS += -DHAS_PVTP
endif

ifdef ADIOS_HOME
  # ADIOS Full Output Module
  OBJS += adiosfull.o
  SRCS += adiosfull.c
  # ADIOS Isosurface Output Module
  OBJS += adiosiso.o
  SRCS += adiosiso.c
endif

### if NETCDF_DIR set then enable netCDF
ifdef NETCDF_DIR
	OBJS += nci.o ncp.o
	SRCS += nci.c ncp.c
endif

### if HDF5_DIR set then enable HDF5 ###
ifdef HDF5_DIR
  OBJS += hdf5i.o hdf5p.o
  SRCS += hdf5i.c hdf5p.c
endif

### End of Output Modules ###

CFLAGS += $(OSNINC)

.PHONY: clean depend

cartiso: $(OSNOBJ) $(OBJS)
	$(CC)  -o $@ $^ $(LDFLAGS) $(LIBS)

clean:
	rm -f *.o cartiso sfc

depend:
	makedepend -- $(CFLAGS) -- -Y $(SRCS)    # ignore warnings about system headers

# Test version of sfc (for testing proper functioning of sfc code)
sfc: sfc.c
	cc $(CFLAGS) -DSFC_SERPENTINE_TEST -o sfc sfc.c

.c.o:
	$(CC) $(CFLAGS) -c $< $(INCLUDE)

# DO NOT DELETE

iso.o: iso.h
sfc.o: sfc.h
cartiso.o: sfc.h iso.h ../timer.h ../osn/open-simplex-noise.h pvti.h pvtp.h
cartiso.o: adiosfull.h adiosiso.h
pvti.o: ../pdirs.h pvti.h
pvtp.o: ../pdirs.h pvtp.h
adiosfull.o: adiosfull.h ../pdirs.h
adiosiso.o: adiosiso.h
