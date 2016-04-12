
include ../Makefile.inc

OBJS = iso.o sfc.o cartiso.o
SRCS = iso.c sfc.c cartiso.c

### Add Output Modules Here ###

# PVTI Full Output Module
OBJS += pvti.o
SRCS += pvti.c
CFLAGS += -DHAS_PVTI

# PVTP Isosurface Output Module
OBJS += pvtp.o
SRCS += pvtp.c
CFLAGS += -DHAS_PVTP

### End of Output Modules ###

CFLAGS += $(OSNINC)

.PHONY: clean depend

cartiso: $(OSNOBJ) $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f *.o cartiso sfc

depend:
	makedepend -- $(CFLAGS) -- -Y $(SRCS)    # ignore warnings about system headers

# Test version of sfc (for testing proper functioning of sfc code)
sfc: sfc.c
	cc $(CFLAGS) -DSFC_SERPENTINE_TEST -o sfc sfc.c

.c.o:
	$(CC) $(CFLAGS) -c $<

# DO NOT DELETE

iso.o: iso.h
sfc.o: sfc.h
cartiso.o: sfc.h iso.h ../timer.h ../osn/open-simplex-noise.h pvti.h pvtp.h
pvti.o: pvti.h
pvtp.o: pvtp.h
