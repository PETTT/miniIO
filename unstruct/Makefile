
include ../Makefile.inc

OBJS = unstruct.o
SRCS = unstruct.c

### Add Output Modules Here ###

# PRZM Output Module
OBJS += przm.o
SRCS += przm.c
CFLAGS += -DHAS_PRZM

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
	$(CC) $(CFLAGS) -c $<

# DO NOT DELETE

unstruct.o: ../osn/open-simplex-noise.h przm.h
przm.o: ../pdirs.h przm.h
