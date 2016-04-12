
include ../Makefile.inc

OBJS = struct.o
SRCS = struct.c

### Add Output Modules Here ###

# Example
#OBJS += mymod.o
#SRCS += mymod.c
#CFLAGS += -DHAS_MYMOD

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
