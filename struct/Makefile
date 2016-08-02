
include ../Makefile.inc

OBJS = struct.o
SRCS = struct.c

### Add Output Modules Here ###
OBJS += adiosstruct.o
SRCS += adiosstruct.c

# Example
#OBJS += mymod.o
#SRCS += mymod.c
#CFLAGS += -DHAS_MYMOD

# ADIOS Full Output Module
CFLAGS += -DHAS_ADIOS -I$(HOME)/ace4/adios/include
LIBS += -L$(HOME)/ace4/adios/lib -L$(HOME)/ace4/mxml/lib -ladios -lmxml

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
