
include ../Makefile.inc

OBJS = unstruct.o
SRCS = unstruct.c

### Add Output Modules Here ###

# Example
#OBJS += mymod.o
#SRCS += mymod.c
#CFLAGS += -DHAS_MYMOD

### End of Output Modules ###

.PHONY: clean depend

unstruct: $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f *.o unstruct

depend:
	makedepend -- $(CFLAGS) -- -Y $(SRCS)    # ignore warnings about system headers

.c.o:
	$(CC) $(CFLAGS) -c $<

# DO NOT DELETE

