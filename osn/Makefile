
include ../Makefile.inc

CFLAGS += -W -Wall -Wextra 

# Flags specific to building the test app (requires serial compiler)
SERIAL_CC = gcc
PNGINC = 
PNGLIBPATH = 

.PHONY: clean

open-simplex-noise.o:	open-simplex-noise.h open-simplex-noise.c Makefile
	$(CC) $(CFLAGS) -c open-simplex-noise.c

open-simplex-noise-test:	open-simplex-noise-test.c open-simplex-noise.o
	$(SERIAL_CC) $(CFLAGS) $(PNGINC) -o open-simplex-noise-test open-simplex-noise.o open-simplex-noise-test.c $(PNGLIBPATH) -lpng

clean:
	rm -f open-simplex-noise.o open-simplex-noise-test test2d.png test3d.png test4d.png

