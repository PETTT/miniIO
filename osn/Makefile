
include ../Makefile.inc

CFLAGS += -W -Wall -Wextra 
SERIAL_CC = gcc

.PHONY: clean

open-simplex-noise.o:	open-simplex-noise.h open-simplex-noise.c Makefile
	$(CC) $(CFLAGS) -c open-simplex-noise.c

open-simplex-noise-test:	open-simplex-noise-test.c open-simplex-noise.o
	$(SERIAL_CC) $(CFLAGS) -o open-simplex-noise-test open-simplex-noise.o open-simplex-noise-test.c -lpng

clean:
	rm -f open-simplex-noise.o open-simplex-noise-test test2d.png test3d.png test4d.png

