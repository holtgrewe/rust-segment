CC = gcc
AR = ar
RANLIB = ranlib

haarseglib_default_libs = -lm

CFLAGS = -g -O2
LIBS = $(haarseglib_default_libs)

.PHONY: clean

default: lib-static

clean:
	rm -f libhaarseg.a HaarSeg.o

lib-static: libhaarseg.a

libhaarseg.a: HaarSeg.o
	@-rm -f $@
	$(AR) -rc $@ $^
	-$(RANLIB) $@

HaarSeg.o: pkg/HaarSeg/src/HaarSeg.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $< $(LIBS)
