CC=gcc

CFLAGS=-g -std=c11
CFLAGS+=-Waddress -Wstrict-aliasing -Wopenmp-simd -Wparentheses
CFLAGS+=  -Wimplicit-function-declaration -Wformat=2

ifeq ($(CC),icc)
	CFLAGS+=-openmp
else
	ifeq ($(CC),gcc)
		CFLAGS+=-fopenmp
	endif
endif

CFLAGS+=-O3

SOURCES=jacobi1d.c test.c
HEADERS=utils.h

.PHONY: clean

LDFLAGS=-lrt

#_______________________________#

jbi1d : jacobi1d.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^

clean:
	rm jbi1d

tar:
	tar -zcf $(SOURCES) $(HEADERS) Makefile