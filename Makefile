ifndef CC
	CC=gcc
endif

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
HEADERS=utils.
OBJECTS=jbi1d

# Profiling

#Profiling
BENCH_RESULT=vtune_profile
VTUNE=amplxe-cl
VTFLAGS=-collect general-exploration -analyze-system
VT_R_DIR=--result-dir $(BENCH_RESULT_DIR)

.PHONY: clean

LDFLAGS=-lrt

#_______________________________#

jbi1d : jacobi1d.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^

clean:
	rm jbi1d

vtune: $(OBJECTS)
	rm -rf $(BENCH_RESULT_DIR)
	$(VTUNE) $(VTFLAGS) $(VT_R_DIR) -- ./jbi1d $(IMAGE) $(RUNS)
	tar -zcf $(BENCH_RESULT).tar $(BENCH_RESULT)

tar:
	tar -zcf $(SOURCES) $(HEADERS) Makefile