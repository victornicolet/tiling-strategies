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
BENCH_RESULT=profile
VTUNE=amplxe-cl
VTFLAGS=-collect general-exploration -analyze-system
VT_R_DIR=--result-dir $(BENCH_RESULT_DIR)
#Hardware counters for profiling
REPORT_FREQ=99
HW_COUNTERS=L1-dcache-loads,L1-dcache-load-misses,L1-dcache-stores
# Valgrind options
VALGRIND_OPTS+= -q --alignment=16
ifndef P_TARGET
	P_TARGET=jbi1d
endif
ifndef P_ARGS
	P_ARGS=10 0010
endif

.PHONY: clean

LDFLAGS=-lrt

#_______________________________#

jbi1d : jacobi1d.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^

clean:
	rm jbi1d cachegrind.out.*

vtune: $(OBJECTS)
	rm -rf $(BENCH_RESULT_DIR)
	$(VTUNE) $(VTFLAGS) $(VT_R_DIR) -- ./$(P_TARGET) $(P_ARGS)
	tar -zcf $(BENCH_RESULT).tar $(BENCH_RESULT)

perfmem: $(P_TARGET)
	perf record -F $(REPORT_FREQ) -e $(HW_COUNTERS) -a -g \
	-- ./$(P_TARGET) $(P_ARGS)

memcheck: $(P_TARGET)
	valgrind --tool=cachegrind $(VALGRIND_OPTS) ./$(P_TARGET) $(P_ARGS)

viewopts: 
	@echo "\t-- Profiling parameters --"
	@echo "P_TARGET\t" $(P_TARGET) 
	@echo "P_ARGS\t\t" $(P_ARGS)
	@echo "VALGRIND_OPTS \t" $(VALGRIND_OPTS)
	@echo "HW_COUNTERS \t" $(HW_COUNTERS)
tar:
	tar -zcf $(SOURCES) $(HEADERS) Makefile
