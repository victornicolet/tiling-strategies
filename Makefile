ifndef CC
	CC=gcc
endif

CFLAGS=-g -std=c11
CFLAGS+=-Waddress -Wstrict-aliasing -Wimplicit-function-declaration -Wformat=2

ifeq ($(CC),icc)
	CFLAGS+=-openmp
else
	ifeq ($(CC),gcc)
		CFLAGS+=-fopenmp -Wopenmp-simd
	else
		WMSG="Compiler unsupported : CC =" $(CC)
		WMSG+="\nYou might want to set CC to gcc or icc"
	endif
endif

CFLAGS+=-O3
LDFLAGS=-lrt -lm

SOURCES=jacobi1d.c test.c
HEADERS=utils.
OBJECTS=jbi1d

# Profiling options ------------------------------------------------------------

BENCH_RESULT=profile
BENCH_RESULT_DIR=$(BENCH_RESULT)/vtune/
# VTune options
VTUNE=amplxe-cl
VTFLAGS=-collect general-exploration -analyze-system
VT_R_DIR=--result-dir $(BENCH_RESULT_DIR)

# Hardware counters for profiling
REPORT_FREQ=99
HW_COUNTERS=L1-dcache-loads,L1-dcache-load-misses,L1-dcache-stores
HW_COUNTERS+=,cache-misses

# Valgrind options
VALGRIND_OPTS+= -q

# Profiling target and application arguments
ifndef P_TARGET
	P_TARGET=jbi1d
endif
ifndef P_ARGS
	P_ARGS=10 0001
endif

#-------------------------------------------------------------------------------

.PHONY: clean

all: $(OBJECTS)

jbi1d : jacobi1d.c
	@echo $(WMSG)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^

clean:
	rm jbi1d cachegrind.out.* perf.data.*

vtune: $(OBJECTS)
	rm -rf $(BENCH_RESULT_DIR)
	$(VTUNE) $(VTFLAGS) $(VT_R_DIR) -- ./$(P_TARGET) $(P_ARGS)
	tar -zcf $(BENCH_RESULT).tar $(BENCH_RESULT)

perfmem: $(P_TARGET)
	perf mem record -- ./$(P_TARGET) $(P_ARGS)

reportmem: perfmem
	perf mem report --sort=mem

memcheck: $(P_TARGET)
	valgrind --tool=cachegrind $(VALGRIND_OPTS) ./$(P_TARGET) $(P_ARGS)

viewopts: 
	@echo "\t-- Profiling parameters --"
	@echo "BENCH_RESULT_DIR \t" $(BENCH_RESULT_DIR)
	@echo "P_TARGET\t" $(P_TARGET) 
	@echo "P_ARGS\t\t" $(P_ARGS)
	@echo "VALGRIND_OPTS \t" $(VALGRIND_OPTS)
	@echo "HW_COUNTERS \t" $(HW_COUNTERS)
tar:
	tar -zcf $(SOURCES) $(HEADERS) Makefile
