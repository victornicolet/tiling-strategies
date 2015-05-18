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
BENCH_RESULT=profile
VTUNE=amplxe-cl
VTFLAGS=-collect general-exploration -analyze-system
VT_R_DIR=--result-dir $(BENCH_RESULT_DIR)
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
	rm jbi1d

vtune: $(OBJECTS)
	rm -rf $(BENCH_RESULT_DIR)
	$(VTUNE) $(VTFLAGS) $(VT_R_DIR) -- ./$(P_TARGET) $(P_ARGS)
	tar -zcf $(BENCH_RESULT).tar $(BENCH_RESULT)

perf: $(P_TARGET)
	@echo "Scheduler and IPC mechanisms benchmarks .."
	perf bench sched $(P_TARGET) $(P_ARGS)
	@echo "Memory access performance benchmark.."
	perf bench mem $(P_TARGET) $(P_ARGS)

tar:
	tar -zcf $(SOURCES) $(HEADERS) Makefile