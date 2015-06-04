ifndef CC
	CC=gcc
endif
ifeq ($(CC),cc)
	CC=gcc
endif

CFLAGS=-g -std=c11 -O3 -Wall
CFLAGS+= -fdiagnostics-color=auto

ifeq ($(CC),icc)
	CFLAGS+=-openmp
else
	CFLAGS+=-fopenmp
endif

# Debuggging and sequential version testing
ifdef DBG
	CFLAGS+=-DDEBUG
endif
ifdef SEQ
	CFLAGS+=-DSEQ
endif

LDFLAGS=-lrt -lm

SRC_DIR=benchmarks
PROGRAM=test
PROGRAMS= $(PROGRAM) jbi1d jbi2d
SOURCES.c=$(SRC_DIR)/jacobi1d.c $(SRC_DIR)/jacobi2d.c utils.c
HEADERS=utils.h $(SRC_DIR)/jacobi1d.h $(SRC_DIR)/jacobi2d.h
OBJECTS=$(SOURCES.c:.c=.o)

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
	P_ARGS=10 000010
endif

#-------------------------------------------------------------------------------

.PHONY: clean

all: $(PROGRAM)

$(PROGRAM) : $(PROGRAM).c $(SOURCES.c) $(HEADERS)
	$(CC) $(CFLAGS) -o $@ $(PROGRAM).c $(LDFLAGS) $(SOURCES.c)

jbi1d : $(SRC_DIR)/jacobi1d.c $(SRC_DIR)/jacobi1d.h
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(SRC_DIR)/jacobi1d.c utils.c

jbi2d : $(SRC_DIR)/jacobi2d.c $(SRC_DIR)/jacobi2d.h
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(SRC_DIR)/jacobi2d.c utils.c

jbi1d_assembly : jbi1d.s

jbi2d_assembly : jbi2d.s

jbi1d.s : $(SRC_DIR)/jacobi1d.c
	$(CC) $(CFLAGS) $(LDFLAGS) -S -o $@ $^

jbi2d.s : $(SRC_DIR)/jacobi2d.c
	$(CC) $(CFLAGS) $(LDFLAGS) -S -o $@ $^

clean:
	rm -f $(PROGRAMS) $(SRC_DIR)/*.o cachegrind.out.* perf.data.* *.s

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
