#FILES = README.md COPYING.md pcubature.c hcubature.c cubature.h clencurt.h vwrapper.h converged.h test.c clencurt_gen.c NEWS.md

CFLAGS = -pg -O3 -fno-inline-small-functions -Wall -ansi -pedantic
CFLAGS = -g -Wall -ansi -pedantic
CFLAGS = -O3 -Wall #-ansi -pedantic

all: run_MC

run_MC: MC_cap.c
	gcc $(CFLAGS) -o $@ MC_cap.c -lgsl -lgslcblas -lm
