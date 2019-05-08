#FILES = README.md COPYING.md pcubature.c hcubature.c cubature.h clencurt.h vwrapper.h converged.h test.c clencurt_gen.c NEWS.md

# CFLAGS = -pg -O3 -fno-inline-small-functions -Wall -ansi -pedantic
# CFLAGS = -g -Wall -ansi -pedantic
CFLAGS = -O3 #-Wall -ansi -pedantic

all: capture 

capture: caprate.c hcubature.c cubature.h converged.h vwrapper.h
	gcc $(CFLAGS) -o $@ caprate.c hcubature.c -lm
