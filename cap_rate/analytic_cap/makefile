#FILES = README.md COPYING.md pcubature.c hcubature.c cubature.h clencurt.h vwrapper.h converged.h test.c clencurt_gen.c NEWS.md

# CFLAGS = -pg -O3 -fno-inline-small-functions -Wall -ansi -pedantic
# CFLAGS = -g -Wall -ansi -pedantic
# CFLAGS = -O3 -Wall #-ansi -pedantic

CFLAGS = -O3

LDFLAGS = -lgsl -lgslcblas -lm

SOURCES=analytic_approx.c NSinterp.c cap_funcs.c
OBJECTS=$(SOURCES:.c=.o)
EXEC=run_analytic


all: $(EXEC)

$(EXEC): $(OBJECTS)
	$(CC) $(CFLAGS)  -o $@ $^  $(LDFLAGS)
	rm *.o

.o:
	$(CC) $(CFLAGS) $@ -c $<

clean:
	rm *.o $(EXEC)
