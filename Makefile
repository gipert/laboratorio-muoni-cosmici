# Makefile per la compilazione

CC      = c++
INCLUDE = $(shell root-config --cflags)
LIB     = $(shell root-config --libs)

EXEC = code/montecarlo/montecarlo code/analysis/lifetimeAnalysis

all : $(EXEC)

code/montecarlo/montecarlo : code/montecarlo/montecarlo.cc
	$(CC) $(INCLUDE) -o $@ $< $(LIB)

code/analysis/lifetimeAnalysis : code/analysis/lifetime_nocalib_bkgr_sub.cc
	$(CC) $(INCLUDE) -o $@ $< $(LIB)

.PHONY : all clean

clean :
	rm code/montecarlo/montecarlo code/analysis/lifetimeAnalysis
