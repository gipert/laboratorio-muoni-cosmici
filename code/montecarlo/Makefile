# Makefile per la compilazione

CC      = c++
INCLUDE = $(shell root-config --cflags)
LIB     = $(shell root-config --libs)

EXEC = montecarlo

all : $(EXEC)

montecarlo : montecarlo.cc
	$(CC) $(INCLUDE) -o $@ $< $(LIB)

.PHONY : all clean

clean :
	rm montecarlo
