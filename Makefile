# Makefile per la compilazione

CC      = c++ -Wall -O3
INCLUDE = $(shell root-config --cflags)
LIB     = $(shell root-config --libs)

EXEC = code/montecarlo/semestre2/montecarlo code/analysis/lifetimeAnalysis code/calibration/openADC code/montecarlo/semestre2/montecarlo_modifiedforbaseline code/montecarlo/semestre2/baselineStart

all : $(EXEC)

code/montecarlo/semestre2/montecarlo : code/montecarlo/semestre2/montecarlo.cc
	$(CC) $(INCLUDE) -o $@ $< $(LIB)

code/montecarlo/semestre2/montecarlo_modifiedforbaseline : code/montecarlo/semestre2/montecarlo_modifiedforbaseline.cc
	$(CC) $(INCLUDE) -o $@ $< $(LIB)

code/analysis/lifetimeAnalysis : code/analysis/lifetime_nocalib_bkgr_sub.cc
	$(CC) $(INCLUDE) -o $@ $< $(LIB)

code/calibration/openADC : code/calibration/openADC.cc
	$(CC) $(INCLUDE) -o $@ $< $(LIB)

code/montecarlo/semestre2/baselineStart : code/montecarlo/semestre2/baselineStart.cc
	$(CC) $(INCLUDE) -o $@ $< $(LIB)

.PHONY : all clean

clean :
	rm $(EXEC)
