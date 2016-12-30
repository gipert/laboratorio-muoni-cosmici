# Makefile per la compilazione

CC      = c++ -Wall -O3
INCLUDE = $(shell root-config --cflags)
LIB     = $(shell root-config --libs)
BARPATH = code/ProgressBar/

EXEC = code/analysis/analisiFinale code/montecarlo/semestre2/main code/analysis/lifetimeAnalysis code/TAC/openADC code/montecarlo/semestre2/montecarlo_modifiedforbaseline code/montecarlo/semestre2/baselineStart code/montecarlo/semestre2/expStart code/montecarlo/semestre2/montecarloSingle code/montecarlo/semestre2/montecarlo1exp latex/muon.pdf

all : $(EXEC)

code/analysis/analisiFinale : code/analysis/analisiFinale.cc
	$(CC) $(INCLUDE) -o $@ $^ $(LIB)

code/montecarlo/semestre2/main : code/montecarlo/semestre2/main.cc code/montecarlo/semestre2/montecarlo.cc $(BARPATH)progressbar.cc
	$(CC) $(INCLUDE) -I$(BARPATH) -o $@ $^ $(LIB)

code/montecarlo/semestre2/montecarloSingle : code/montecarlo/semestre2/montecarlo_singleNSim.cc $(BARPATH)progressbar.cc
	$(CC) $(INCLUDE) -I$(BARPATH) -o $@ $^ $(LIB)

code/montecarlo/semestre2/montecarlo1exp : code/montecarlo/semestre2/montecarlo_1exp.cc $(BARPATH)progressbar.cc
	$(CC) $(INCLUDE) -I$(BARPATH) -o $@ $^ $(LIB)

code/montecarlo/semestre2/montecarlo_modifiedforbaseline : code/montecarlo/semestre2/montecarlo_modifiedforbaseline.cc $(BARPATH)progressbar.cc
	$(CC) $(INCLUDE) -I$(BARPATH) -o $@ $^ $(LIB)

code/analysis/lifetimeAnalysis : code/analysis/lifetime_nocalib_bkgr_sub.cc
	$(CC) $(INCLUDE) -o $@ $^ $(LIB)

code/TAC/openADC : code/TAC/openADC.cc
	$(CC) $(INCLUDE) -o $@ $^ $(LIB)

code/montecarlo/semestre2/baselineStart : code/montecarlo/semestre2/baselineStart.cc $(BARPATH)progressbar.cc
	$(CC) $(INCLUDE) -I$(BARPATH) -o $@ $^ $(LIB)

code/montecarlo/semestre2/expStart : code/montecarlo/semestre2/fitexpStart.cc $(BARPATH)progressbar.cc
	$(CC) $(INCLUDE) -I$(BARPATH) -o $@ $^ $(LIB)

latex/muon.pdf : latex/muon.tex latex/imgEff.tex latex/draw.tex latex/analApp.tex latex/Scemi.tex
	cd latex && \
	mkdir -p log && \
	pdflatex -interaction=nonstopmode -output-directory=log muon.tex && \
	cp log/muon.pdf . && \
	cd ..

.PHONY : all clean

clean :
	rm $(EXEC)
