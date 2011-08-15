LATEXFILES = *.aux\
        *.bbl\
        *.blg\
        *.ilg\
        *.log\
        *.nlo\
        *.nls\
        *.toc\
        *.aux\
        mzR.tex\
        Rplots.pdf\
	mzR-plotspectrum.pdf\
        *.dvi\
        *.map\
        *.out\

all:
	make vignette
	make clean

setvars:
ifeq (${R_HOME},)
R_HOME= $(shell R RHOME)
endif

vignette:
	"$(R_HOME)/bin/R" CMD Sweave mzR.Rnw
	pdflatex mzR.tex
	bibtex mzR
	pdflatex mzR.tex
	pdflatex mzR.tex

clean:
	rm -f $(LATEXFILES)
	rm -f *~

