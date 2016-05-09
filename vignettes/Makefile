# Makefile to use knitr for package vignettes

# put all PDF targets here, separated by spaces
PDFS= glmmADMB.pdf

all: $(PDFS)

clean:
	rm -Rf *.tex *.bbl *.blg *.aux *.out *.log *.spl *tikzDictionary unnamed-* *.dvi *.idv *.css *.4tc *.4ct glmmadmb-figs/ Rplots.pdf *.Rout

## %.pdf: %.Rnw
##	$(R_HOME)/bin/Rscript -e "library(knitr); knit2pdf('$*.Rnw')"
