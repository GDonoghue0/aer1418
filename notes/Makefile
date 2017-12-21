UNAME = $(shell uname)

PDFLATEX = pdflatex

all: pdf

pdf: 
	$(PDFLATEX) main.tex
	$(PDFLATEX) main.tex

clean:
	rm -f *.bbl *.blg *.ps *.dvi *.aux *.toc \
	      *.idx *.ind *.ilg *.log *.out \
	      *.nav *.snm 
