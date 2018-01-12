UNAME = $(shell uname)

PDFLATEX = pdflatex

all: pdf

pdf: 
	$(PDFLATEX) main.tex
	$(PDFLATEX) main.tex
	cp main.pdf assignment1.pdf

clean:
	rm -f *.bbl *.blg *.ps *.dvi *.aux *.toc \
	      *.idx *.ind *.ilg *.log *.out \
	      *.nav *.snm 
