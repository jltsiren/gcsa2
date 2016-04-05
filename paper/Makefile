SOURCES=$(wildcard *.tex)

all: paper.pdf

paper.pdf: cc-by.pdf lipics-logo-bw.pdf paper.bbl $(SOURCES)
	pdflatex paper
	pdflatex paper

bib:
	pdflatex paper
	bibtex paper

paper.bbl: paper.bib
	pdflatex paper
	bibtex paper

clean:
	rm -f *.aux *.bak *.bbl *.blg *.log *.out *.vtc paper.pdf
