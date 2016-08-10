SOURCES=$(wildcard *.tex)
FIGURES=gcsa2_text_indexes.pdf gcsa2_graph_dbg.pdf gcsa2_pruned_index.pdf

all: paper.pdf

paper.pdf: paper.bbl $(FIGURES) $(SOURCES)
	pdflatex paper
	pdflatex paper

force: paper.bbl $(FIGURES) $(SOURCES)
	rm -f paper.pdf
	pdflatex paper
	pdflatex paper

bib: $(FIGURES)
	pdflatex paper
	bibtex paper

paper.bbl: $(FIGURES) paper.bib
	pdflatex paper
	bibtex paper

gcsa2_text_indexes.pdf: gcsa2_text_indexes.ipe
	ipetoipe -pdf $<

gcsa2_graph_dbg.pdf: gcsa2_graph_dbg.ipe
	ipetoipe -pdf $<

gcsa2_pruned_index.pdf: gcsa2_pruned_index.ipe
	ipetoipe -pdf $<

clean:
	rm -f *.aux *.bak *.bbl *.blg *.log *.out *.vtc paper.pdf $(FIGURES)
