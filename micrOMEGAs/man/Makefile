
all: manual_6.2.pdf 

manual_6.2.pdf:manual_6.2.tex
	pdflatex  manual_6.2.tex
	bibtex manual_6.2
	pdflatex   manual_6.2.tex
	pdflatex   manual_6.2.tex
clean: 
	rm -f  *.dvi *.log *.aux *.toc *.bbl *.out *.blg .DS_Store