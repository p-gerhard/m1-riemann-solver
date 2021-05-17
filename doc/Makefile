FILE=plot_latex_pgfplot 

all:
	pdflatex $(FILE).tex
	rm $(FILE:=.aux) $(FILE:=.log)
clean:
	rm $(FILE:=.aux) $(FILE:=.log)