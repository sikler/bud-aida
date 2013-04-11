INTER = -interaction batchmode

all : pdf

dvi : report.dvi
pdf : report.pdf

%.pdf : %.tex $(wildcard ../eps/*.pdf)
	pdflatex $(INTER) $<
	pdflatex $(INTER) $<

%.dvi : %.tex $(wildcard ../eps/*.eps)
	latex  $<
	latex  $<
	killall --signal SIGUSR1 xdvi.bin &
