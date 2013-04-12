pdf : $(patsubst %.eps,%.pdf, $(wildcard *.eps))

%.pdf : %.eps
	epstopdf $<

clean :
	rm -f *.eps *.pdf
