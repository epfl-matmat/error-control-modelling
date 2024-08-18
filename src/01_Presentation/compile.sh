#!/bin/sh

pdflatex presentation.tex && \
	pdflatex presentation.tex && \
	pdflatex presentation.tex && \
	rm presentation.log presentation.aux presentation.nav presentation.out presentation.snm presentation.toc
