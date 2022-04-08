#!/bin/bash

docker run ljocha/fragmass Rscript \
	-e "require(rcdk)" \
	-e "require(purrr)" \
	-e "elem=list(c('C',0,50),c('H',0,50),c('N',0,50),c('O',0,50),c('S',0,50),c('Si',0,10))" \
	-e "mass=map( generate.formula($1, validation=TRUE, elements=elem), function(x) x@mass)" \
	-e "cat(paste(mass),'\n')"
