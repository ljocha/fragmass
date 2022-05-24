#!/usr/bin/env Rscript


suppressMessages(require(Rfast))
suppressMessages(require(rcdk))
suppressMessages(require(purrr))
suppressMessages(require(argparser))
suppressMessages(require(MsBackendMsp))

p <- arg_parser("fragmass.R")
p <- add_argument(p,
	c("--rtdev","--mzdev","--mzrange","input"),
	default=list(0.15,1.,0.15,"fragmassin.msp"),
	help=c("standard deviation in RT (s)", "standard deviation in m/z (ppm)", "range around nominal mass","input file")
)

argv <- parse_args(p, commandArgs(trailingOnly = TRUE)) 

elem <- list(c('C',0,50),c('H',0,50),c('N',0,50),c('O',0,50),c('S',0,50),c('Si',0,10))

sp <- Spectra(argv$input,source=MsBackendMsp())


mmax = 2.5
intlog = 7.5
intdev = 2.
for (i in 1:length(sp)) {
	masses <- c()
	peaks <- sp$mz[i][[1]]
	if (length(peaks) > 0) {
		for (j in 1:length(peaks)) {
			mass <- map( generate.formula(peaks[j], window=argv$mzrange, validation=TRUE, elements=elem), function(x) x@mass)
			masses <- c(masses,unlist(mass))
		}
		rmasses <- round(masses / argv$mzdev * 1000.)
		umasses <- sort_unique(rmasses / 1000. * argv$mzdev)
	
		rt = as.double(sp$RETENTION_TIME[i])
		for (m in umasses) {
			cat(paste(list('unk','unk','unk','unk','unk',
				m,42,14,m-mmax*argv$mzdev,m+mmax*argv$mzdev,
				rt, argv$rtdev, rt-mmax*argv$rtdev, rt+mmax*argv$rtdev,
				intlog, intdev, intlog - intdev, intlog + intdev
			),collapse="	"),'\n')
		}
		}
}

#cat(paste(masses),'\n')
#cat(paste(rmasses),'\n')
#cat(paste(umasses),'\n')





# odhady z dilution series
# mereni je po 0.0005 Da, peaky siroke celkove 0.005, tj. stdev 0.001 
# slabsi peaky 1s siroke, tj. stdev 0.15
# intenzita silne peaky 1G, slabe 20M, sum 500k
