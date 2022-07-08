#!/usr/bin/env Rscript


suppressMessages(require(doParallel))
suppressMessages(require(Rfast))
suppressMessages(require(rcdk))
suppressMessages(require(purrr))
suppressMessages(require(argparser))
suppressMessages(require(MsBackendMsp))

p <- arg_parser("fragmass.R")
p <- add_argument(p,
	c("--rtdev","--mzdev","--mzrange","input","output"),
	default=list(0.15,1.,0.15,"fragmassin.msp","fragmassout.parquet"),
	help=c("standard deviation in RT (s)", "standard deviation in m/z (ppm)", "range around nominal mass","input file","output file")
)

argv <- parse_args(p, commandArgs(trailingOnly = TRUE)) 

elem <- list(c('C',0,50),c('H',0,50),c('N',0,50),c('O',0,50),c('S',0,50),c('Si',0,10))

sp <- Spectra(argv$input,source=MsBackendMsp())
sp.mz <- sp$mz
sp.RETENTION_TIME <- sp$RETENTION_TIME

mmax = 2.5
intlog = 7.5
intdev = 2.
plus.h.mass = 1.007276

cores = parallel::detectCores() / 2 # FIXME: HT
cl <- parallel::makeCluster(cores[1])
doParallel::registerDoParallel(cl)

out <- foreach(i=1:length(sp),.combine=rbind,.packages=c('rcdk','purrr','Rfast')) %dopar% {
#for (i in 1:length(sp)) {
	mass <- c()
	m.z <- c()
	mz_min <- c()
	mz_max <- c()
	RT_mean <- c()
	RT_min <- c()
	RT_max <- c()
	
	peaks <- sp.mz[i][[1]]
	if (length(peaks) > 0) {
		masses <- c()
		for (j in 1:length(peaks)) {
			mass.gen <- map( generate.formula(peaks[j], window=argv$mzrange, validation=TRUE, elements=elem), function(x) x@mass)
			masses <- c(masses,unlist(mass.gen))
		}
		rmasses <- round(masses / argv$mzdev * 1000.)
		umasses <- sort_unique(rmasses / 1000. * argv$mzdev)
	
		rt = as.double(sp.RETENTION_TIME[i])
		for (m in umasses) {
			mass <- c(mass,m)
			mplus <- m + plus.h.mass
			m.z <- c(m.z,mplus)
			mz_min <- c(mz_min,mplus-mmax*argv$mzdev)
			mz_max <- c(mz_max,mplus+mmax*argv$mzdev)
			RT_mean <- c(RT_mean,rt)
			RT_min <- c(RT_min,rt-mmax*argv$rtdev)
			RT_max <- c(RT_max,rt+mmax*argv$rtdev)
#cat(paste("m:",m,"\n"))
#cat(paste("mass:",length(mass),"\n"))
#cat(paste("m.z:",length(m.z),"\n"))
		}
	}
	
	l <- length(mass)
	ion.type <- rep('M+H',l)
	RT_sd <- rep(argv$rtdev,l)
	na <- rep(NA,l)

#print(length(na))
#print(length(mass))
#print(length(m.z))
#print(length(mz_min))
#print(length(mz_max))
#print(length(RT_mean))
#print(length(RT_min))
#print(length(RT_max))
#print(length(ion.type))
#print(length(RT_sd))
#print("uff")

	data.frame(
		chemical_formula=na,
		HMDB_ID=na,
		KEGG_compound_ID=na,
		mass=mass,
		ion.type=ion.type,
		m.z=m.z,
		Number_profiles_processed=na,
		Percent_found=na,
		mz_min=mz_min,
		mz_max=mz_max,
		RT_mean=RT_mean,
		RT_sd=RT_sd,
		RT_min=RT_min,
		RT_max=RT_max,
		int_mean=na,
		int_sd=na,
		int_min=na,
		int_max=na
	)
}

parallel::stopCluster(cl)

colnames(out) <- c("chemical_formula", "HMDB_ID", "KEGG_compound_ID", "mass", "ion.type",
    "m.z", "Number_profiles_processed", "Percent_found", "mz_min", "mz_max",
    "RT_mean", "RT_sd", "RT_min", "RT_max", "int_mean(log)", "int_sd(log)",
    "int_min(log)", "int_max(log)")

arrow::write_parquet(out,argv$output)



# odhady z dilution series
# mereni je po 0.0005 Da, peaky siroke celkove 0.005, tj. stdev 0.001 
# slabsi peaky 1s siroke, tj. stdev 0.15
# intenzita silne peaky 1G, slabe 20M, sum 500k
