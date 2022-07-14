#!/usr/bin/env Rscript


suppressMessages(require(doParallel))
suppressMessages(require(Rfast))
suppressMessages(require(rcdk))
suppressMessages(require(purrr))
suppressMessages(require(argparser))
suppressMessages(require(MsBackendMsp))

p <- arg_parser("fragmass.R")
p <- add_argument(p,
	c("--rtdev","--mzdev","--mzrange","--known_table","input","output"),
	help=c("standard deviation in RT (s)", "relative standard deviation in m/z", "range around nominal mass (Da)","known features table","input file","output file")
)

argv <- parse_args(p, commandArgs(trailingOnly = TRUE)) 

elem <- list(c('C',0,50),c('H',0,50),c('N',0,50),c('O',0,50),c('S',0,50),c('Si',0,10))

sp <- Spectra(argv$input,source=MsBackendMsp())
sp.mz <- sp$mz
sp.RETENTION_TIME <- sp$RETENTION_TIME

rtdev <- argv$rtdev
mzdev <- argv$mzdev
mzrange <- argv$mzrange

# minima/maximal m/z assumed within mean plus/minus 2.5 * stdev
rel.min.max = 2.5

if (! is.na(argv$known_table)) {
	knt <- arrow::read_parquet(argv$known_table)

	rtdev <- mean(knt$RT_sd)
	mzdev <- mean((knt$mz_max - knt$mz_min)/knt$m.z) / (2 * rel.min.max)
	cat(paste("RT_sd from known table: ",as.character(rtdev),"\n"))
	cat(paste("m/z relative stdev from knwon table: ",as.character(mzdev),"\n"))
}

if (is.na(rtdev)) { rtdev <- 0.15 }
if (is.na(mzdev)) { mzdev <- 0.001 / (2 * rel.min.max) }

if (is.na(mzrange)) { mzrange <- 0.15 }


cores = parallel::detectCores() / 2 # FIXME: HT
cl <- parallel::makeCluster(cores[1])
doParallel::registerDoParallel(cl)

out <- foreach(i=1:length(sp),.combine=rbind,.packages=c('rcdk','purrr','Rfast')) %dopar% {
#for (i in 1:length(sp)) {

	cat(paste("processing ",as.character(i),"\n"))
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
			mass.gen <- map( generate.formula(peaks[j], window=mzrange, validation=TRUE, elements=elem), function(x) x@mass)
			masses <- c(masses,unlist(mass.gen))
		}
		rmasses <- round(masses / mzdev)
		umasses <- sort_unique(rmasses * mzdev)
	
		rt = as.double(sp.RETENTION_TIME[i])
		for (m in umasses) {
			m.z <- c(m.z,m)
			mz_min <- c(mz_min,m * (1. - mzdev*rel.min.max))
			mz_max <- c(mz_max,m * (1. + mzdev*rel.min.max))
			RT_mean <- c(RT_mean,rt)
			RT_min <- c(RT_min,rt-rel.min.max*rtdev)
			RT_max <- c(RT_max,rt+rel.min.max*rtdev)
		}
	}
	
	l <- length(m.z)
	RT_sd <- rep(rtdev,l)
	na <- rep(NA,l)

	data.frame(
		chemical_formula=na,
		HMDB_ID=na,
		KEGG_compound_ID=na,
		mass=na,
		ion.type=na,
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

# z apLCMS mz_max-mz_min vychazi 0.000150, to je podezrele malo

# pomery namerene hmotnosti proti teoreticke v rcx_gc-orbitrap_metabolites jsou v radu 1.000004 max
