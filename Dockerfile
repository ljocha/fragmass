FROM rocker/r-ver

RUN apt update 
RUN apt install -y openjdk-17-jre-headless
RUN Rscript -e 'install.packages("rcdk")'
RUN Rscript -e 'install.packages("purrr")'

RUN R CMD javareconf

RUN Rscript -e 'install.packages("BiocManager")' 
RUN Rscript -e 'BiocManager::install("Spectra")'
RUN Rscript -e 'BiocManager::install("MsBackendMsp")'

RUN Rscript -e 'install.packages("argparser")'
RUN Rscript -e 'install.packages("Rfast")'
RUN apt install -y libgsl23
RUN Rscript -e 'install.packages("doParallel")'
RUN Rscript -e 'install.packages("arrow")'
