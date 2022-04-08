FROM rocker/r-ver

RUN apt update 
RUN apt install -y openjdk-17-jre-headless
RUN Rscript -e 'install.packages("rcdk")'
RUN Rscript -e 'install.packages("purrr")'

RUN R CMD javareconf
