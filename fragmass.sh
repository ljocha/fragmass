#!/bin/bash

sed 's/COMPOUND_NAME/NAME/' "$1" >/tmp/fragmass_$$.msp

docker run -u $(id -u) -v /tmp:/tmp -v $PWD:/R -w /tmp ljocha/fragmass Rscript /R/fragmass.R /tmp/fragmass_$$.msp

rm -f /tmp/fragmass_$$.msp
