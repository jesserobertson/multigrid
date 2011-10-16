#!/usr/bin/env bash
# Jess Robertson, 2010-11-04

# Generate plots
../python/ncPlot.py                       
                       
# Concatenate files 
mkdir temp_1x35
../python/joinPDF.py --output="temp_1x35/velocity.pdf" *_velocity.pdf
../python/joinPDF.py --output="temp_1x35/log_residual.pdf" *_log_residual.pdf
../python/joinPDF.py --output="temp_1x35/strain_rate.pdf" *_strain_rate.pdf

# Delete other pdf files
rm *pdf                 

# Move all files to current folder and open
mv temp_1x35/*.pdf . && rm -R temp_1x35
open *.pdf