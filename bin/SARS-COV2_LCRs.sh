#!/bin/bash
# sample use ./SARS-COV2_LCRs.sh
#
#
window=12
locut=1.9
hicut=2.1
# INVOKE PIPELINE
ls -1 ../data/Raw_database/ | sed 's/.gb//g' | while read variante; do ./BD_seg.sh $variante $window $locut $hicut ; done
ls -1 ../data/Raw_database/ | sed 's/.gb//g' | while read variante; do  Rscript Analyze_LCRs.R $variante ; done
ls -1 ../data/Raw_database/ | sed 's/.gb//g' | while read variante; do  Rscript Plot_proteome_complexity.R $variante ; done

