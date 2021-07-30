#!/bin/bash
# sample use ./SARS-COV2_LCRs.sh
#
echo "seg parameters"
echo
read -p "type window length and hit enter: " window
echo
read -p "type locut and hit enter: " locut
echo
read -p "type hicut length and hit enter: " hicut
# INVOKE PIPELINE
ls -1 ../data/Raw_database/ | sed 's/.gb//g' | while read variante; do ./BD_seg.sh $variante $window $locut $hicut ; done
ls -1 ../data/Raw_database/ | sed 's/.gb//g' | while read variante; do  Rscript Analyze_LCRs.R $variante ; done
ls -1 ../data/Raw_database/ | sed 's/.gb//g' | while read variante; do  Rscript Plot_proteome_complexity.R $variante ; done

