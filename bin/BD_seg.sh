#!/bin/bash
# sample use ./genomic_complexplot.sh $1 12 1.9 2.1
#
#
# DELETING PREVIOUS INTERMEDIATE FILES, ONLY TO PREVENT WORKFLOW ISSUES
#
rm ../results/Proteomic_fasta/$1_proteins.faa 2> /dev/null
rm ../results/GenFeatures_locations/$1_features_locations.csv 2> /dev/null
rm ../results/seg/$1_complexity.csv 2> /dev/null
rm ../results/seg/$1_high_complexity.faa 2> /dev/null
rm ../results/seg/$1_low_complexity.faa 2> /dev/null
#
# INVOKE PIPELINE
#
perl BD_Genbank_to_peptide_locations.pl $1.gb
perl BD_Gb_2_faa.pl $1.gb
./Just_a_seg_envelope.sh $1 $2 $3 $4
Rscript Filter_LCR_x_aa.R $1
sed -i 's/"//g' ../results/seg/$1_complexity.csv


