#!/bin/bash
# use ./Just_a_seg_envelope.sh Nidovirales 12 1.9 2.1
# output in ../results/seg/
mkdir -p ../results/seg/
seg ../results/Proteomic_fasta/$1_proteins.faa -h $2 $3 $4 > ../results/seg/$1_high_complexity.faa
seg ../results/Proteomic_fasta/$1_proteins.faa -l $2 $3 $4 > ../results/seg/$1_low_complexity.faa
cat ../results/seg/$1_high_complexity.faa > ../results/seg/$1.tmp 
cat ../results/seg/$1_low_complexity.faa >> ../results/seg/$1.tmp 
echo "viral_species,ncbi_taxid,origin_beg,origin_end,codon_start,product,seg_begin,seg_end,complexity,window,locut,hicut,secuence" > ../results/seg/$1_complexity.csv
sed s'/|/,/g' ../results/seg/$1.tmp | sed -e ':a' -e 'N' -e '$!ba' -e 's/)\n/),/g' | sed -e ':a' -e 'N' -e '$!ba' -e 's/\n>/>/g' | awk -F, '{gsub("to",",",$3)}1' | awk '{gsub("from","",$3)}1' | awk '{gsub("id_","",$2)}1' | awk '{gsub("start","",$4)}1' | awk '{gsub("\\(",",",$5)}1' | awk '{gsub("-",",",$5)}1' | awk '{gsub("\\)",",",$5)}1' | sed 's/ complexity=//g' | awk '{gsub("\\(","",$6)}1' | awk '{gsub("\\/",",",$6)}1' | awk '{gsub("\\)","",$6)}1' | sed 's/ /,/g' | sed -e ':a' -e 'N' -e '$!ba' -e 's/\n//g; s/>/\n/g'  >> ../results/seg/$1_complexity.csv
sed -i '/^[[:space:]]*$/d' ../results/seg/$1_complexity.csv
rm ../results/seg/$1.tmp 
