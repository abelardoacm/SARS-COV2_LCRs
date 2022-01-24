# SARS-COV2_LCRs
SARS-COV2_LCRs is a short pipeline to plot low complexity regions (LCRs) prevalence from SARS-CoV2 genbank files, using [NCBI seg](https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&ved=2ahUKEwjXhqnZxPLwAhXIm-AKHXmoAp4QFjAAegQIBBAD&url=ftp%3A%2F%2Fftp.ncbi.nlm.nih.gov%2Fpub%2Fseg%2Fseg%2F&usg=AOvVaw2s1FT-lfX5HmgPegjJk2tB) output. It's main objective is to automatize seg analysis and easily locate LCRs of interest in SARS-CoV2 proteomes.

This scripts were used for LCRs detection in the article ["Two short low complexity regions (LCRs) are hallmark sequences of the Delta SARS-CoV-2 variant spike protein"](https://www.nature.com/articles/s41598-022-04976-8).


>### Dependencies
These are SARS-COV2_LCRs dependencies, links to their installing instructions, commands for installation, and used versions:
 | Software name (version used) 	| Terminal 	| Installation |
|-	|-	|-	|
| [tidyverse](https://www.tidyverse.org/) (1.3.0) 	| R (4.0.3) 	| `install.packages("tidyverse")` 	|
| [viridis](https://cran.r-project.org/web/packages/viridis/viridis.pdf)(0.6.1) | R (4.0.3)	| `install.packages("viridis")` |
| [scales](https://cran.r-project.org/web/packages/scales/index.html)(1.1.1) | R (4.0.3)	| `install.packages("scales")` |
| [grid](https://cran.r-project.org/web/packages/grid/index.html)(4.1.0) | R (4.0.3)	| `install.packages("grid")` |
| [NCBI seg](https://www.biostars.org/p/424116/) | bash 	| download via ftp, compile and export permanently to *$PATH* |


It also assumes properly functioning [**perl**](https://www.perl.org/) (tested with v5.30.0), and [**seg**](https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&ved=2ahUKEwjXhqnZxPLwAhXIm-AKHXmoAp4QFjAAegQIBBAD&url=ftp%3A%2F%2Fftp.ncbi.nlm.nih.gov%2Fpub%2Fseg%2Fseg%2F&usg=AOvVaw2s1FT-lfX5HmgPegjJk2tB) working under any location. For instance, this is the expected output of entering `seg` at commnand line:

```
Usage: seg <file> <window> <locut> <hicut> <options>
         <file>   - filename containing fasta-formatted sequence(s) 
         <window> - OPTIONAL window size (default 12) 
         <locut>  - OPTIONAL low (trigger) complexity (default 2.2) 
         <hicut>  - OPTIONAL high (extension) complexity (default locut + 0.3) 
	 <options> 
            -x  each input sequence is represented by a single output 
                sequence with low-complexity regions replaced by 
                strings of 'x' characters 
            -c <chars> number of sequence characters/line (default 60)
            -m <size> minimum length for a high-complexity segment 
                (default 0).  Shorter segments are merged with adjacent 
                low-complexity segments 
            -l  show only low-complexity segments (fasta format) 
            -h  show only high-complexity segments (fasta format) 
            -a  show all segments (fasta format) 
            -n  do not add complexity information to the header line 
            -o  show overlapping low-complexity segments (default merge) 
            -t <maxtrim> maximum trimming of raw segment (default 100) 
            -p  prettyprint each segmented sequence (tree format) 
            -q  prettyprint each segmented sequence (block format)
```
>### Cloning this repo
To clone this repo via command line git, enter the following commands:
```
git clone https://github.com/abelardoacm/SARS-COV2_LCRs.git
cd SARS-COV2_LCRs/bin/
chmod +x *
```
if all requirements are met you should be ready to go

>### Repo tree

``` 
.
├── bin <-  Where you have to be to invoke scripts
│
├── data
│   └── Raw_database <- Place here genomic genbank files 
│
└── results
```
>### Using SARS-COV2_LCRs

This repository contains all the scripts used in the research article titled: "Two Short Low Complexity Regions (LCRs) are Hallmark Sequences of the Delta SARS-CoV- 2 Variant Spike Protein", [currently published in Scientific Reports (Becerra et al., 2022)](https://t.co/oYl7vinB49?amp=1). To access the complete database, please check the supplementary files for the list of NCBI accession IDs of analyzed files for each variant. A sample dataset containing Beta variant genomes is shared [here](/data/Raw_database/).

To reproduce the analysis, from [/bin/](/bin/) type the following command:
```
./SARS-COV2_LCRs.sh
```
you will be prompted for the **seg** parameters one at a time (we used, window = 12, locut = 1.9 and hicut = 2.1).

Note that the previous command performs the full analysis for each **.gb** file contained in [/data/Raw_database/](/data/Raw_database/). It begins invoking [BD_seg.sh](bin/BD_seg.sh), which reads each **.gb** file to build a fasta file with genomic coordinates on headers. Then NCBI seg is used via [Just_a_seg_envelope.sh](bin/Just_a_seg_envelope.sh). Output is saved in a csv file suitable for exploration with R. The pipeline continues with the script [Analyze_LCRs.R](bin/Analyze_LCRs.R) that computes counts of the SARS-CoV2 LCRs of interest and exports **.tiff** barplots to [results](/results/). Finally, [Plot_proteome_complexity.R](bin/Plot_proteome_complexity.R) graphs the estimated average complexity for each position in analyzed proteomes per file.
