# SARS-COV2_LCRs
SARS-COV2_LCRs is a short pipeline to plot low complexity regions (LCRs) prevalence from SARS-CoV2 genbank files, using [NCBI seg](https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&ved=2ahUKEwjXhqnZxPLwAhXIm-AKHXmoAp4QFjAAegQIBBAD&url=ftp%3A%2F%2Fftp.ncbi.nlm.nih.gov%2Fpub%2Fseg%2Fseg%2F&usg=AOvVaw2s1FT-lfX5HmgPegjJk2tB) output. It's main objective is to automatize seg analysis and easily locate LCRs of interest in SARS-CoV2 proteomes.

>### Dependencies
These are SARS-COV2_LCRs dependencies, links to their installing instructions, commands for installation, and used versions:
 | Software name (version used) 	| Terminal 	| Installation *debian based distros 	|
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
    ├── GenFeatures_locations <- Contains csv files with genomic features extracted from genbank file
    ├── Proteomic_fasta <- Fasta aminoacid files to serve as seg input
    ├── seg <- NCBI seg analysis output for high (.faa), low (.faa) and both (.csv) type of sequences
```
>### Using SARS-COV2_LCRs
