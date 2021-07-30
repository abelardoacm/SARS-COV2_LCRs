#!/usr/local/bin/perl
#Gb_2_faa.pl
#
#Author: Abelardo Aguilar Camara
#
#Task performed:
#    This script opens a file with genbank files
#    and reads line by line to write a unique contatenated
#    fasta proteomic. Info in headers allows to locate
#    simple sequences with NCBI SEG
#
#INPUT: Concatenated full genbanks of anytaxon in ../data/Raw_database/anytaxon.gb (e.g. Geminiviridae.gb)
#
#OUTPUT: Unique fasta aminoacid proteomic file of anytaxon (*.faa)
#
#
#
##########################################
#
#Gb_2_faa.pl
#
##########################################
#
#Output location
system ("mkdir -p ../results/Proteomic_fasta/");
#
#Reach input file
$archivo = $ARGV[0]; # Reads first positional argument. Expected $ARGV[0] is anytaxon.gb
chomp ($archivo); # Removes trailing endlines
$carpeta = $archivo;
$proc_file = $archivo;
$proc_file =~ s/.gb//g;
$taxon = $proc_file;
$proc_file = "$proc_file"."_proteins.faa";
open (FILE, "../data/Raw_database/$archivo"); # Opens input
##########################################
#Read input file line by line
open (OUT, ">>../results/Proteomic_fasta/$proc_file");
while ($linea = <FILE>) { #Read by line until end of file, defines $linea variable as the current line
	chomp ($linea);
	#################################
	if (index($linea, "collection_date") != -1) {
    		$collection_date = $linea;
    		$collection_date =~ s/.*collection_date//;
		$collection_date =~ s/\,.*//;
	}
	if (index($linea, "country") != -1) {
    		$country = $linea;
    		$country =~ s/.*country//;
		$country =~ s/\,.*//;
	}
	#Define source organism from "DEFINITION" line
	if ($linea =~ /^ACCESSION\s+/) {
		$organismo = $linea;
		$organismo =~ s/ACCESSION\s+//;
		$organismo =~ s/\s+/-/g;
		$organismo =~ s/:/-/g;
		$organismo =~ s/;/-/g;
		$organismo =~ s/_/-/g;
		$organismo =~ s/\,//g;
		$organismo =~ s/\.//g;
	}
	################################
	#Recover accession identifier from "ACCESSION" line
        if ($linea =~ /^ACCESSION\s+/) {
		$accession = $linea; #Accession equals ACCESSION line
		#Remove undesired text and special characters
		$accession =~ s/ACCESSION\s+//;
		$accession =~ s/\s+/_/g;
		$accession =~ s/\,//g;
		$accession =~ s/\.//g;
		$forname = "$accession"; #Variable to write in fasta header
	}
	################################
	#Recover taxid from db_xref in source field
        if ($linea =~ /\s+\/db_xref\="taxon/) {
                $TAXID = $linea; #Taxonomic ID defined as line next to dr_xref pattern
                #Remove undesired text and special characters
                $TAXID =~ s/\s+//;
		$TAXID =~ s/\/db_xref//;
		$TAXID =~ s/taxon//;
		$TAXID =~ s/://g;
		$TAXID =~ s/=//g;
		$TAXID =~ s/"//g;
		$taxid = "taxid"; #Variable to concatenate in header
		#########################
		#Define and give proper format to individual output name
		$NombreOutput = "$organismo"."_"."$forname"."_"."$taxid"."_"."$TAXID"."_"."$faa";
		#Remove spaces and undesired text
		$NombreOutput =~ s/--/-/g;
		$NombreOutput =~ s/\(//g;
		$NombreOutput =~ s/\)//g;
		$NombreOutput =~ s/\[//g;
		$NombreOutput =~ s/\]//g;
		$NombreOutput =~ s/\{//g;
		$NombreOutput =~ s/\}//g;
		$NombreOutput =~ s/\//-/g;
		#open (OUT, ">>$NombreOutput");
        }  
        #################################
        #Save protein_id to variable $id
        if ($linea =~ /\s+\/protein_id\=/) { #when pattern "protein_id" is found...
                $id = $linea;
                $id =~ s/\s+//; # Removes spaces
		$id =~ s/\/protein_id=//; # Removes tag
		$id =~ s/"//g; #Variable containing only protein_id
        }
        #################################
        #Save taxid to variable $ID
        if ($linea =~ /\s+\/db_xref\=/) { #when pattern "db_xref" is found...
                $ID = $linea;
                $ID =~ s/\s+//; # Removes spaces
                $ID =~ s/\/db_xref=//; # Removes tag
                $ID =~ s/"//g; # Variable containing only taxonomic id
        }
        #################################
        #Save GI to variable $GI
        if ($linea =~ /\s+\/db_xref\="GI/) { #when pattern "db_xref="GI" is found...
                $GI = $linea;
                #Remove spaces and undesired text
                $GI =~ s/\s+//;
                $GI =~ s/\/db_xref\="GI//;
                $GI =~ s/"//g;
		$G = "GI";
		$conc = $G.$GI; #Variable to concatenate in protein fasta header
	}
	#################################
        #Save locus_tag to variable $locus
        if ($linea =~ /\s+\/locus_tag\=/) { #when pattern "db_xref="GI" is found...
                $locus = $linea;
                #Remove spaces and undesired text
                $locus =~ s/\s+//;
                $locus =~ s/\/locus_tag\=//;
                $locus =~ s/"//g;
		$l = "locus-";
                $con = $l.$locus; #Variable to concatenate in protein fasta header
        }
        #################################
        #Save protein name to variable $proteina
	if ($linea =~ /\s+\/product\=/) { #when pattern "db_xref="GI" is found...
		$proteina = $linea;
		#Remove spaces and undesired text
		$proteina =~ s/\s+\/product\=//;
		$proteina =~ s/"//g;
	}
	if ($linea =~ /\s+\/codon_start\=/) { #when pattern "db_xref="GI" is found...
		$codon_start = $linea;
		#Remove spaces and undesired text
		$codon_start =~ s/\s+\/codon_start\=//;
		$codon_start =~ s/"//g;
		$codon_start =~ s/\s+//;
	}
	if ($linea =~ /\s+\/db_xref\="taxon/) {
                $TAXID = $linea; #Taxonomic ID defined as line next to dr_xref pattern
                #Remove undesired text and special characters
                $TAXID =~ s/\s+//;
		$TAXID =~ s/\/db_xref//;
		$TAXID =~ s/taxon//;
		$TAXID =~ s/://g;
		$TAXID =~ s/=//g;
		$TAXID =~ s/"//g;
		$taxid = "taxid"; #Variable to concatenate in header
	}
	if ($linea =~ /\s+CDS\s+/ || $linea =~ /\s+gene\s+/ ) { #when pattern "CDS" is found...
                $POSICION = $linea;
                $POSICION =~ s/\s+//; # Removes spaces
		$POSICION =~ s/CDS//; # Removes tag
		$POSICION =~ s/"//g;
		$POSICION =~ s/\s//g;
		$POSICION =~ s/ //g;
		$POSICION =~ s/\,/_/g;
        }
        if (index($POSICION, $isjoin) != -1) {
		$chunk1 = $POSICION;
		$chunk2 = $POSICION;
		$chunk1 =~ s/.*join\(//;
		$chunk1 =~ s/_.*//;
		$chunk2 =~ s/.*_//;
		$chunk2 =~ s/\).*//;
		$chunk1 =~ s/[^0-9.]+//g;
		$chunk2 =~ s/[^0-9.]+//g;
		$chunk1begin = "BEGIN"."$chunk1";
		$chunk1end = "$chunk1"."END";
		$chunk2begin = "BEGIN"."$chunk2";
		$chunk2end = "$chunk2"."END";
		$chunk1begin =~ s/.*BEGIN//;
		$chunk1begin =~ s/\.\..*//;
		$chunk1end =~ s/.*\.\.//;
		$chunk1end =~ s/END.*//;
		$chunk2begin =~ s/.*BEGIN//;
		$chunk2begin =~ s/\.\..*//;
		$chunk2end =~ s/.*\.\.//;
		$chunk2end =~ s/END.*//;
	} else {
		$POSICION =~ s/[^0-9.]+//g; #remove everything but numbers and dot
		$chunk1begin = "BEGIN"."$POSICION";
		$chunk1end = "$POSICION"."END";
		$chunk1begin =~ s/.*BEGIN//;
		$chunk1begin =~ s/\.\..*//;
		$chunk2end =~ s/.*\.\.//;
		$chunk2end =~ s/END.*//;
	}
	$proteina =~ s/ /_/;
	$proteina =~ s/\s+/_/;
	$proteina =~ s/[^a-zA-Z0-9,]/_/g;
	$proteina =~ s/,//g;
	$organismo =~ s/,//g;
	$TAXID =~ s/,//g;
	$chunk1begin =~ s/,//g;
	$chunk2end =~ s/,//g;
	$codon_start =~ s/,//g;
	$virusname = "$organismo$country$collection_date";
	$virusname =~ s/[^a-zA-Z0-9]/_/g;
	#################################
	#Write fasta header of protein
	if ($linea =~ /\s+\/translation\=/) { #If $linea contains "translation" tag, recover lines corresponding to the aminoacid sequence of a protein
		$numero = 1; #Line meets condition of being part of aminacids sequence
		#print OUT "> $organismo id_$TAXID | $id | $proteina\n"; # Write fasta header of protein
		print OUT ">$virusname"."|id_$TAXID|from$chunk1begin"."to$chunk2end|start$codon_start|$proteina\n"; # Write fasta header of protein
		#print OUT ">$TAXID|$chunk1begin"."to$chunk2end|$codon_start|$proteina\n";
		#EXAMPLE: >YP_007974221.1|GeneID-15486949|GI-498907898|locus-L677_gp1|protein_AV2-protein_|
	}
	#################################
	#Line is not part of proteomic fasta when...
	#Lines to ignore (1st set)
	if ($linea =~ /\s+gene\s+/ || $linea =~ /^ORIGIN/ || $linea =~ /\s+CDS\s+/ || $linea =~ /\s+polyA_site\s+/ || $linea =~ /\s+repeat_region\s+/ || $linea =~ /\s+polyA_signal\s+/ || $linea =~ /\s+rep_origen\s+/ || $linea =~ /\s+promoter\s+/ || $linea =~ /\s+sig_peptide\s+/ || $linea =~ /\s+misc_feature\s+/ || $linea =~ /\s+5'UTR\s+/) {
		$numero = 0; # Line does not correspond to aminoacids sequence, condition variable set to 0.
	}
	#Lines to ignore (2nd set)
        if ($linea =~ /\s+stem_loop\s+/ || $linea =~ /\s+3'UTR\s+/ || $linea =~ /\s+CDS\s+complement/ || $linea =~ /\s+regulatory\s+/ || $linea =~ /\/product\=/ || $linea =~ /\s+intron\s+/ || $linea =~ /\s+unsure\s+/ ||  $linea =~ /\/gene\=/ || $linea =~ /\s+tRNA\s+/ || $linea =~ /\s+exon\s+/ || $linea =~ /\s+ncRNA\s+/ || $linea =~ /\s+variation\s+/ || $linea =~ /\TATA_signal\s+/ || $linea =~ /\/locus_tag\=/ || $linea =~ /\/note\=/) {
                $numero = 0;
	}
	#Lines to ignore (3rd set)
 	if ($linea =~ /\s+STS\s+/ || $linea =~ /\/standard_name\=/ || $linea =~ /mat_peptide/ || $linea =~ /\/db_xref\=/ || $linea =~ /\s+rep_origin\s+/ || $linea =~ /\s+misc_recomb\s+/ || $linea =~ /\s+protein_bind\s+/ ) {  
		 $numero = 0;
 	}
 	#################################
 	#Line is part of proteomic fasta file when...
	if ($numero == 1) { #Variable $linea contains "translation" tag
		#Remove undesired spaces and text
		$linea =~ s/\s+//;
		$linea =~ s/\/translation\=//;
		$linea =~ s/"//g;
		print OUT "$linea\n"; #Print line in proteomic fasta file
	} #END OF INDIVIDUAL ENTRY
} #END OF CONCATENATED GENBANK
close (FILE); #Close output at input end
