#!/usr/local/bin/perl
########################
#######################
######################
####################
# CREO QUE DEBERIA AÑADIR LA COLUMNA CON TRANSLATION CUANDO LO ENCUENTRE, LUEGO GENERAR UN IDENTIFICADOR UNICO PARA LOS CDS (QUE SON LOS QUE TIENEN TRANSLATION)
# CUANDO HAGA EL INPUT PARA SECUENCIAS SIMPLES HARE QUE EL HEADER DEL FASTA INPUT CONTENGA ESE IDENTIFICADOR, Y LA POSICION DE ESE CDS EN SU GENOMA
#
# LA FINALIDAD ES QUE PUEDA SABER LA POSICION DE LAS SECUECNIAS SIMPLES Y QUE ESTAS PUEDAN CONVERTIRSE EN UNA FILA MAS DE LA BASE DE DATOS QUE GENERA ESTE PROGRAMA
#
#
#PROBABLEMENTE TENGA QUE DIVIDIR LAS FILAS QUE DIGAN JOIN EN TANTAS FILAS COMO SEA NECESARIO, ESTO PARA PODER HACER UNA COLUMNA QUE SEA "INICIO" y otra "FIN"
####################
####################
#WEAKPOINTS
####################
#If there are more segments of a feature, for example, multiple joints in $position variable, then program would only recover the first two
system ("mkdir -p ../results/GenFeatures_locations");
$archivo = $ARGV[0]; # Reads first positional argument. Expected $ARGV[0] is anytaxon.gb
chomp ($archivo); # Removes trailing endlines
$carpeta = $archivo;
$proc_file = $archivo;
$proc_file =~ s/.gb//g;
$taxon = $proc_file;
$proc_file = "$proc_file"."_peptide_locations.csv";
open (FILE, "../data/Raw_database/$archivo"); # Opens input
##########################################
#Read input file line by line
$numero = 0; #Variable to indicate favorable (1) or absent (0) condition
open (OUT, ">>../results/GenFeatures_locations/$proc_file");
while ($linea = <FILE>) { #Read by line until end of file, defines variable $linea as the current line
	chomp ($linea);
	#################################
	#Following lines are found once by gb file
	#################################
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
	if ($linea =~ /^FEATURES\s+/) {
		$numero = 1;
	}
	if ($linea =~ /^ORIGIN\s+/) {
		$numero = 0;
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
		#########################
        }  
        #################################
        #LAS QUE ME INTERESAN
        #################################
        #5'UNTRANSLATED REGION
        if ($linea =~ /\s+5'UTR\s+/ ) { #when pattern "5'UTR" is found...
                $FIVE_UTR = $linea;
                $FIVE_UTR =~ s/\s+//; # Removes spaces
		$FIVE_UTR =~ s/5'UTR//; # Removes tag
		$FIVE_UTR =~ s/"//g;
		$FIVE_UTR =~ s/\s//g;
		$FIVE_UTR =~ s/ //g;
		$FIVE_UTR =~ s/\,/_/g;
		print OUT "\nORGANISMO$organismo,TAXID$TAXID,CARACTERISTICA5_UTR,POSICION$FIVE_UTR";
        }
        #3'UNTRANSLATED REGION
        if ($linea =~ /\s+3'UTR\s+/ ) { #when pattern "3'UTR" is found...
                $THREE_UTR = $linea;
                $THREE_UTR =~ s/\s+//; # Removes spaces
		$THREE_UTR =~ s/3'UTR//; # Removes tag
		$THREE_UTR =~ s/"//g;
		$THREE_UTR =~ s/\s//g;
		$THREE_UTR =~ s/ //g;
		$THREE_UTR =~ s/\,/_/g;
		print OUT "\nORGANISMO$organismo,TAXID$TAXID,CARACTERISTICA3_UTR,POSICION$THREE_UTR";
        }
        #gene
        if ($linea =~ /\s+gene\s+/ ) { #when pattern "gene" is found...
                $gene = $linea;
                $gene =~ s/\s+//; # Removes spaces
		$gene =~ s/gene//; # Removes tag
		$gene =~ s/"//g;
		$gene =~ s/\s//g;
		$gene =~ s/ //g;
		$gene =~ s/\,/_/g;
		print OUT "\nORGANISMO$organismo,TAXID$TAXID,CARACTERISTICAgene,POSICION$gene";
        }
        #CDS
        if ($linea =~ /\s+CDS\s+/ ) { #when pattern "CDS" is found...
                $CDS = $linea;
                $CDS =~ s/\s+//; # Removes spaces
		$CDS =~ s/CDS//; # Removes tag
		$CDS =~ s/"//g;
		$CDS =~ s/\s//g;
		$CDS =~ s/ //g;
		$CDS =~ s/\,/_/g;
		print OUT "\nORGANISMO$organismo,TAXID$TAXID,CARACTERISTICACDS,POSICION$CDS";
        }
        #mat_peptide
        if ($linea =~ /\s+mat_peptide\s+/ ) { #when pattern "mat_peptide" is found...
                $mat_peptide = $linea;
                $mat_peptide =~ s/\s+//; # Removes spaces
		$mat_peptide =~ s/mat_peptide//; # Removes tag
		$mat_peptide =~ s/"//g;
		$mat_peptide =~ s/\s//g;
		$mat_peptide =~ s/ //g;
		$mat_peptide =~ s/\,/_/g;
		print OUT "\nORGANISMO$organismo,TAXID$TAXID,CARACTERISTICAmat_peptide,POSICION$mat_peptide";
        }
        #misc_feature
        if ($linea =~ /\s+misc_feature\s+/ ) { #when pattern "misc_feature" is found...
                $misc_feature = $linea;
                $misc_feature =~ s/\s+//; # Removes spaces
		$misc_feature =~ s/misc_feature//; # Removes tag
		$misc_feature =~ s/"//g;
		$misc_feature =~ s/\s//g;
		$misc_feature =~ s/ //g;
		$misc_feature =~ s/\,/_/g;
		print OUT "\nORGANISMO$organismo,TAXID$TAXID,CARACTERISTICAmisc_feature,POSICION$misc_feature";
        }
        #misc_RNA
        if ($linea =~ /\s+misc_RNA\s+/ ) { #when pattern "misc_RNA" is found...
                $misc_RNA = $linea;
                $misc_RNA =~ s/\s+//; # Removes spaces
		$misc_RNA =~ s/misc_RNA//; # Removes tag
		$misc_RNA =~ s/"//g;
		$misc_RNA =~ s/\s//g;
		$misc_RNA =~ s/ //g;
		$misc_RNA =~ s/\,/_/g;
		print OUT "\nORGANISMO$organismo,TAXID$TAXID,CARACTERISTICAmisc_RNA,POSICION$misc_RNA";
        }
        #misc_structure
        if ($linea =~ /\s+misc_structure\s+/ ) { #when pattern "misc_structure" is found...
                $misc_structure = $linea;
                $misc_structure =~ s/\s+//; # Removes spaces
		$misc_structure =~ s/misc_structure//; # Removes tag
		$misc_structure =~ s/"//g;
		$misc_structure =~ s/\s//g;
		$misc_structure =~ s/ //g;
		$misc_structure =~ s/\,/_/g;
		print OUT "\nORGANISMO$organismo,TAXID$TAXID,CARACTERISTICAmisc_structure,POSICION$misc_structure";
        }
        #mRNA
        if ($linea =~ /\s+mRNA\s+/ ) { #when pattern "mRNA" is found...
                $mRNA = $linea;
                $mRNA =~ s/\s+//; # Removes spaces
		$mRNA =~ s/mRNA//; # Removes tag
		$mRNA =~ s/"//g;
		$mRNA =~ s/\s//g;
		$mRNA =~ s/ //g;
		$mRNA =~ s/\,/_/g;
		print OUT "\nORGANISMO$organismo,TAXID$TAXID,CARACTERISTICAmRNA,POSICION$mRNA";
        }
        #regulatory
        if ($linea =~ /\s+regulatory\s+/ ) { #when pattern "regulatory" is found...
                $regulatory = $linea;
                $regulatory =~ s/\s+//; # Removes spaces
		$regulatory =~ s/regulatory//; # Removes tag
		$regulatory =~ s/"//g;
		$regulatory =~ s/\s//g;
		$regulatory =~ s/ //g;
		$regulatory =~ s/\,/_/g;
		print OUT "\nORGANISMO$organismo,TAXID$TAXID,CARACTERISTICAregulatory,POSICION$regulatory";
        }
        #sig_peptide
        if ($linea =~ /\s+sig_peptide\s+/ ) { #when pattern "sig_peptide" is found...
                $sig_peptide = $linea;
                $sig_peptide =~ s/\s+//; # Removes spaces
		$sig_peptide =~ s/sig_peptide//; # Removes tag
		$sig_peptide =~ s/"//g;
		$sig_peptide =~ s/\s//g;
		$sig_peptide =~ s/ //g;
		$sig_peptide =~ s/\,/_/g;
		print OUT "\nORGANISMO$organismo,TAXID$TAXID,CARACTERISTICAsig_peptide,POSICION$sig_peptide";
        }
        #source
        if ($linea =~ /\s+source\s+/ ) { #when pattern "source" is found...
                $source = $linea;
                $source =~ s/\s+//; # Removes spaces
		$source =~ s/source//; # Removes tag
		$source =~ s/"//g;
		$source =~ s/\s//g;
		$source =~ s/ //g;
		$source =~ s/\,/_/g;
		print OUT "\nORGANISMO$organismo,TAXID$TAXID,CARACTERISTICAsource,POSICION$source";
        }
        #stem_loop
        if ($linea =~ /\s+stem_loop\s+/ ) { #when pattern "stem_loop" is found...
                $stem_loop = $linea;
                $stem_loop =~ s/\s+//; # Removes spaces
		$stem_loop =~ s/stem_loop//; # Removes tag
		$stem_loop =~ s/"//g;
		$stem_loop =~ s/\s//g;
		$stem_loop =~ s/ //g;
		$stem_loop =~ s/\,/_/g;
		print OUT "\nORGANISMO$organismo,TAXID$TAXID,CARACTERISTICAstem_loop,POSICION$stem_loop";
        }
       if ($linea !~ /\s+\/translation\=/) {
		if ($numero == 1) {
			if ($linea =~ /\s+\//) {
				$to_cat = $linea;
				$to_cat =~ s/\s+//;
				$to_cat =~ s/\,/_/;
				print OUT ",$to_cat";
			}			
		}
	}
} #END OF CONCATENATED GENBANK
close (FILE); #Close output at input end
system ("mv ../results/GenFeatures_locations/$taxon\_peptide_locations.csv ../results/GenFeatures_locations/$taxon.tmp");
system ("awk 'NF' ../results/GenFeatures_locations/$taxon.tmp > ../results/GenFeatures_locations/$taxon\_peptide_locations.csv");
$archivo = $ARGV[0]; # Reads first positional argument. Expected $ARGV[0] is anytaxon.gb
chomp ($archivo); # Removes trailing endlines
$carpeta = $archivo;
$proc_file = $archivo;
$proc_file =~ s/.gb//g;
$csvout = "$proc_file"."_features_locations.csv";
$proc_file = "$proc_file"."_peptide_locations.csv";
open (FILE, "../results/GenFeatures_locations/$proc_file"); # Opens input
##########################################
#Read input file line by line
$numero = 0; #Variable to indicate favorable (1) or absent (0) condition
open (OUT, ">>../results/GenFeatures_locations/$csvout");
print OUT  "Viral_species,NCBI_taxid,Region_feature_class,Beginning,End,locus_tag,db_xref,gene,product,protein_id,codon_start,note,mol_type,host,country,collection_date,strain,ribosomal_slippage,isolate,ribosomal_slippage,inference,experiment,function,isolation_source,old_locus_tag,regulatory_class,gene_synonym,collected_by,lab_host,genotype,pseudo,lat_lon,cell_line,acronym,transl_except,serotype,segment,exception,citation,altitude,subgroup,standard_name,identified_by,group,culture_collection,clone\n";
while ($linea = <FILE>) { #Read by line until end of file, defines variable $linea as the current line
	chomp ($linea);
	#################
	#Each posible feature is saved in a variable
	$POSICION = "POSICION";
	if (index($linea, $POSICION) != -1) {
    		$POSICION = $linea;
    		$POSICION =~ s/.*POSICION//;
		$POSICION =~ s/\,.*//;
		#$POSICION =~ s/\.\..*/,/;
	} else {
    		$POSICION = "";
	}	
	$acronym = "acronym";
	if (index($linea, $acronym) != -1) {
    		$acronym = $linea;
    		$acronym =~ s/.*acronym//;
		$acronym =~ s/\,.*//;
	} else {
    		$acronym = "";
	}
	$altitude = "altitude";
	if (index($linea, $altitude) != -1) {
    		$altitude = $linea;
    		$altitude =~ s/.*altitude//;
		$altitude =~ s/\,.*//;
	} else {
    		$altitude = "";
	}
	$cell_line = "cell_line";
	if (index($linea, $cell_line) != -1) {
    		$cell_line = $linea;
    		$cell_line =~ s/.*cell_line//;
		$cell_line =~ s/\,.*//;
	} else {
    		$cell_line = "";
	}
	$citation = "citation";
	if (index($linea, $citation) != -1) {
    		$citation = $linea;
    		$citation =~ s/.*citation//;
		$citation =~ s/\,.*//;
	} else {
    		$citation = "";
	}
	$clone = "clone";
	if (index($linea, $clone) != -1) {
    		$clone = $linea;
    		$clone =~ s/.*clone//;
		$clone =~ s/\,.*//;
	} else {
    		$clone = "";
	}
	$codon_start = "codon_start";
	if (index($linea, $codon_start) != -1) {
    		$codon_start = $linea;
    		$codon_start =~ s/.*codon_start//;
		$codon_start =~ s/\,.*//;
	} else {
    		$codon_start = "";
	}
	$collected_by = "collected_by";
	if (index($linea, $collected_by) != -1) {
    		$collected_by = $linea;
    		$collected_by =~ s/.*collected_by//;
		$collected_by =~ s/\,.*//;
	} else {
    		$collected_by = "";
	}
	$collection_date = "collection_date";
	if (index($linea, $collection_date) != -1) {
    		$collection_date = $linea;
    		$collection_date =~ s/.*collection_date//;
		$collection_date =~ s/\,.*//;
	} else {
    		$collection_date = "";
	}
	$country = "country";
	if (index($linea, $country) != -1) {
    		$country = $linea;
    		$country =~ s/.*country//;
		$country =~ s/\,.*//;
	} else {
    		$country = "";
	}
	$culture_collection = "culture_collection";
	if (index($linea, $culture_collection) != -1) {
    		$culture_collection = $linea;
    		$culture_collection =~ s/.*culture_collection//;
		$culture_collection =~ s/\,.*//;
	} else {
    		$culture_collection = "";
	}
	$db_xref = "db_xref";
	if (index($linea, $db_xref) != -1) {
    		$db_xref = $linea;
    		$db_xref =~ s/.*db_xref//;
		$db_xref =~ s/\,.*//;
	} else {
    		$db_xref = "";
	}
	$exception = "exception";
	if (index($linea, $exception) != -1) {
    		$exception = $linea;
    		$exception =~ s/.*exception//;
		$exception =~ s/\,.*//;
	} else {
    		$exception = "";
	}
	$experiment = "experiment";
	if (index($linea, $experiment) != -1) {
    		$experiment = $linea;
    		$experiment =~ s/.*experiment//;
		$experiment =~ s/\,.*//;
	} else {
    		$experiment = "";
	}
	$function = "function";
	if (index($linea, $function) != -1) {
    		$function = $linea;
    		$function =~ s/.*function//;
		$function =~ s/\,.*//;
	} else {
    		$function = "";
	}
	$gene = "gene";
	if (index($linea, $gene) != -1) {
    		$gene = $linea;
    		$gene =~ s/.*gene//;
		$gene =~ s/\,.*//;
	} else {
    		$gene = "";
	}
	$gene_synonym = "gene_synonym";
	if (index($linea, $gene_synonym) != -1) {
    		$gene_synonym = $linea;
    		$gene_synonym =~ s/.*gene_synonym//;
		$gene_synonym =~ s/\,.*//;
	} else {
    		$gene_synonym = "";
	}
	$genotype = "genotype";
	if (index($linea, $genotype) != -1) {
    		$genotype = $linea;
    		$genotype =~ s/.*genotype//;
		$genotype =~ s/\,.*//;
	} else {
    		$genotype = "";
    		$inference = "";
	}
	$isolate = "isolate";
	if (index($linea, $isolate) != -1) {
    		$isolate = $linea;
    		$isolate =~ s/.*isolate//;
		$isolate =~ s/\,.*//;
	} else {
    		$isolate = "";
	}
	$isolation_source = "isolation_source";
	if (index($linea, $isolation_source) != -1) {
    		$isolation_source = $linea;
    		$isolation_source =~ s/.*isolation_source//;
		$isolation_source =~ s/\,.*//;
	} else {
    		$isolation_source = "";
	}
	$lab_host = "lab_host";
	if (index($linea, $lab_host) != -1) {
    		$lab_host = $linea;
    		$lab_host =~ s/.*lab_host//;
		$lab_host =~ s/\,.*//;
	} else {
    		$lab_host = "";
	}
	$lat_lon = "lat_lon";
	if (index($linea, $lat_lon) != -1) {
    		$lat_lon = $linea;
    		$lat_lon =~ s/.*lat_lon//;
		$lat_lon =~ s/\,.*//;
	} else {
    		$lat_lon = "";
	}
	$locus_tag = "locus_tag";
	if (index($linea, $locus_tag) != -1) {
    		$locus_tag = $linea;
    		$locus_tag =~ s/.*locus_tag//;
		$locus_tag =~ s/\,.*//;
	} else {
    		$locus_tag = "";
	}
	$mol_type = "mol_type";
	if (index($linea, $mol_type) != -1) {
    		$mol_type = $linea;
    		$mol_type =~ s/.*mol_type//;
		$mol_type =~ s/\,.*//;
	} else {
    		$mol_type = "";
	}
	$note = "note";
	if (index($linea, $note) != -1) {
    		$note = $linea;
    		$note =~ s/.*note//;
		$note =~ s/\,.*//;
	} else {
    		$note = "";
	}
	$old_locus_tag = "old_locus_tag";
	if (index($linea, $old_locus_tag) != -1) {
    		$old_locus_tag = $linea;
    		$old_locus_tag =~ s/.*old_locus_tag//;
		$old_locus_tag =~ s/\,.*//;
	} else {
    		$old_locus_tag = "";
	}
	$organism = "organism";
	if (index($linea, $organism) != -1) {
    		$organism = $linea;
    		$organism =~ s/.*organism//;
		$organism =~ s/\,.*//;
	} else {
    		$organism = "";
	}
	$product = "product";
	if (index($linea, $product) != -1) {
    		$product = $linea;
    		$product =~ s/.*product//;
		$product =~ s/\,.*//;
	} else {
    		$product = "";
	}
	$protein_id = "protein_id";
	if (index($linea, $protein_id) != -1) {
    		$protein_id = $linea;
    		$protein_id =~ s/.*protein_id//;
		$protein_id =~ s/\,.*//;
	} else {
    		$protein_id = "";
	}
	$regulatory_class = "regulatory_class";
	if (index($linea, $regulatory_class) != -1) {
    		$regulatory_class = $linea;
    		$regulatory_class =~ s/.*regulatory_class//;
		$regulatory_class =~ s/\,.*//;
	} else {
    		$regulatory_class = "";
	}
	$segment = "segment";
	if (index($linea, $segment) != -1) {
    		$segment = $linea;
    		$segment =~ s/.*segment//;
		$segment =~ s/\,.*//;
	} else {
    		$segment = "";
	}
	$serotype = "serotype";
	if (index($linea, $serotype) != -1) {
    		$serotype = $linea;
    		$serotype =~ s/.*serotype//;
		$serotype =~ s/\,.*//;
	} else {
    		$serotype = "";
	}
	$standard_name = "standard_name";
	if (index($linea, $standard_name) != -1) {
    		$standard_name = $linea;
    		$standard_name =~ s/.*standard_name//;
		$standard_name =~ s/\,.*//;
	} else {
    		$standard_name = "";
	}
	$strain = "strain";
	if (index($linea, $strain) != -1) {
    		$strain = $linea;
    		$strain =~ s/.*strain//;
		$strain =~ s/\,.*//;
	} else {
    		$strain = "";
	}
	$subgroup = "subgroup";
	if (index($linea, $subgroup) != -1) {
    		$subgroup = $linea;
    		$subgroup =~ s/.*subgroup//;
		$subgroup =~ s/\,.*//;
	} else {
    		$subgroup = "";
	}
	$transl_except = "transl_except";
	if (index($linea, $transl_except) != -1) {
    		$transl_except = $linea;
    		$transl_except =~ s/.*transl_except//;
		$transl_except =~ s/\,.*//;
	} else {
    		$transl_except = "";
	}
	$ORGANISMO = "ORGANISMO";
	if (index($linea, $ORGANISMO) != -1) {
    		$ORGANISMO = $linea;
    		$ORGANISMO =~ s/.*ORGANISMO//;
		$ORGANISMO =~ s/\,.*//;
	} else {
    		$ORGANISMO = "";
	}
	$TAXID = "TAXID";
	if (index($linea, $TAXID) != -1) {
    		$TAXID = $linea;
    		$TAXID =~ s/.*TAXID//;
		$TAXID =~ s/\,.*//;
	} else {
    		$TAXID = "";
	}
	$CARACTERISTICA = "CARACTERISTICA";
	if (index($linea, $CARACTERISTICA) != -1) {
    		$CARACTERISTICA = $linea;
    		$CARACTERISTICA =~ s/.*CARACTERISTICA//;
		$CARACTERISTICA =~ s/\,.*//;
	} else {
    		$CARACTERISTICA = "";
	}
	$group = "group";
	if (index($linea, $group) != -1) {
    		$group = $linea;
    		$group =~ s/.*group//;
		$group =~ s/\,.*//;
	} else {
    		$group = "";
	}
	$host = "host";
	if (index($linea, $host) != -1) {
    		$host = $linea;
    		$host =~ s/.*host//;
		$host =~ s/\,.*//;
	} else {
    		$host = "";
	}
	$identified_by = "identified_by";
	if (index($linea, $identified_by) != -1) {
    		$identified_by = $linea;
    		$identified_by =~ s/.*identified_by//;
		$identified_by =~ s/\,.*//;
	} else {
    		$identified_by = "";
	}
	$inference = "inference";
	if (index($linea, $inference) != -1) {
    		$inference = $linea;
    		$inference =~ s/.*inference//;
		$inference =~ s/\,.*//;
	} else {
    		$inference = "";
	}
	$isolate = "isolate";
	if (index($linea, $isolate) != -1) {
    		$isolate = $linea;
    		$isolate =~ s/.*isolate//;
		$isolate =~ s/\,.*//;
	} else {
    		$isolate = "";
	}
	$isolation_source = "isolation_source";
	if (index($linea, $isolation_source) != -1) {
    		$isolation_source = $linea;
    		$isolation_source =~ s/.*isolation_source//;
		$isolation_source =~ s/\,.*//;
	} else {
    		$isolation_source = "";
	}
	$lab_host = "lab_host";
	if (index($linea, $lab_host) != -1) {
    		$lab_host = $linea;
    		$lab_host =~ s/.*lab_host//;
		$lab_host =~ s/\,.*//;
	} else {
    		$lab_host = "";
	}
	$lat_lon = "lat_lon";
	if (index($linea, $lat_lon) != -1) {
    		$lat_lon = $linea;
    		$lat_lon =~ s/.*lat_lon//;
		$lat_lon =~ s/\,.*//;
	} else {
    		$lat_lon = "";
	}
	$locus_tag = "locus_tag";
	if (index($linea, $locus_tag) != -1) {
    		$locus_tag = $linea;
    		$locus_tag =~ s/.*locus_tag//;
		$locus_tag =~ s/\,.*//;
	} else {
    		$locus_tag = "";
	}
	$mol_type = "mol_type";
	if (index($linea, $mol_type) != -1) {
    		$mol_type = $linea;
    		$mol_type =~ s/.*mol_type//;
		$mol_type =~ s/\,.*//;
	} else {
    		$mol_type = "";
	}
	$note = "note";
	if (index($linea, $note) != -1) {
    		$note = $linea;
    		$note =~ s/.*note//;
		$note =~ s/\,.*//;
	} else {
    		$note = "";
	}
	$old_locus_tag = "old_locus_tag";
	if (index($linea, $old_locus_tag) != -1) {
    		$old_locus_tag = $linea;
    		$old_locus_tag =~ s/.*old_locus_tag//;
		$old_locus_tag =~ s/\,.*//;
	} else {
    		$old_locus_tag = "";
	}
	$organism = "organism";
	if (index($linea, $organism) != -1) {
    		$organism = $linea;
    		$organism =~ s/.*organism//;
		$organism =~ s/\,.*//;
	} else {
    		$organism = "";
	}
	$product = "product";
	if (index($linea, $product) != -1) {
    		$product = $linea;
    		$product =~ s/.*product//;
		$product =~ s/\,.*//;
	} else {
    		$product = "";
	}
	$protein_id = "protein_id";
	if (index($linea, $protein_id) != -1) {
    		$protein_id = $linea;
    		$protein_id =~ s/.*protein_id//;
		$protein_id =~ s/\,.*//;
	} else {
    		$protein_id = "";
	}
	$regulatory_class = "regulatory_class";
	if (index($linea, $regulatory_class) != -1) {
    		$regulatory_class = $linea;
    		$regulatory_class =~ s/.*regulatory_class//;
		$regulatory_class =~ s/\,.*//;
	} else {
    		$regulatory_class = "";
	}
	$segment = "segment";
	if (index($linea, $segment) != -1) {
    		$segment = $linea;
    		$segment =~ s/.*segment//;
		$segment =~ s/\,.*//;
	} else {
    		$segment = "";
	}
	$serotype = "serotype";
	if (index($linea, $serotype) != -1) {
    		$serotype = $linea;
    		$serotype =~ s/.*serotype//;
		$serotype =~ s/\,.*//;
	} else {
    		$serotype = "";
	}
	$standard_name = "standard_name";
	if (index($linea, $standard_name) != -1) {
    		$standard_name = $linea;
    		$standard_name =~ s/.*standard_name//;
		$standard_name =~ s/\,.*//;
	} else {
    		$standard_name = "";
	}
	$strain = "strain";
	if (index($linea, $strain) != -1) {
    		$strain = $linea;
    		$strain =~ s/.*strain//;
		$strain =~ s/\,.*//;
	} else {
    		$strain = "";
	}
	$subgroup = "subgroup";
	if (index($linea, $subgroup) != -1) {
    		$subgroup = $linea;
    		$subgroup =~ s/.*subgroup//;
		$subgroup =~ s/\,.*//;
	} else {
    		$subgroup = "";
	}
	$transl_except = "transl_except";
	if (index($linea, $transl_except) != -1) {
    		$transl_except = $linea;
    		$transl_except =~ s/.*transl_except//;
		$transl_except =~ s/\,.*//;
	} else {
    		$transl_except = "";
	}
	$ORGANISMO = "ORGANISMO";
	if (index($linea, $ORGANISMO) != -1) {
    		$ORGANISMO = $linea;
    		$ORGANISMO =~ s/.*ORGANISMO//;
		$ORGANISMO =~ s/\,.*//;
	} else {
    		$ORGANISMO = "";
	}
	$TAXID = "TAXID";
	if (index($linea, $TAXID) != -1) {
    		$TAXID = $linea;
    		$TAXID =~ s/.*TAXID//;
		$TAXID =~ s/\,.*//;
	} else {
    		$TAXID = "";
	}
	$CARACTERISTICA = "CARACTERISTICA";
	if (index($linea, $CARACTERISTICA) != -1) {
    		$CARACTERISTICA = $linea;
    		$CARACTERISTICA =~ s/.*CARACTERISTICA//;
		$CARACTERISTICA =~ s/\,.*//;
	} else {
    		$CARACTERISTICA = "";
	}
	$isjoin = "join";
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
		print OUT "$ORGANISMO,$TAXID,$CARACTERISTICA,$chunk1begin,$chunk2end,$locus_tag,$db_xref,$gene,$product,$protein_id,$codon_start,$note,$mol_type,$host,$country,$collection_date,$strain,$ribosomal_slippage,$isolate,$ribosomal_slippage,$inference,$experiment,$function,$isolation_source,$old_locus_tag,$regulatory_class,$gene_synonym,$collected_by,$lab_host,$genotype,$pseudo,$lat_lon,$cell_line,$acronym,$transl_except,$serotype,$segment,$exception,$citation,$altitude,$subgroup,$standard_name,$identified_by,$group,$culture_collection,$clone\n";
	} else {
		$POSICION =~ s/[^0-9.]+//g;
		$chunk1begin = "BEGIN"."$POSICION";
		$chunk1end = "$POSICION"."END";
		$chunk1begin =~ s/.*BEGIN//;
		$chunk1begin =~ s/\.\..*//;
		$chunk1end =~ s/.*\.\.//;
		$chunk1end =~ s/END.*//;
    		print OUT "$ORGANISMO,$TAXID,$CARACTERISTICA,$chunk1begin,$chunk1end,$locus_tag,$db_xref,$gene,$product,$protein_id,$codon_start,$note,$mol_type,$host,$country,$collection_date,$strain,$ribosomal_slippage,$isolate,$ribosomal_slippage,$inference,$experiment,$function,$isolation_source,$old_locus_tag,$regulatory_class,$gene_synonym,$collected_by,$lab_host,$genotype,$pseudo,$lat_lon,$cell_line,$acronym,$transl_except,$serotype,$segment,$exception,$citation,$altitude,$subgroup,$standard_name,$identified_by,$group,$culture_collection,$clone\n";
	}
} #END OF CONCATENATED GENBANK
close (FILE); #Close output at input end
system ("rm ../results/GenFeatures_locations/$taxon\_peptide_locations.csv");
system ("rm ../results/GenFeatures_locations/$taxon.tmp");
system ("sed -i 's/\=//g' ../results/GenFeatures_locations/$taxon\_features_locations.csv");
system ("sed -i 's/\"//g' ../results/GenFeatures_locations/$taxon\_features_locations.csv");
system ("Rscript Tab_correction.R $taxon");
