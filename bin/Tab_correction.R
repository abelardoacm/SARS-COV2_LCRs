#!/usr/bin/env Rscript
#
# task: this script mask my lack of perl skills, it just deletes the first cell 
# of anytaxon_features_locations.csv and shift cells up. Then repeats and appends
# the last taxid
#
#
# Tab_correction.R
#
# reading input
PosArgs <- as.character(commandArgs(trailingOnly = TRUE))
anytaxon = PosArgs[1]
tab_location <- paste("../results/GenFeatures_locations/",anytaxon,"_features_locations.csv", sep='')
tabla_features <- read.csv(tab_location, header = TRUE)
# deleting first empty cell and shift cells up (hardcoded)
taxids <- tabla_features$NCBI_taxid[-1]
taxids[(length(taxids)+1)] <- taxids[length(taxids)]
# restoring original but corrected table
tabla_features$NCBI_taxid <- taxids
# delete duplicates
tabla_features <- tabla_features[!duplicated(tabla_features), ]
# saving output
write.csv(tabla_features, file = tab_location, row.names = FALSE)

