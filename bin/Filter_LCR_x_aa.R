library("stringr")
PosArgs <- as.character(commandArgs(trailingOnly = TRUE))
#PosArgs <- "kappas"
dataset_name = PosArgs[1]
infile <- paste("../results/seg/",dataset_name,"_complexity.csv",sep="")
seg_2_correct <- read.table(infile, sep = ',', header = TRUE, na.strings = "", fill=TRUE,quote = "")
condition_one <- grepl("x",seg_2_correct$secuence, ignore.case = TRUE)
condition_two <- str_detect(seg_2_correct$secuence, "^[:lower:]+$")
decision.df <- data.frame(c1=condition_one,c2=condition_two,amiout="")
for (i in 1:length(condition_one)){
  if(condition_one[i]==TRUE || condition_two[i]==TRUE){
    if(condition_one[i]==condition_two[i]){
      decision.df$amiout[i] <- TRUE
    }
  }
}
corrected_seg <- seg_2_correct[which(decision.df$amiout != TRUE),]
write.csv(corrected_seg,infile,row.names = FALSE)


