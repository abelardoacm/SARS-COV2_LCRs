library("tidyverse")
library("viridis")

PosArgs <- as.character(commandArgs(trailingOnly = TRUE))
#PosArgs <- "alphas"
variante = PosArgs[1]
infile <- paste("../results/seg/",variante,"_complexity.csv",sep="")

#object with seg output
var_comp.df <- read.table(infile, sep = ',', header = TRUE, na.strings = "", fill=TRUE,quote = "")

#estimate absolute end and begining positions of segs
var_comp.df$length_in_aa <- var_comp.df$seg_end - var_comp.df$seg_begin
var_comp.df$length_in_nc <- var_comp.df$length_in_aa * 3
var_comp.df$abs_nc_beg <- var_comp.df$origin_beg + ((var_comp.df$seg_begin*3)-3)
var_comp.df$abs_nc_end <- var_comp.df$abs_nc_beg+var_comp.df$length_in_nc
var_comp.df$abs_aa_beg <- round(var_comp.df$abs_nc_beg/3)
var_comp.df$abs_aa_end <- var_comp.df$abs_aa_beg + var_comp.df$length_in_aa

#remove duplicated sequences for the same virus
# Remove duplicated rows based on 
# Sepal.Length and Petal.Width
f_var_comp.df <- var_comp.df %>% distinct(secuence, viral_species, .keep_all = TRUE)

#sample size
sasi<-length(f_var_comp.df$viral_species)

#aminoacid positions vector
aa_pvec <- c(1:max(f_var_comp.df$abs_aa_end))

#define function to get complexity value of a position
getComplexity_moda <- function(x, data) {
  tmp <- data %>%
    filter(abs_aa_beg <= x, x <= abs_aa_end)
  return(mean(tmp$complexity))
}

getProduct_moda <- function(x, data) {
  tmp <- data %>%
    filter(abs_aa_beg <= x, x <= abs_aa_end)
  return(names(sort(summary(as.factor(tmp$product)), decreasing=T)[1]))
}

x<-aa_pvec
positional_complexity_moda <- sapply(aa_pvec,getComplexity_moda,data=f_var_comp.df)
positional_product_moda <- sapply(aa_pvec,getProduct_moda,data=f_var_comp.df)

#fill null values
meancomp <- mean(var_comp.df$complexity)
for(i in 1:length(x)){
  if (length(positional_complexity_moda[[i]])==0){
    positional_complexity_moda[[i]]<-meancomp
  }
  if (length(positional_product_moda[[i]])==0){
    positional_product_moda[[i]]<-"Noncoding"
  }
}

#create data frame with positions, moda complexity and product
pos_comp.df <- data.frame(Posicion_aa=x)

for (i in 1:length(x)){
  pos_comp.df$Complexity[i] <- positional_complexity_moda[[i]][1]
  pos_comp.df$Product[i]<-positional_product_moda[[i]][1]
}

pos_comp.df$Complexity<-as.numeric(pos_comp.df$Complexity)
pos_comp.df$Product<-factor(pos_comp.df$Product, levels=unique(pos_comp.df$Product))

#PLOTIME
proteinas <- unique(pos_comp.df$Product)
if("Noncoding" %in% proteinas){
  proteinas <- proteinas[-which(proteinas=="Noncoding")]
}
nombres <-c()
empiezan <- c()
acaban <- c()
for (p in proteinas){
  nombres<-append(nombres,p,after = length(nombres))
  empiezan<-append(empiezan,min(pos_comp.df$Posicion_aa[which(pos_comp.df$Product==p)]),after = length(empiezan))
  acaban<-append(acaban,max(pos_comp.df$Posicion_aa[which(pos_comp.df$Product==p)]),after = length(acaban))
}
porproteina.df <- data.frame(Protein=nombres, P_start=as.numeric(empiezan),P_end=as.numeric(acaban))
porproteina.df$Protein <- factor(porproteina.df$Protein, levels=unique(porproteina.df$Protein))

baseplot <- ggplot() +
  geom_line(data=pos_comp.df, aes(x=Posicion_aa, y=Complexity,group=1))+
  scale_y_continuous(limits=c(-2,6),name="seg estimated complexity")+
  geom_rect(data=porproteina.df, mapping=aes(xmin=P_start, xmax=P_end, ymin=-2, ymax=0, fill=Protein), alpha=1)+
  #geom_text(data=porproteina.df, aes(x=(P_start + P_end)/2, y=-1, label=str_wrap(Protein, width = 1)), size=2.5, angle =90)+
  geom_segment(data = porproteina.df , aes(x = P_start, y = 0, xend = P_start, yend = 6), color = "gray",size=0.3)+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none",
        axis.line = element_line(colour = "gray"))




system("mkdir -p ../results/Mean_complots/")
complot <- baseplot + geom_text(data=porproteina.df, aes(x=(P_start + P_end)/2, y=-1, label=str_wrap(Protein, width = 1)), size=2.5, angle =90)
outcomplot <- paste("../results/Mean_complots/",variante,"_proteome_complot.tiff",sep="")
tiff(outcomplot, units="in", width=18, height=7, res=300) #tiff resolution parameters
complot
dev.off()


# Splitted by protein
for (p in proteinas){
  p2p = p
  outprotcomplot <- paste("../results/Mean_complots/",variante,"_",p2p,"_complot.tiff",sep="")
  protcomplot <- baseplot + scale_x_continuous(limits=c(min(which(pos_comp.df$Product==p2p)),max(which(pos_comp.df$Product==p2p)))) + geom_text(data=porproteina.df, aes(x=(P_start + P_end)/2, y=-1, label=str_wrap(Protein, width = 1)), size=2.5)
  tiff(outprotcomplot, units="in", width=18, height=7, res=300) #tiff resolution parameters
  print(protcomplot)
  dev.off()
}







### Descomponer las LCRs en elementos y ver cuantos comparten, pueden ser kmeros
