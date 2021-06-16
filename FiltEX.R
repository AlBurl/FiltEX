

####################################################################


#****************FILTEX*******************

#Written by Alice Burleigh
#June 2021 


####################################################################

#COMPLETE THIS SECTION BEFORE USE 

#master folder location
master.folder <- 'Enter/folder/master/folder/location/here'

#input csv file name (without .csv)
input.file <- 'Enter csv file name here'
  
#filter out synonymous variants? Enter Y/N. Default = Y.
syn.filt <- 'Y' 

#filter out all variants with population allele frequency greater than? Enter value between 0-1. Default = 0.01. 
MAF<- 0.01

#use custom gene panel? Enter Y/N. Default = N. 
use.panel <- 'N'

#NB if custom gene panel being used this must in csv format in master folder.
#if custom gene panel being used enter file name here (without .csv) 
panel.name <- 'Enter panel csv file name here'


#use Exomiser? Enter Enter Y/N. Default = N.
use.exomiser <- 'N'

#if using Exomiser, enter working directory of Exomiser folder here
exomiser.wd <- 'Enter/exomiser/wd/path/here'

#if using Exomiser, enter the yml file name here (without .yml)
exomiser.input <- 'Enter exact yml file name here'

##################################################################


#LOAD PACKAGES 
#NB if packages not used before, install using install.packages()

library(dplyr)
library(writexl)
library(stringr)
library(tidyr)
library(readr)


##################################################################

#IMPORT WANNOVAR FILE AND TIDY DATA 

#user input file
input.path <- paste(master.folder,'/input/',input.file,".csv",sep = "")

#import WANNOVAR CSV file and assign to full_data
full.data <- read.csv(input.path,fill=TRUE,header=TRUE,stringsAsFactors = FALSE)

#remove unneccessary columns and assign to tidy_data
tidy.data <- full.data[c(1:7,9:11,17,46,52,58,110,127,132,136,137)]

#change all Blanks and . to '0' in gene frequency columns (1000G)
dots.1000G<-tidy.data %>% mutate(X1000G_ALL = replace(X1000G_ALL,X1000G_ALL == ".", 0))
zeros.1000G<-dots.1000G %>% mutate(X1000G_ALL = replace(X1000G_ALL,X1000G_ALL == "", 0))

#change all Blanks and . to '0' in gene frequency columns (ExAC)
dots.ExAC<-zeros.1000G %>% mutate(ExAC_Freq = replace(ExAC_Freq,ExAC_Freq == ".", 0))
zeros.ExAC<-dots.ExAC %>% mutate(ExAC_Freq = replace(ExAC_Freq,ExAC_Freq == "", 0))

#change all Blanks and . to '0' in gene frequency columns (gnomad)
dots.gnomad<-zeros.ExAC %>% mutate(gnomAD_exome_ALL = replace(gnomAD_exome_ALL,gnomAD_exome_ALL == ".", 0))
zeros.gnomad<-dots.gnomad %>% mutate(gnomAD_exome_ALL = replace(gnomAD_exome_ALL,gnomAD_exome_ALL == "", 0))

##################################################################

#FILTER OUT SYNONYMOUS VARIANTS

if (syn.filt == 'Y'){
  NS.filt.data <- filter(zeros.gnomad, ExonicFunc.refGene!="synonymous SNV")
}else{
  NS.filt.data <- zeros.gnomad
}

##################################################################

#FILTER OUT COMMON VARIANTS

#Change frequency columns to numeric type 
NS.filt.numeric <- transform(NS.filt.data, X1000G_ALL = as.numeric(X1000G_ALL), ExAC_Freq = as.numeric(ExAC_Freq), gnomAD_exome_ALL = as.numeric(gnomAD_exome_ALL))

#filter out variants with population frequency higher than defined % (in databases 1000G, ExAC, gnomad)
All.rare <- filter(NS.filt.numeric, (X1000G_ALL <= MAF & ExAC_Freq <= MAF & gnomAD_exome_ALL <= MAF))


##################################################################

#EXTRACT CORRECT DNA AND AMINO ACID VARIANTS USING HGNC DATABASE FOR TRANSCRIPT


#download HGNC database from online
HGNC.download <- 'http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_types/gene_with_protein_product.txt'
HGNC.name <- paste(master.folder,'/HGNC.txt',sep='')
download.file(HGNC.download,HGNC.name,mode='wb')

#import HGNC database
HGNC <- read.delim(HGNC.name,fill=TRUE,header=TRUE,stringsAsFactors = FALSE)

#match tables based on gene name 
All.rare.HGNC <- left_join(All.rare, HGNC, by =c('Gene.refGene' = 'symbol'))

#removed unwanted columns and rename columns 
tidy.HGNC <- All.rare.HGNC[c(1:21,42)]

tidy.HGNC.2 <- tidy.HGNC %>%
  rename(HGNC_Match = name, HGNC_Ref_seq = refseq_accession)

#add columns with positions of correct transcript
ref.positions <- str_locate(tidy.HGNC.2$AAChange.refGene, tidy.HGNC.2$HGNC_Ref_seq)
HGNC.2 <- cbind(tidy.HGNC.2, ref.positions)

#add comma to end of variant column 
HGNC.2$AAChange.refGene <- paste(HGNC.2$AAChange.refGene, ',', sep="")

#extract transcript and variant into new column 'transcript_variant_int'
HGNC.3 <- HGNC.2 %>% 
  mutate(transcript.variant.int = str_sub(HGNC.2$AAChange.refGene, HGNC.2$start, -1))

#remove comma from the end of transcript variant 
HGNC.4 <- HGNC.3 %>%
  mutate(transcript.variant = str_sub(HGNC.3$transcript.variant.int, 1, (regexpr(',', HGNC.3$transcript.variant.int) -1)))

#extract DNA variant into new column 
DNA.variant.df <- HGNC.4 %>%
  mutate(DNA_variant = str_sub(HGNC.4$transcript.variant, (regexpr(':c.', HGNC.4$transcript.variant)+3), (regexpr(':p.', HGNC.4$transcript.variant) -1)))

#extract AA variant into new column 
AA.variant.df <- DNA.variant.df %>%
  mutate(AA_variant = str_sub(HGNC.4$transcript.variant, (regexpr(':p.', HGNC.4$transcript.variant)+3), -1))

#remove unwanted columns
DNA_AA_variant <- AA.variant.df[-c(22:24)]

#################################################################

#SORT INTO HETEROZYGOUS AND HOMOZYGOUS LISTS

#create homozygous list
Hom.rare <- filter(DNA_AA_variant, Otherinfo=="hom")

#create heterozygous list
Het.rare <- filter(DNA_AA_variant, Otherinfo=="het")

#split heterozygous list into single and mutliple het lists
Het.count <- Het.rare %>% add_count(Gene.refGene, name = "Gene.count")
Het.single.rare <- filter(Het.count, Gene.count=="1")
Het.multiple.rare <- filter(Het.count, Gene.count!= "1")


#################################################################


#FILTER USING CUSTOM GENE PANEL 


if (use.panel == 'Y') {
  
  #import gene panel
 panel.path <- paste(master.folder,'/',panel.name,'.csv',sep='')
 panel.gene.list <- read.csv(panel.path,fill=TRUE,header=FALSE,stringsAsFactors = FALSE)
 
 #check variants in filtered gene list against panel list
 panel.check <- DNA_AA_variant %>% mutate(Gene.present.in.panel = (DNA_AA_variant$Gene.refGene %in% panel.gene.list$V1))
 
 #filter for panel genes 
 panel.filt <- filter(panel.check, Gene.present.in.panel== "TRUE")

} else {
  
  panel.filt <- DNA_AA_variant
}

#################################################################


#EXOMISER

if (use.exomiser == 'Y') {
  
  #set working directory
  setwd(exomiser.wd)

  #user input file name
  input.yml <- paste(exomiser.input, ".yml", sep = "")

  #create exomiser command
  exomiser.command <- paste('java -Xms4g -Xmx6g -jar exomiser-cli-12.1.0.jar --analysis examples/',input.yml, sep='')

  #execute exomiser
  system(command="cmd.exe", input = exomiser.command)

  #create results file names 
  AD.results.file <- paste(exomiser.input, '_AD.genes.tsv', sep='')
  AR.results.file <- paste(exomiser.input, '_AR.genes.tsv', sep='')
  MT.results.file <- paste(exomiser.input, '_MT.genes.tsv', sep='')
  XD.results.file <- paste(exomiser.input, '_XD.genes.tsv', sep='')
  XR.results.file <- paste(exomiser.input, '_XR.genes.tsv', sep='')
  
  AD.results.path <- paste(exomiser.wd,'/results/',AD.results.file,sep='')
  AR.results.path <- paste(exomiser.wd,'/results/',AR.results.file,sep='')
  MT.results.path <- paste(exomiser.wd,'/results/',MT.results.file,sep='')
  XD.results.path <- paste(exomiser.wd,'/results/',XD.results.file,sep='')
  XR.results.path <- paste(exomiser.wd,'/results/',XR.results.file,sep='')
  
  #import results tsv file
  AD.exomiser.results <-read_tsv(AD.results.path)
  AR.exomiser.results <-read_tsv(AR.results.path)
  MT.exomiser.results <-read_tsv(MT.results.path)
  XD.exomiser.results <-read_tsv(XD.results.path)
  XR.exomiser.results <-read_tsv(XR.results.path)
  
  #remove hashtag from column header
  header_change_func <-function(x){
    x %>%
      rename('GENE_SYMBOL' = '#GENE_SYMBOL')
  }
  AD.exomiser.results.1 <- header_change_func(AD.exomiser.results)
  AR.exomiser.results.1 <- header_change_func(AR.exomiser.results)
  MT.exomiser.results.1 <- header_change_func(MT.exomiser.results)
  XD.exomiser.results.1 <- header_change_func(XD.exomiser.results)
  XR.exomiser.results.1 <- header_change_func(XR.exomiser.results)
  
  
  #add inheritance column and combine
  AD.exomiser.results.2 <- AD.exomiser.results.1 %>%
    mutate(Inheritance='AD')
  AR.exomiser.results.2 <- AR.exomiser.results.1 %>%
    mutate(Inheritance='AR')
  MT.exomiser.results.2 <- MT.exomiser.results.1 %>%
    mutate(Inheritance='MT')
  XD.exomiser.results.2 <- XD.exomiser.results.1 %>%
    mutate(Inheritance='XD')
  XR.exomiser.results.2 <- AD.exomiser.results.1 %>%
    mutate(Inheritance='XR')
  
  all.results<- rbind(AD.exomiser.results.2,AR.exomiser.results.2,MT.exomiser.results.2,XD.exomiser.results.2,XR.exomiser.results.2)
  
  
  #truncate the Exomiser results tables 
  AD.trunc <- select(AD.exomiser.results.2, GENE_SYMBOL, EXOMISER_GENE_COMBINED_SCORE)
  AR.trunc <- select(AR.exomiser.results.2, GENE_SYMBOL, EXOMISER_GENE_COMBINED_SCORE)
  MT.trunc <- select(MT.exomiser.results.2, GENE_SYMBOL, EXOMISER_GENE_COMBINED_SCORE)
  XD.trunc <- select(XD.exomiser.results.2, GENE_SYMBOL, EXOMISER_GENE_COMBINED_SCORE)
  XR.trunc <- select(XR.exomiser.results.2, GENE_SYMBOL, EXOMISER_GENE_COMBINED_SCORE)
  
  #append and rename the Exomiser combined scores for each inheritance
  AD.append <- left_join(zeros.gnomad, AD.trunc, by =c('Gene.refGene' = 'GENE_SYMBOL'))
  AD.append<- rename(AD.append, AD_EXOMISER_SCORE = EXOMISER_GENE_COMBINED_SCORE)
  
  AR.append <- left_join(AD.append, AR.trunc, by =c('Gene.refGene' = 'GENE_SYMBOL'))
  AR.append<- rename(AR.append, AR_EXOMISER_SCORE = EXOMISER_GENE_COMBINED_SCORE)
  
  MT.append <- left_join(AR.append, MT.trunc, by =c('Gene.refGene' = 'GENE_SYMBOL'))
  MT.append<- rename(MT.append, MT_EXOMISER_SCORE = EXOMISER_GENE_COMBINED_SCORE)
  
  XD.append <- left_join(MT.append, XD.trunc, by =c('Gene.refGene' = 'GENE_SYMBOL'))
  XD.append<- rename(XD.append, XD_EXOMISER_SCORE = EXOMISER_GENE_COMBINED_SCORE)
  
  XR.append <- left_join(XD.append, XR.trunc, by =c('Gene.refGene' = 'GENE_SYMBOL'))
  XR.append<- rename(XR.append, XR_EXOMISER_SCORE = EXOMISER_GENE_COMBINED_SCORE)
  
  #add exomiser rationale from AD dataframe
  rationale.trunc <- select(AD.exomiser.results.2, GENE_SYMBOL, HUMAN_PHENO_EVIDENCE, MOUSE_PHENO_EVIDENCE, FISH_PHENO_EVIDENCE, HUMAN_PPI_EVIDENCE, MOUSE_PPI_EVIDENCE, FISH_PPI_EVIDENCE)
  exomiser.rationale <- left_join(XR.append, rationale.trunc, by =c('Gene.refGene' = 'GENE_SYMBOL'))
  
  #filter out genes with exomiser score
  exomiser.filtered <- filter(exomiser.rationale, AD_EXOMISER_SCORE != '')
  
  #order exomiser results 
  exomiser.ordered <- arrange(exomiser.filtered, -AD_EXOMISER_SCORE)
  exomiser.ordered.1 <- exomiser.ordered[,c(1:7,20:30,8:19)]
  
  
  
  
  ###############################################################
  
  #CREATE TABLE OF VARIANT NUMBERS WITH EXOMISER 
  
  Number_of_variants <- c(nrow(tidy.data),nrow(NS.filt.data),nrow(DNA_AA_variant),nrow(Het.rare),nrow(Het.single.rare),nrow(Het.multiple.rare),nrow(Hom.rare),nrow(panel.filt), nrow(exomiser.ordered.1))
  Filter_step <- c('All','All non-synonymous','All rare, non-synonymous','All heterozygous rare','Single heterozygous rare','Multiple heterozygous rare','All homozygous rare','Panel filtered','Exomiser')
  Stats<-data.frame(Filter_step,Number_of_variants)
  
  ###############################################################
  
  #EXPORT RESULTS EXCEL WORKBOOK WITH EXOMISER
  
  sheets <- list("Unfiltered"=full.data, "FILTAR Stats"=Stats, "All rare"=DNA_AA_variant, "Hom rare"=Hom.rare, "Het rare"=Het.rare, "Het single rare"=Het.single.rare, "Het multiple rare"=Het.multiple.rare, "Panel filtered"=panel.filt, 'Exomiser'=exomiser.ordered.1) 
  output.folder <- paste(master.folder,'/output/',sep='')
  output.path <- paste(output.folder, "FILTAR_", input.file, ".xlsx")
  write_xlsx(sheets, output.path)
  
} else {
  
  
###############################################################

#CREATE TABLE OF VARIANT NUMBERS 

Number_of_variants <- c(nrow(tidy.data),nrow(NS.filt.data),nrow(DNA_AA_variant),nrow(Het.rare),nrow(Het.single.rare),nrow(Het.multiple.rare),nrow(Hom.rare),nrow(panel.filt))
Filter_step <- c('All','All non-synonymous','All rare, non-synonymous','All heterozygous rare','Single heterozygous rare','Multiple heterozygous rare','All homozygous rare','Panel filtered')
Stats<-data.frame(Filter_step,Number_of_variants)

###############################################################

#EXPORT RESULTS EXCEL WORKBOOK

sheets <- list("Unfiltered"=full.data, "FILTAR Stats"=Stats, "All rare"=DNA_AA_variant, "Hom rare"=Hom.rare, "Het rare"=Het.rare, "Het single rare"=Het.single.rare, "Het multiple rare"=Het.multiple.rare, "Panel filtered"=panel.filt) 
output.folder <- paste(master.folder,'/output/',sep='')
output.path <- paste(output.folder, "FILTAR_", input.file, ".xlsx")
write_xlsx(sheets, output.path)}

###############################################################


