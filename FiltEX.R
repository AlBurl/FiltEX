library(dplyr)
library(writexl)
library(stringr)
library(tidyr)
library(readr)


#user input file
input.file <- readline(prompt="Enter exact wANNOVAR file name: ")
input.folder <- "S:/ICH_IIIP_IIR_GeneHunting/FILTAR/input/"
input.path <- paste(input.folder, input.file, ".csv", sep = "")

#import WANNOVAR CSV file and assign to full_data
full_data <- read.csv(input.path,fill=TRUE,header=TRUE,stringsAsFactors = FALSE)

#remove unneccessary columns and assign to tidy_data
tidy_data <- full_data[c(1:7,9:11,17,46,52,58,110,127,132,136,137)]


#change all Blanks and . to '0' in gene frequency columns (1000G)
dots_1000G<-tidy_data %>% mutate(X1000G_ALL = replace(X1000G_ALL,X1000G_ALL == ".", 0))
zeros_1000G<-dots_1000G %>% mutate(X1000G_ALL = replace(X1000G_ALL,X1000G_ALL == "", 0))

#change all Blanks and . to '0' in gene frequency columns (ExAC)
dots_ExAC<-zeros_1000G %>% mutate(ExAC_Freq = replace(ExAC_Freq,ExAC_Freq == ".", 0))
zeros_ExAC<-dots_ExAC %>% mutate(ExAC_Freq = replace(ExAC_Freq,ExAC_Freq == "", 0))

#change all Blanks and . to '0' in gene frequency columns (gnomad)
dots_gnomad<-zeros_ExAC %>% mutate(gnomAD_exome_ALL = replace(gnomAD_exome_ALL,gnomAD_exome_ALL == ".", 0))
zeros_gnomad<-dots_gnomad %>% mutate(gnomAD_exome_ALL = replace(gnomAD_exome_ALL,gnomAD_exome_ALL == "", 0))


#filter out synonymous variants 
NS_filt_data <- filter(zeros_gnomad, ExonicFunc.refGene!="synonymous SNV")

#Change frequency columns to numeric type 
NS_filt_numeric <- transform(NS_filt_data, X1000G_ALL = as.numeric(X1000G_ALL), ExAC_Freq = as.numeric(ExAC_Freq), gnomAD_exome_ALL = as.numeric(gnomAD_exome_ALL))




#EXOMISER

#set working directory
setwd('C:/Users/sejjag8/Exomiser/exomiser-cli-12.1.0')

#user input file name
input.name <- readline(prompt="Enter exact Exomiser yml file name:")
input.yml <- paste(input.name, ".yml", sep = "")

#create exomiser command
exomiser.command <- paste('java -Xms4g -Xmx6g -jar exomiser-cli-12.1.0.jar --analysis examples/',input.yml, sep='')

#execute exomiser
system(command="cmd.exe", input = exomiser.command)

#create results file names 
#add 1 here if necessary
AD.results.file <- paste(input.name, '_AD.genes.tsv', sep='')
AR.results.file <- paste(input.name, '_AR.genes.tsv', sep='')
MT.results.file <- paste(input.name, '_MT.genes.tsv', sep='')
XD.results.file <- paste(input.name, '_XD.genes.tsv', sep='')
XR.results.file <- paste(input.name, '_XR.genes.tsv', sep='')

AD.results.path <- paste('C:/Users/sejjag8/Exomiser/exomiser-cli-12.1.0/results/', AD.results.file, sep='')
AR.results.path <- paste('C:/Users/sejjag8/Exomiser/exomiser-cli-12.1.0/results/', AR.results.file, sep='')
MT.results.path <- paste('C:/Users/sejjag8/Exomiser/exomiser-cli-12.1.0/results/', MT.results.file, sep='')
XD.results.path <- paste('C:/Users/sejjag8/Exomiser/exomiser-cli-12.1.0/results/', XD.results.file, sep='')
XR.results.path <- paste('C:/Users/sejjag8/Exomiser/exomiser-cli-12.1.0/results/', XR.results.file, sep='')

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
AD.append <- left_join(zeros_gnomad, AD.trunc, by =c('Gene.refGene' = 'GENE_SYMBOL'))
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








#filter out non-rare (freq>1%) variants (1000G, ExAC, gnomad)
All_rare <- filter(NS_filt_numeric, (X1000G_ALL <= 0.01 & ExAC_Freq <= 0.01 & gnomAD_exome_ALL <= 0.01)| Otherinfo.5=='rs3743930'| Otherinfo.5=='rs121908147'| Otherinfo.5=='rs35947132'| Otherinfo.5=='rs4149584'| Otherinfo.5=='rs146885082')


#Using HGNC reference sequences to extract DNA and AA variant from correct transcript 

#import HGNC file 
HGNC <- read.delim("S:/ICH_IIIP_IIR_GeneHunting/FILTAR/databases/HGNC.txt",fill=TRUE,header=TRUE,stringsAsFactors = FALSE)

#match tables based on gene name 
All_rare_HGNC <- left_join(All_rare, HGNC, by =c('Gene.refGene' = 'Approved.symbol'))

#removed unwanted columns and rename columns 
tidy_HGNC <- All_rare_HGNC[c(1:21,27)]

tidy_HGNC_2 <- tidy_HGNC %>%
  rename(HGNC_Match = Approved.name, HGNC_Ref_seq = RefSeq.IDs)

#add columns with positions of correct transcript
ref_positions <- str_locate(tidy_HGNC_2$AAChange.refGene, tidy_HGNC_2$HGNC_Ref_seq)
HGNC_2 <- cbind(tidy_HGNC_2, ref_positions)



#add comma to end of variant column 
HGNC_2$AAChange.refGene <- paste(HGNC_2$AAChange.refGene, ',', sep="")

#extract transcript and variant into new column 'transcript_variant_int'
HGNC_3 <- HGNC_2 %>% 
  mutate(transcript_variant_int = str_sub(HGNC_2$AAChange.refGene, HGNC_2$start, -1))

#remove comma from the end of transcript variant 
HGNC_4 <- HGNC_3 %>%
  mutate(transcript_variant = str_sub(HGNC_3$transcript_variant_int, 1, (regexpr(',', HGNC_3$transcript_variant_int) -1)))

#extract DNA variant into new column 
DNA_variant_df <- HGNC_4 %>%
  mutate(DNA_variant = str_sub(HGNC_4$transcript_variant, (regexpr(':c.', HGNC_4$transcript_variant)+3), (regexpr(':p.', HGNC_4$transcript_variant) -1)))

#extract AA variant into new column 
AA_variant_df <- DNA_variant_df %>%
  mutate(AA_variant = str_sub(HGNC_4$transcript_variant, (regexpr(':p.', HGNC_4$transcript_variant)+3), -1))

#remove unwanted columns
DNA_AA_variant <- AA_variant_df[-c(22:24)]




#Using OMIM to annotate the genes further 

#import OMIM file
OMIM <- read.delim("S:/ICH_IIIP_IIR_GeneHunting/FILTAR/databases/OMIM.txt",sep='\t', skip=3, header=TRUE, stringsAsFactors = FALSE)

#join data with OMIM based on common gene symbol
OMIM_append <- left_join(DNA_AA_variant, OMIM, by =c('Gene.refGene' = 'Gene.Symbols'))

#rename columns 
OMIM_All_rare <- OMIM_append %>%
  rename(OMIM_Chromosome=X..Chromosome, OMIM_Gene_start_position=Genomic.Position.Start, OMIM_Gene_end_position=Genomic.Position.End, OMIM_Cyto.Location=Cyto.Location, OMIM_Computed_Cyto.location=Computed.Cyto.Location, OMIM_number=MIM.Number, OMIM_Gene_name=Gene.Name, OMIM_Comments=Comments, OMIM_Phenotypes=Phenotypes, OMIM_Mouse_Symbol_ID=Mouse.Gene.Symbol.ID)



#separating into hom, het, etc

#create homozygous list
Hom_rare <- filter(OMIM_All_rare, Otherinfo=="hom")

#create heterozygous list
Het_rare <- filter(OMIM_All_rare, Otherinfo=="het")

#split heterozygous list into single and mutliple het lists
Het_count <- Het_rare %>% add_count(Gene.refGene, name = "Gene_count")
Het_single_rare <- filter(Het_count, Gene_count=="1")
Het_multiple_rare <- filter(Het_count, Gene_count!= "1")




#filter list using virtual gene panel 

#import gene list of virtual panel 
IP_gene_list <- read.csv("S:/ICH_IIIP_IIR_GeneHunting/FILTAR/panel_gene_lists/IP_5.csv",fill=TRUE,header=FALSE,stringsAsFactors = FALSE)
VIP_gene_list <- read.csv("S:/ICH_IIIP_IIR_GeneHunting/FILTAR/panel_gene_lists/VIP_5_info.csv",fill=TRUE,header=TRUE,stringsAsFactors = FALSE)
NIP_gene_list <- read.csv("S:/ICH_IIIP_IIR_GeneHunting/FILTAR/panel_gene_lists/NIP.csv",fill=TRUE,header=FALSE,stringsAsFactors = FALSE)

#add column to 'All_rare' dataframe to check gene presence in panels 
IP_panel_check<-OMIM_All_rare %>% mutate(Gene_present_in_IP = (OMIM_All_rare$Gene.refGene %in% IP_gene_list$V1))
VIP_panel_check<-OMIM_All_rare %>% mutate(Gene_present_in_VIP = (OMIM_All_rare$Gene.refGene %in% VIP_gene_list$Gene))
NIP_panel_check<-OMIM_All_rare %>% mutate(Gene_present_in_NIP = (OMIM_All_rare$Gene.refGene %in% NIP_gene_list$V1))


#filter variants present in panel gene lists
IP_filt<- filter(IP_panel_check, Gene_present_in_IP== "TRUE")
VIP_filt<- filter(VIP_panel_check, Gene_present_in_VIP== "TRUE")
NIP_filt<- filter(NIP_panel_check, Gene_present_in_NIP== "TRUE")

#Add info to VIP filtered list 
VIP_info <- left_join(VIP_filt, VIP_gene_list, by =c('Gene.refGene' = 'Gene'))

#Rearrange columns in VIP filtered list
RoR_VIP <-VIP_info[,c(7,40:43,45,24,25,9,16,12:14,10,11,15,1:7,17:21)]


#Create stats table 
Number_of_variants <- c(nrow(tidy_data),nrow(NS_filt_data),nrow(OMIM_All_rare),nrow(Het_rare),nrow(Het_single_rare),nrow(Het_multiple_rare),nrow(Hom_rare),nrow(IP_filt),nrow(VIP_filt),nrow(NIP_filt), nrow(exomiser.filtered))
Filter_step <- c('All','All non-synonymous','All rare, non-synonymous','All heterozygous rare','Single heterozygous rare','Multiple heterozygous rare','All homozygous rare','IP filtered','VIP filtered','NIP filtered','Exomiser')
Stats<-data.frame(Filter_step,Number_of_variants)



#export master excel sheet
sheets <- list("Unfiltered"=full_data, "FILTAR Stats"=Stats, "All rare"=OMIM_All_rare, "Hom rare"=Hom_rare, "Het rare"=Het_rare, "Het single rare"=Het_single_rare, "Het multiple rare"=Het_multiple_rare, "IP (all) filtered"=IP_filt, "VIP filtered"=VIP_filt, 'VIP_RoR'=RoR_VIP, "NIP filtered"=NIP_filt, 'Exomiser'=exomiser.ordered.1) 
output.folder <- "S:/ICH_IIIP_IIR_GeneHunting/FILTAR/output/"
output.path <- paste(output.folder, "FILTAR_", input.file, ".xlsx")
write_xlsx(sheets, output.path)

