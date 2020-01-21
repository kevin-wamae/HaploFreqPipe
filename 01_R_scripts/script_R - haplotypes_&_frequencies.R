##############################################################################################################
# An R.script for generating haplotype frequencies from amino acid sequence alignments
# The script takes as input an alignment in Phylip format (line 22)
##############################################################################################################
#clear environment
rm(list = ls())
#-------------------------------------------------------------------------------------------------------------
#load packages
library(data.table) # fast data importation using fread
library(tidyverse)  # data aggregation

#load functions
source("02_R_functions/function_R - identical_bases_per_locus.R")   #checking for polymorphic loci
source("02_R_functions/function_R - merge_data_frames.R")           #merging dataframes
source("02_R_functions/function_R - snpFreq_per_locus.R")           #snp freq per column/locus

#*************************************************************************************************************
#load data 
#*************************************************************************************************************
#skip line 1 as well as empty lines and select colums 1 and 2
fastaIN <- fread("03_data_input/example_alignment_file.phy", 
                 skip = 1, #skip first line
                 header = F, #no headers
                 blank.lines.skip = T) #skip any blank lines

#split sequence column into multiple columns
fastaIN_seqs <- fastaIN[,tstrsplit(V2, "")]  #require data.table

#select sequence names
fastaIN_names <- fastaIN[,.(V1)]

#bind the two datasets and remove temporary files
fastaIN <- cbind(fastaIN_names, fastaIN_seqs); rm(fastaIN_names,fastaIN_seqs)

#remame column 1 as it matches column 2
names(fastaIN)[1] <- c("V0")

#convert to data.frame
fastaIN <- as.data.frame(fastaIN)

#-------------------------------------------------------------------------------------------------------------
#replace all gaps with NA to exclude them from the analysis 
fastaIN[fastaIN == "-"] <- NA

#create new dataframe with just the sequence IDs, to be used in the subsequent loop
fastaIN_v1 <- as.data.frame(fastaIN[1])

#keep only columns that represent segregating sites
for(colN in 2:dim(fastaIN)[2]) {
  if(dim(table(fastaIN[colN])) != 1 & dim(table(fastaIN[colN])) != 0) { # !=1 means identical bases, != 0 means loci with NA values
    fastaIN_v1 <- as.data.frame(cbind(fastaIN_v1, fastaIN[colN]))
  }
}

#create an empty dataframe, to be used in subsequent loop
freqStats_v1 <- data.frame(freq = as.numeric(),
                        percentage = as.numeric(),
                        base = as.character(),
                        locus = as.character())

#compute sample size and Func_SnpFrequency for each nucleotide per locus
for(i in 2:dim(fastaIN_v1)[2]) {
  freqStats_v2 <- Func_SnpFreq(fastaIN_v1[,i])
  freqStats_v2 <- as.data.frame(freqStats_v2)
  freqStats_v2 <- freqStats_v2 %>% mutate(base = row.names(freqStats_v2),
                                          locus = colnames(fastaIN_v1[i]))
  freqStats_v1 <- rbind(freqStats_v1, freqStats_v2)
} 

#clear environment except snpFreq and fastaData_v1
rm(list=setdiff(ls(), c("fastaIN_v1", "freqStats_v1")))

#remove sequences that are not represented across all locus
fastaIN_v1 <- fastaIN_v1[complete.cases(fastaIN_v1), ]

#-------------------------------------------------------------------------------------------------------------
#summary
#-------------------------------------------------------------------------------------------------------------
# a table of haplotype frequencies
hapFreq <- fastaIN_v1 %>% 
  unite("haplotypes", 2:dim(fastaIN_v1)[2], sep = "", remove = TRUE) %>% #merge residues into a string of haplotypes
  group_by(haplotypes) %>%
  dplyr::summarize(n=n()) %>%
  mutate(freq = (n / sum(n)*100)) %>%
  mutate(freq = round(freq,2)) %>%
  mutate(freq_perc = paste(n, " [",freq,"]", sep = "")) %>%
  select(haplotypes, freq_perc, freq) %>%
  arrange(desc(freq)) %>%
  rename(Haplotypes = haplotypes,
         "n (Freq)" = freq_perc,
         Freq = freq)

#write data
fwrite(hapFreq, "04_data_output/example_haplotype_frequencies.csv")


