library(dplyr)
library("data.table") 
library(stringr)
library(ggplot2)
library(tidyverse)
library("htmlwidgets")

###########################
#   Reading FASTQ files   #
###########################

filelist <- list.files(pattern="*.fastq",full.names = T)
txtlist <- list()

for (i in 1:length(filelist))
{
  txtlist [[i]]<-read.table(filelist[i], quote = "\"", comment.char = "")
}

View(filelist)

## Manually inspected top 2 haplotypes in each sample
hap_freq15

## Patient annotation data
EWS.sample.info
View(EWS.sample.info)

## Merge two information
haplotype.patient.info <- merge(hap_freq15, EWS.sample.info, by.x = "sampleID", by.y = "preferredsampleid")
dim(haplotype.patient.info)
dim(hap_freq15)

###########################
#    Halotype Frequecy    #
###########################

### Haplotype frequency Chart
unique.hap <- haplotype.patient.info %>%
  count(sequence) %>%
  mutate(prop = prop.table(n)) %>%
  arrange(desc(n))
View(unique.hap)

unique.hap$number <-1:50
unique.hap$number <- sub("^", "hap", unique.hap$number)
## Assiginig the level
unique.hap$number_level = factor(unique.hap$number, levels=c(paste0('hap',1:50)))

### number of samples containing each 50 frequent haplotype
df_total <- data.frame()
for (i in 1:length(fastqlist))
{
  df.new <-fastqlist [[i]][seq(2, nrow(txtlist[[i]]), 4), ]
  df.new <- as.data.frame(df.new)
  setnames(df.new, "sequence") 
  df.new <- na.omit(df.new)
  tmp.new <- df.new %>%
    count(sequence) %>%
    mutate(prop = prop.table(n)) %>%
    arrange(desc(n))
  tmp.new$sampleID <- names(fastqlist[i])
  
  df_total <- rbind(df_total,tmp.new)
}
View(df_total)
fastq.total <- df_total

############################################################################
# Please contact Olivia Lee (olivia.lee@nih.gov) if you have any questions about the codes and details.
############################################################################