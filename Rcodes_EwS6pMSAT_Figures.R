library(dplyr)
library("data.table") 
library(stringr)
library(ggplot2)
library(tidyverse)
library("htmlwidgets")

#################
#    Figures    #
#################

ggplot(cohort.pattern.count, aes(x = hap, y = Freq, fill = color, group = rank))+ 
  geom_bar(stat = "identity") +
  ggtitle ("Microsatellites")+
  xlab("microsatellite") +
  ylab("microsatellite Pattern") +
  geom_text(position = position_stack(vjust=0.5), aes(label = Freq), size=4) +
  scale_fill_manual(values = c("AGAA" = "#ccebc5", 
                               "GAAG" = "#fddaec", 
                               "GAG_" = "#b3cde3" ,
                               "GAGG" = "#fbb4ae" ,
                               "GGA_" = "#fed9a6" ,
                               "GGAA" = "#fff2ae" ,
                               "GGAG" = "#decbe4" ,
                               "GGGA" = "#cccccc")) +
  scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

############### Flipped ###############
## Setting levles
cohort.pattern.count.v2 <- cohort.pattern.count
cohort.pattern.count.v2$hap <- factor(cohort.pattern.count.v2$hap, levels = unique(rev(as.character(cohort.pattern.count.v2$hap))))
cohort.pattern.count.v2$rank <- factor(cohort.pattern.count.v2$rank, levels = unique(rev(as.character(cohort.pattern.count.v2$rank))))

ggplot(cohort.pattern.count.v3, aes(x = hap, y = Freq, fill = color, group = rank))+ 
  geom_bar(stat = "identity") +
  #  ggtitle ("Microsatellites")+
  xlab("Microsatellite") +
  ylab("Microsatellite Allele Length (bp)") +
  geom_text(position = position_stack(vjust=0.5), aes(label = Freq/4), size=4) +
  scale_fill_manual(values = c("AGAA" = "#ccebc5", 
                               "GAAG" = "#fddaec", 
                               "GAG_" = "#b3cde3" ,
                               "GAGG" = "#fbb4ae" ,
                               "GGA_" = "#fed9a6" ,
                               "GGAA" = "#fff2ae" ,
                               "GGAG" = "#decbe4" ,
                               "GGGA" = "#cccccc")) +
  scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.major = element_blank())

dev.off()

# ****************************  #
#     Motif Expansion Label     #
# ****************************  #

ggplot(cohort.pattern.count.v4, aes(x = hap, y = Freq, fill = expansion, group = rank))+ 
  geom_bar(stat = "identity") +
  xlab("Microsatellite") +
  ylab("Microsatellite Allele Length (bp)") +
  #  geom_text(position = position_stack(vjust=0.5), aes(label = Freq), size=2) +
  scale_fill_manual(values = c("GGAA expansion 1" = "#fff2b3", 
                               "GGAA expansion 2" = "#ffee99", 
                               "GGAA expansion 3" = "#ffea80" ,
                               "GGAA expansion 4" = "#ffe14d" ,
                               "AGAA expansion 1" = "#def2d9" ,
                               "AGAA expansion 2" = "#ccebc5" ,
                               "AGAA expansion 3" = "#9bd88d" ,
                               "GGGA expansion 1" = "#cccccc",
                               "GGGA expansion 2" = "#cccccc",
                               "other motifs" = "#fbb4ae")) +
  scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.major = element_blank())

dev.off()

# ************************ #
#     Grouped Bar Plot     #
# ************************ #
input.stackedplot6 <- mSat.information.v3 %>% 
  select(mSat_order, sequence, ews_allele, control_allele, mSat_length) %>% 
  arrange(desc(mSat_length)) %>% 
  select(mSat_order, ews_allele, control_allele)
dat.m3 <- melt(input.stackedplot6, id.vars = "mSat_order")
dat.m3$mSat_order = factor(dat.m3$mSat_order, input.order)
dat.m3$variable = factor(dat.m3$variable, c("control_allele", "ews_allele"))

ggplot(dat.m3, aes(fill=variable, y=value, x=mSat_order)) + 
  geom_bar(position="dodge", stat="identity", alpha = 1.0) +
  scale_fill_brewer(palette="Paired") +
  scale_y_continuous(expand = expansion(mult = c(0.01,0.01)), limits = c(0,100)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank() ) +
  #  ggtitle ("Frequency of EwS case and control")+
  xlab("Microsatellite") +
  ylab("Allele Count") 

# ************************************************** #
#   Length of microsatellite with the odds of EWS    #
# ************************************************** #
### Using geom_point (better) -USE THIS
df <- mSat.information.v3 %>% select(mSat_length, ews_prop, ews_allele, mSat_status, total_allele, GGAAcounts, GGAA_block34)
View(df)
df$ews_percent <- df$ews_prop*100
df$ews_allele_percent <- (df$ews_allele/542)*100
sum(df$total_allele)
df$total_allele_percent <- (df$total_allele/770)*100
df$odds <- (df$ews_allele)/(df$total_allele-df$ews_allele)

### To add a trend line
require(stats)
reg <-lm(ews_percent ~ mSat_length, data = df)
reg
#reg2 <-lrm(ews_percent ~ mSat_length, data = df)
coeff=coefficients(reg)
# Equation of the line : 
eq = paste0("y = ", round(coeff[2],1), "*x + ", round(coeff[1],1))

### Plotting
ggplot(df, aes(mSat_length, ews_percent, size = ews_allele_percent*5)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank()) +
  xlab("Length of Microsatellie (bp)") +
  ylab("Percetage of EwS in Allele") +
  labs(size = "Frequency of Allele")
# scale_color_gradient(low="#4eb3d3", high="#081d58")


# *********************************** #
#    mSat only in case and control    #
# *********************************** #
msat.control <- setdiff(overall.control$mSat_order, overall.case$mSat_order)
msat.case <- setdiff(overall.case$mSat_order, overall.control$mSat_order)

##### CASE
View(cohort.pattern.count.v2)
cohort.pattern.count.v2$hap <- gsub("mSat", "Allele ", cohort.pattern.count.v2$hap)
msat.case <- gsub("mSat", "Allele ", msat.case)
msat.case

cohort.pattern.count.v2$hap <- factor(cohort.pattern.count.v2$hap, levels = unique(rev(as.character(cohort.pattern.count.v4$hap))))
cohort.pattern.count.v2$rank <- factor(cohort.pattern.count.v2$rank, levels = unique(rev(as.character(cohort.pattern.count.v4$rank))))

EWS.only.pattern.count <- cohort.pattern.count.v2 %>% filter(hap %in% msat.case)

ggplot(EWS.only.pattern.count, aes(x = `hap_Allele_w/o space`, y = Freq, fill = color, group = rank))+ 
  geom_bar(stat = "identity") +
  ggtitle ("Microsatellite alleles found only in EwS cases")+
  xlab("Microsatellite") +
  ylab("Microsatellite Length") +
  geom_text(position = position_stack(vjust=0.5), aes(label = Freq), size=4) +
  scale_fill_manual(values = c("AGAA" = "#ccebc5", 
                               "GAAG" = "#fddaec", 
                               "GAG_" = "#b3cde3" ,
                               "GAGG" = "#fbb4ae" ,
                               "GGA_" = "#fed9a6" ,
                               "GGAA" = "#fff2ae" ,
                               "GGAG" = "#decbe4" ,
                               "GGGA" = "#cccccc")) +
  scale_y_continuous(limits = c(0, 160), expand = expansion(mult = c(0,0.1))) +
  theme_bw()+
  coord_flip() +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.major = element_blank()) 

##### CONTROL
msat.control <- gsub("mSat", "Allele ", msat.control)
control.only.pattern.count <- cohort.pattern.count.v2 %>% filter(hap %in% msat.control)
View(control.only.pattern.count)

ggplot(control.only.pattern.count, aes(x = `hap_Allele_w/o space`, y = Freq, fill = color, group = rank))+ 
  geom_bar(stat = "identity") +
  ggtitle ("Microsatellite alleles found only in controls")+
  xlab("Microsatellite") +
  ylab("Microsatellite Length") +
  geom_text(position = position_stack(vjust=0.5), aes(label = Freq), size=4) +
  scale_fill_manual(values = c("AGAA" = "#ccebc5", 
                               "GAAG" = "#fddaec", 
                               "GAG_" = "#b3cde3" ,
                               "GAGG" = "#fbb4ae" ,
                               "GGA_" = "#fed9a6" ,
                               "GGAA" = "#fff2ae" ,
                               "GGAG" = "#decbe4" ,
                               "GGGA" = "#cccccc")) +
  scale_y_continuous(limits = c(0, 160), expand = expansion(mult = c(0,0.1))) +
  theme_bw()+
  coord_flip() +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.major = element_blank())


# *********************************** #
#    Gviz - Visualize genomic data    #
# *********************************** #
library(Gviz)
library(GenomicRanges)

grtrack <- GeneRegionTrack(hg19.ncbiRefSeq.chr6.v6, genome = test.gen,
                           chromosome = chr, name = "Gene Model", 
                           transcriptAnnotation = "symbol",
                           background.panel = "#FFFEDB",
                           background.title = "darkblue")

### Adding DNA bases
library(BSgenome.Hsapiens.UCSC.hg19)
strack <- SequenceTrack(Hsapiens, chromosome = chr)
plotTracks(list(test.itrack, gtrack, grtrack, strack),  from = 4000000, to = 9000000, cex = 0.8)

final.fig <- plotTracks(list(test.itrack, gtrack, grtrack, strack, dTrack),  from = 6300000, to = 7500000, cex = 0.8, type = "histogram")


# ******************************************************* #
#    Suuplemetary Figure: Each motif count in each mSat   #
# ******************************************************* # 
View(mSat.information.v3)

# ***************#
#   GGAA count   #
# ***************#

### only GGAA counts
input.stackedplot.GGAA <- mSat.information.v3 %>% 
  select(mSat_order, mSat_length, GGAAcounts) %>% 
  arrange(desc(mSat_length)) %>%
  select(!mSat_length)

input.stackedplot.GGAA$mSat_order <- gsub("mSat", "Allele ", input.stackedplot.GGAA$mSat_order)
input.order <- gsub("mSat", "Allele ", input.order)
input.stackedplot.GGAA$mSat_order = factor(input.stackedplot.GGAA$mSat_order, input.order)
View(input.stackedplot.GGAA)

ggplot(input.stackedplot.GGAA, aes(y=GGAAcounts,  x=reorder(mSat_order, desc(mSat_order)))) +
  geom_bar(stat="identity", alpha = 1.0) +
  scale_y_continuous(expand = expansion(mult = c(0.01,0.01)), limits = c(0,35)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.major = element_blank()) +
  coord_flip() +
  xlab("Microsatellite") +
  ylab("Number of GGAA motif")

### only AGAA counts
input.stackedplot.AGAA <- mSat.information.v3 %>% 
  select(mSat_order, AGAAcounts) 

input.stackedplot.AGAA$mSat_order <- gsub("mSat", "Allele ", input.stackedplot.AGAA$mSat_order)
input.stackedplot.AGAA$mSat_order = factor(input.stackedplot.AGAA$mSat_order, input.order)
View(input.stackedplot.AGAA)

ggplot(input.stackedplot.AGAA, aes(y=AGAAcounts,  x=reorder(mSat_order, desc(mSat_order)))) + 
  geom_bar(stat="identity", alpha = 1.0) +
  scale_y_continuous(expand = expansion(mult = c(0.01,0.01)), limits = c(0,35)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.major = element_blank()) +
  coord_flip() +
  xlab("Microsatellite") +
  ylab("Number of AGAA motif")

### only GGGA counts
input.stackedplot.GGGA <- mSat.information.v3 %>% 
  select(mSat_order, GGGAcounts) 

input.stackedplot.GGGA$mSat_order <- gsub("mSat", "Allele ", input.stackedplot.GGGA$mSat_order)
input.stackedplot.GGGA$mSat_order = factor(input.stackedplot.GGGA$mSat_order, input.order)
View(input.stackedplot.GGGA)

ggplot(input.stackedplot.GGGA, aes(y=GGGAcounts,  x=reorder(mSat_order, desc(mSat_order)))) + 
  geom_bar(stat="identity", alpha = 1.0) +
  scale_y_continuous(expand = expansion(mult = c(0.01,0.01)), limits = c(0,35)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.major = element_blank()) +
  coord_flip() +
  xlab("Microsatellite") +
  ylab("Number of GGGA motif")

######################################################
#     Supp Fig# - 3 SNPs haplotype and GGAA count    #
######################################################
### Counts with GGAA freq
df.aligned.patternCounts.v9$rep_pattern <- gsub("_", "", df.aligned.patternCounts.v9$pattern1)
ggplot(df.aligned.patternCounts.v9, aes(x=rep_pattern, y=freq_GGAA)) + 
  geom_violin(trim = FALSE) +
  geom_jitter(position=position_jitter(0.05), alpha = 0.7, aes(colour = Disease.Group)) +
  scale_color_manual(values=c("#ff6666", "#00CC99")) +
  xlab("Haplotype") +
  ylab("GGAA motif count") +
  labs(color = "EwS status") +
  ylim(5,40) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.major = element_blank())

df.aligned.patternCounts.v10$rep_pattern <- factor(df.aligned.patternCounts.v9$rep_pattern,
                                                   levels = c("ATT", "ATC", "AGT", "GTT", "GTC", "GGC", "GGT"))

ggplot(df.aligned.patternCounts.v10, aes(x=rep_pattern, y=freq_GGAA)) + 
  geom_violin(trim = FALSE) +
  geom_jitter(position=position_jitter(0.05), alpha = 0.7, aes(colour = Disease.Group)) +
  scale_color_manual(values=c("#ff6666", "#00CC99")) +
  xlab("Haplotype") +
  ylab("GGAA motif count") +
  labs(color = "EwS status") +
  ylim(5,40) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.major = element_blank())

############################################################################
# Please contact Olivia Lee (olivia.lee@nih.gov) if you have any questions about the codes and details.
############################################################################
