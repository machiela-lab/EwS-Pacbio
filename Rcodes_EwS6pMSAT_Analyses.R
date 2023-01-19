library(stringr)
library("dplyr")
library("ggpubr")
library("rms")
library(broom)

###### All samples with 3 microsatellite alleles
####################################
#   Total Cohort (case+control)   #
####################################
EWS.case <- haplotype.patient.info.doubled %>% filter(Affection.Status== "A")
EWS.control <- haplotype.patient.info.doubled %>% filter(Affection.Status == "U")
EWS.case %>% group_by(GGAAcounts) %>% tally()
EWS.control %>% group_by(GGAAcounts) %>% tally()

summary(EWS.case$GGAAcounts)
summary(EWS.control$GGAAcounts)
sum(EWS.case$GGAAcounts)
sum(EWS.control$GGAAcounts)
summary(EWS.case$AGAAcounts)
summary(EWS.control$AGAAcounts)
summary(EWS.case$GGGAcounts)
summary(EWS.control$GGGAcounts)
summary(EWS.case$motifcounts)
summary(EWS.control$motifcounts)


#**********************************************#
#   Normalization test - shapiro wilk test    #
#**********************************************#
shapiro.test(EWS.case$GGAAcounts)

#******************************************************************#
#   non-parametric test - Mann Whitney test (two sample t-test)   #
#****************************************************************#
wilcox.test(EWS.case$GGAAcounts,EWS.control$GGAAcounts)
wilcox.test(EWS.case$AGAAcounts,EWS.control$AGAAcounts)
wilcox.test(EWS.case$GGGAcounts,EWS.control$GGGAcounts)

#**************#
#   T test    #
#*************#
t.test(EWS.case$GGAAcounts,EWS.control$GGAAcounts)
t.test(EWS.case$AGAAcounts,EWS.control$AGAAcounts)
t.test(EWS.case$GGGAcounts,EWS.control$GGGAcounts)

haplotype.patient.info.doubled %>% 
  group_by(Affection.Status, mSat_status) %>%
  tally()

tbl <- matrix(c(90,452,24,204), ncol =2)
chisq.test(tbl)
prop.test(x = c(90, 24), n = c((90+452), (24+204)))

###############################
#   Microsatellite Analysis   #
###############################
####################################################
#   logistic regression - mSat length increase    #
####################################################
## Question 1: the correlation between the total length of mSat and the EWS risk
df.glm <- overall.df.v2 %>% select(mSat_length.x, Affection.Status)
df.glm$Affection.Status = factor(df.glm$Affection.Status, levels = c("U", "A"))
fit <- glm(Affection.Status ~ mSat_length.x, data=df.glm, family=binomial())
summary(fit)

mod1b <- lrm(Affection.Status ~ mSat_length.x, data=df.glm)
print(mod1b)

## Question 2: the correlation between the total number of GGAA repeat and the EWS risk
df.glm <- overall.df.v2 %>% select(GGAAcounts, Affection.Status)
df.glm$Affection.Status = factor(df.glm$Affection.Status, levels = c("U", "A"))
fit <- glm(Affection.Status ~ GGAAcounts, data=df.glm, family=binomial())
summary(fit) 

## Question 2b: the correlation between the total number of GGAA repeat and the EWS risk in the joined model
df.glm <- overall.df.v2 %>% select(GGAAcounts, AGAAcounts, GGGAcounts, Affection.Status)
df.glm$Affection.Status = factor(df.glm$Affection.Status, levels = c("U", "A"))
fit <- glm(Affection.Status ~ GGAAcounts+AGAAcounts+GGGAcounts, data=df.glm, family=binomial())
summary(fit) 

########## CORRELATION ##########
cor(overall.df.v2$GGAAcounts, overall.df.v2$AGAAcounts)
cor(overall.df.v2$GGAAcounts, overall.df.v2$GGGAcounts)
cor(overall.df.v2$AGAAcounts, overall.df.v2$GGGAcounts)

## Question 3: the correlation between the length of consecutive GGAA repeat (in each mSat)and the EWS risk
df.glm <- overall.df.v2 %>% select(consecutive_GGAA_max, Affection.Status)
df.glm$Affection.Status = factor(df.glm$Affection.Status, levels = c("U", "A"))
fit <- glm(Affection.Status ~ consecutive_GGAA_max, data=df.glm, family=binomial())
summary(fit) 

mod1b <- lrm(Affection.Status ~ consecutive_GGAA_max, data=df.glm)
print(mod1b)

########## AGAA block in the pattern ##########
df.glm <- overall.df.v2 %>% select(AGAAcounts, Affection.Status)
df.glm$Affection.Status = factor(df.glm$Affection.Status, levels = c("U", "A"))
fit <- glm(Affection.Status ~ AGAAcounts, data=df.glm, family=binomial())
summary(fit) 

mod1b <- lrm(Affection.Status ~ AGAAcounts, data=df.glm)
print(mod1b)

########## GGGA block in the pattern ##########
df.glm <- overall.df.v2 %>% select(GGGAcounts, Affection.Status)
df.glm$Affection.Status = factor(df.glm$Affection.Status, levels = c("U", "A"))
fit <- glm(Affection.Status ~ GGGAcounts, data=df.glm, family=binomial())
summary(fit) # display results

mod1b <- lrm(Affection.Status ~ GGGAcounts, data=df.glm)
print(mod1b)

########## Last AGAA block in the pattern ##########
df.glm <- overall.df.v2 %>% select(lastAGAA, Affection.Status)
df.glm$Affection.Status = factor(df.glm$Affection.Status, levels = c("U", "A"))
fit <- glm(Affection.Status ~ lastAGAA, data=df.glm, family=binomial())
summary(fit)

########## the number of GGAA second block and the location of second AGAA ##########
overall.df.v2 <- overall.df.v2 %>%  
  mutate(secondGGAA20 = case_when(GGAA_block2 == 5  ~ "yes",
                                  GGAA_block2 != 5  ~ "no"))

df.glm <- overall.df.v2 %>% select(secondGGAA20, Affection.Status)
df.glm$Affection.Status = factor(df.glm$Affection.Status, levels = c("U", "A"))
fit <- glm(Affection.Status ~ secondGGAA20, data=df.glm, family=binomial())
summary(fit) 

##### long mSat and short mSat
df.glm <- overall.df.v2 %>% select(mSat_class, Affection.Status)
df.glm$Affection.Status = factor(df.glm$Affection.Status, levels = c("U", "A"))
df.glm$mSat_class = factor(df.glm$mSat_class, levels = c("short", "long"))
fit <- glm(Affection.Status ~ mSat_class, data=df.glm, family=binomial())
summary(fit) # display results

mod1b <- lrm(Affection.Status ~ mSat_class, data=df.glm)
print(mod1b)

######################################
#    The motif length comparison    #
######################################

EWS.case <- overall.df.v2 %>% filter(Affection.Status== "A")
EWS.control <- overall.df.v2 %>% filter(Affection.Status == "U")
mean(EWS.case$consecutive_GGAA_max)
mean(EWS.control$consecutive_GGAA_max)
mean(EWS.case$GGAAcounts)
mean(EWS.control$GGAAcounts)
mean(EWS.case$mSat_length.x)
mean(EWS.control$mSat_length.x)


EWS.long.mSat <- overall.df.v2 %>% filter(mSat_class== "long")
EWS.short.mSat <- overall.df.v2 %>% filter(mSat_class == "short")
mean(EWS.long.mSat$consecutive_GGAA_max)
mean(EWS.short.mSat$consecutive_GGAA_max)

###################################
#    EWS risk with three SNPs    #
##################################
### logistic regression rs17142617 (chr6:6838025)
df.glm <- overall.df.v3 %>% select(`6838025`, Affection.Status)
df.glm$Affection.Status = factor(df.glm$Affection.Status, levels = c("U", "A"))
fit <- glm(Affection.Status ~ `6838025`, data=df.glm, family=binomial())
summary(fit) 

### logistic regression rs74781311 (chr6:6839193 T>G)
df.glm <- overall.df.v3 %>% select(`6839193`, Affection.Status)
df.glm$Affection.Status = factor(df.glm$Affection.Status, levels = c("U", "A"))
fit <- glm(Affection.Status ~ `6839193`, data=df.glm, family=binomial())
summary(fit) 

### logistic regression rs2876045 (chr6:6839330 T>C)
df.glm <- overall.df.v3 %>% select(`6839330`, Affection.Status)
df.glm$Affection.Status = factor(df.glm$Affection.Status, levels = c("U", "A"))
fit <- glm(Affection.Status ~ `6839330`, data=df.glm, family=binomial())
summary(fit) 


#############################################
#    Correlation with the length of mSat    #
#############################################
### logistic regression rs2876045 (chr6:6839330 T>C)
df.glm <- overall.df.no.SC518934 %>% select(`6839330`, mSat_status)
df.glm$mSat_status = factor(df.glm$mSat_status, levels = c("short", "long"))
fit <- glm(mSat_status ~ `6839330`, data=df.glm, family=binomial())
summary(fit) # display results

### logistic regression rs17142617 (chr6:6838025)
df.glm <- overall.df.no.SC518934 %>% select(`6838025`, mSat_status)
df.glm$mSat_status = factor(df.glm$mSat_status, levels = c("short", "long"))
fit <- glm(mSat_status ~ `6838025`, data=df.glm, family=binomial())
summary(fit)

### logistic regression rs74781311 (chr6:6839193 T>G)
df.glm <- overall.df.no.SC518934 %>% select(`6839193`, mSat_status)
df.glm$mSat_status = factor(df.glm$mSat_status, levels = c("short", "long"))
fit <- glm(mSat_status ~ `6839193`, data=df.glm, family=binomial())
summary(fit) 


################################################################
##   dichotomous or categorical groups based on the counts?   # 
################################################################
# **************** #
#   dichotomize   #
# **************** #
# median 18
overall.df.v4 <- overall.df.v2 %>%  
  mutate(GGAAmedian = case_when(GGAAcounts < 18  ~ 0,
                                GGAAcounts >= 18  ~ 1))
View(overall.df.v4)

### logistic regression
df.glm <- overall.df.v4 %>% select(GGAAmedian, Affection.Status)
df.glm$Affection.Status = factor(df.glm$Affection.Status, levels = c("U", "A"))
df.glm$GGAAmedian = factor(df.glm$GGAAmedian, levels = c("0", "1"))
fit <- glm(Affection.Status ~ GGAAmedian, data=df.glm, family=binomial())
summary(fit) 

### AUC
s=gof(fit, g=10, plotROC=TRUE)
s$auc
### PSeudo R square
mod1b <- lrm(Affection.Status ~ GGAAmedian, data=df.glm)
print(mod1b)

# ********************************* #
#   categorical value - quartile   #
# ********************************* #
overall.df.v4 <- overall.df.v4 %>%  
  mutate(GGAAquartile = case_when(GGAAcounts < 17  ~ 0,
                                  GGAAcounts >= 17 & GGAAcounts <18  ~ 1,
                                  GGAAcounts >= 18 & GGAAcounts <19  ~ 2,
                                  GGAAcounts >= 19  ~ 3))
### logistic regression
df.glm <- overall.df.v4 %>% select(GGAAquartile, Affection.Status)
df.glm$Affection.Status = factor(df.glm$Affection.Status, levels = c("U", "A"))
df.glm$GGAAquartile = factor(df.glm$GGAAquartile, levels = c("0", "1", "2", "3"))
fit <- glm(Affection.Status ~ GGAAquartile, data=df.glm, family=binomial())
summary(fit)

# ***************************** #
#   categorical value - decile  #
# ***************************** #
overall.df.v4 <- overall.df.v4 %>%  
  mutate(GGAAdecile = case_when(GGAAcounts < 15.7  ~ 1,
                                GGAAcounts >= 15.7 & GGAAcounts <17  ~ 2,
                                GGAAcounts >= 17 & GGAAcounts <17  ~ 3,
                                GGAAcounts >= 17 & GGAAcounts <17  ~ 4,
                                GGAAcounts >= 17 & GGAAcounts <18  ~ 5,
                                GGAAcounts >= 18 & GGAAcounts <18.2  ~ 6,
                                GGAAcounts >= 18.2 & GGAAcounts <19  ~ 7,
                                GGAAcounts >= 19 & GGAAcounts <20  ~ 8,
                                GGAAcounts >= 20 & GGAAcounts <23  ~ 9,
                                GGAAcounts >= 23  ~ 10))
overall.df.v4 <- overall.df.v4 %>% select(-GGAAdecile2)

### logistic regression
df.glm <- overall.df.v4 %>% select(GGAAdecile, Affection.Status)
df.glm$Affection.Status = factor(df.glm$Affection.Status, levels = c("U", "A"))
df.glm$GGAAdecile = factor(df.glm$GGAAdecile, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))
fit <- glm(Affection.Status ~ GGAAdecile, data=df.glm, family=binomial())
summary(fit) 
confint(fit)

# ************************************ #
#      decile for the total length     #
# ************************************ #
overall.df.v5 <- overall.df.v4 %>%  
  mutate(mSat_length_decile = case_when(mSat_length.x < 92  ~ 1,
                                        mSat_length.x >= 92 & mSat_length.x <100  ~ 2,
                                        mSat_length.x >= 100 & mSat_length.x <104  ~ 3,
                                        mSat_length.x >= 104 & mSat_length.x <104  ~ 4,
                                        mSat_length.x >= 104 & mSat_length.x <104  ~ 5,
                                        mSat_length.x >= 104 & mSat_length.x <104  ~ 6,
                                        mSat_length.x >= 18.2 & mSat_length.x <19  ~ 7,
                                        mSat_length.x >= 104 & mSat_length.x <108  ~ 8,
                                        mSat_length.x >= 108 & mSat_length.x <112  ~ 9,
                                        mSat_length.x >= 136  ~ 10))

### logistic regression
df.glm <- overall.df.v5 %>% select(mSat_length_decile, Affection.Status)
df.glm$Affection.Status = factor(df.glm$Affection.Status, levels = c("U", "A"))
df.glm$mSat_length_decile = factor(df.glm$mSat_length_decile, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))
fit <- glm(Affection.Status ~ mSat_length_decile, data=df.glm, family=binomial())
summary(fit) 
confint(fit)

ggplot(mSat_length_OR_data,aes(mSatLength.group,estimate,ymin=estimate-std.error,ymax=estimate+std.error))+
  geom_errorbar()+
  geom_pointrange()+
  ggtitle("Odds Ratios for Ewing Sarcoma among GGAA-motif Quartile Groups")+
  xlab("Quartile Group")+
  ylab("Ewing Sarcoma Odds Ratio")+
  theme_classic()

######################################
#     correlation matrix heatmap     #
######################################
library(ComplexHeatmap)
library(corrplot)

overall.df.v2.corr <- overall.df.v2 %>% select(GGAAcounts, AGAAcounts, GGGAcounts, Affection.Status)
overall.df.v2.corr <- overall.df.v2.corr %>% mutate(EWS = case_when(Affection.Status == "A" ~ 1,
                                                                    Affection.Status == "U" ~ 0)) %>%
  select(-Affection.Status)

colnames(res) <-c("GGAA", "AGAA", "GGGA", "EwS")
rownames(res) <-c("GGAA", "AGAA", "GGGA", "EwS")

cormat<-rquery.cormat(res, graphType="heatmap")
corrplot(res, method="number", tl.col = "black", tl.srt = 45,  order = c("original"))

############################################################################
# Please contact Olivia Lee (olivia.lee@nih.gov) if you have any questions about the codes and details.
############################################################################
