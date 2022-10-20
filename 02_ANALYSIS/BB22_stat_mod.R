################################################
# Statistical Models
# by Simonetta Selva
#
# Created: October 20, 2022
# Project: Bumblebee 2022
################################################

rm(list=ls())

#  !!!!!!!!!!!!!!!!!!!!! BODY POLLEN ONLY !!!!!!!!!!!!!!!!!!!!!


#load libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)

# set working directory to main repository
input <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/01_DATA"
output <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/03_OUTPUT"

# load data
setwd(input)
BBtot <- read.csv("BBtot.csv")

BB.pasc <- BBtot %>% 
  filter(bbspecies == "B.pascuorum")
BB.lapi <- BBtot %>% 
  filter(bbspecies == "B.lapidarius")


#### B.pascuorum ####
# look at data
hist(BB.pasc$Shannon) #not normally distributed
hist(BB.pasc$NrSpecies) #not normally distributed


library(psych)
describeBy(BB.pasc$Shannon, BB.pasc$landscape)

ggplot(BB.pasc, aes(x=landscape, y=Shannon)) + 
  geom_boxplot(notch = T)+ theme_bw()

# Wilcoxon-Test
# w.test <- wilcox.test(Shannon~landscape, BB.pasc, alternative = "two.sided"); w.test
library(rstatix)
w.test <- wilcox_test(BB.pasc,Shannon~landscape); w.test
qnorm(w.test$p/2) # z score = -9.000839
w.test$p # p value = 2.24e-19

ggplot(BB.pasc, aes(x=landscape, y=Shannon)) + 
  geom_boxplot(notch = T)+ theme_bw() + labs(subtitle = get_test_label(w.test, detailed = TRUE))


# correlation analysis
library(ggbiplot)

pairs.panels(BB.pasc[,c(2,3,8:16)], 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
             ) 

# PCA for finding correlations
pca <- prcomp(BB.pasc[,c(2,3,8:16)], center = TRUE,scale. = TRUE)
summary(pca)
ggbiplot(pca, groups=BB.pasc$location, ellipse = TRUE) + #PC1 and PC2
  ggtitle("PCA")+
  theme_minimal()+ theme(aspect.ratio=1)
#  first 2 axis explain 0.6684 of variance
ggbiplot(pca, groups=BB.pasc$landscape, ellipse = TRUE) + #PC1 and PC2
  ggtitle("PCA")+
  theme_minimal()+ theme(aspect.ratio=1)

#  fit GLMM
#load the libraries
library(lme4)
library(nlme)
library(arm)

#first a random intercept model
mod_lme1<-lme(Shannon~glossa,data=BB.pasc,random=~1|Beach)
mod_lmer1<-lmer(Richness~NAP+(1|Beach),data=data)
#then a random slope plus intercept model
mod_lme2<-lme(Richness~NAP,data=data,random=NAP|Beach)
mod_lmer2<-lmer(Richness~NAP+(NAP|Beach),data=data)
#Poisson model
mod_glmer1<-glmer(Richness~NAP+(1|Beach),data=data,family="poisson")
#nested and crossed random effect??

#### B.lapidarius ####
# look at data
hist(BB.lapi$Shannon) #not normally distributed
hist(BB.lapi$NrSpecies) #not normally distributed


library(psych)
describeBy(BB.lapi$Shannon, BB.lapi$landscape)

ggplot(BB.lapi, aes(x=landscape, y=Shannon)) + 
  geom_boxplot(notch = T)+ theme_bw()

# Wilcoxon-Test
# w.test <- wilcox.test(Shannon~landscape, BB.lapi, alternative = "two.sided"); w.test
library(rstatix)
w.test <- wilcox_test(BB.lapi,Shannon~landscape); w.test
qnorm(w.test$p/2) # z score = -6.909606
w.test$p # p value = 4.86e-12

ggplot(BB.lapi, aes(x=landscape, y=Shannon)) + 
  geom_boxplot(notch = T)+ theme_bw() + labs(subtitle = get_test_label(w.test, detailed = TRUE))


