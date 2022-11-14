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
library(rstatix)
w.test <- wilcox_test(BB.pasc,Shannon~landscape); w.test
qnorm(w.test$p/2) # z score = -9.000839
w.test$p # p value = 2.24e-19

ggplot(BB.pasc, aes(x=landscape, y=Shannon)) + 
  geom_boxplot(notch = T)+ theme_bw() + labs(subtitle = get_test_label(w.test, detailed = TRUE))


# correlation analysis

pairs.panels(BB.pasc[,c(3,4,9:17)], 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
             ) 
library(ellipse)
plotcorr(cor(BB.pasc[,c(3,4,9:17)]))

# glossa, prementum and fore wings lenght highly correlated to proboscis lenght --> removed from intial model

# PCA for finding correlations
library(ggbiplot)
pca <- prcomp(BB.pasc[,c(3,4,9:17)], center = TRUE, scale. = TRUE)
summary(pca)
ggbiplot(pca, groups=BB.pasc$location, ellipse = TRUE) + #PC1 and PC2
  ggtitle("PCA")+
  theme_minimal()+ theme(aspect.ratio=1)
#  first 2 axis explain 0.6684 of variance
ggbiplot(pca, groups=BB.pasc$site, ellipse = TRUE) + #PC1 and PC2
  ggtitle("PCA")+
  theme_minimal()+ theme(aspect.ratio=1)

#  fit GLMM
#load the libraries
library(lme4)
library(nlme)
library(arm)


#first a random intercept model
mod_lme1<-lmer(log(Shannon+1)~intertegular_distance + proboscis_length + proboscis_ratio 
              + fore_wing_ratio + corbicula_length + corbicula_ratio + (1|landscape),
              data=BB.pasc) # add constant???
summary(mod_lme1)

plot(mod_lme1)
qqnorm(residuals(mod_lme1)) # Short-Tailed but ok??
qqline(residuals(mod_lme1))
hist(log(BB.pasc$Shannon+1))

library(car)
vif(mod_lme1) #cut-off of five (???) to check for collinearity among our explanatory variables

mod_lme1.1<-lmer(log(Shannon+1)~intertegular_distance + proboscis_length +
               fore_wing_ratio  + corbicula_length + corbicula_ratio + (1|landscape),
               data=BB.pasc) # add constant???
summary(mod_lme1.1)
vif(mod_lme1.1)

mod_lme1.2<-lmer(log(Shannon+1)~intertegular_distance + proboscis_length
                 + fore_wing_ratio  + corbicula_length + (1|landscape),
                 data=BB.pasc) # add constant???
summary(mod_lme1.2)
vif(mod_lme1.2) #looks ok



# Variable selection: Stepwise regression model
library(lmerTest)

step.model <- drop1(mod_lme1, test="Chisq")
summary(step.model)

#then a random slope plus intercept model
mod_lme2<-lme(log(Shannon+1)~prementum,data=BB.pasc,random=~prementum|site)
summary(mod_lme2)

BB.pasc$Shannon + 1



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


