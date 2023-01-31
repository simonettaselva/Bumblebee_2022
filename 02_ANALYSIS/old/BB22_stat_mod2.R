################################################
# Statistical Models script 2
# by Simonetta Selva
#
# Created: November 28, 2022
# Project: Bumblebee 2022
################################################


rm(list=ls())

#  !!!!!!!!!!!!!!!!!!!!! BODY POLLEN ONLY !!!!!!!!!!!!!!!!!!!!!


#load libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(glmmTMB)
library(DHARMa)

# set working directory to main repository
input <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/01_DATA"
output <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/03_OUTPUT"

# load data
setwd(input)
BB22.metrics.body <- read_csv("BB22.metrics.csv")%>%
  filter(bborgan == "B")%>%
  mutate(ID = substring(ID,1, nchar(ID)-1))
BB22.bb.traits <- read_csv("BB22.bb.traits.csv")

BB22.metrics.traits <- merge(BB22.metrics.body, BB22.bb.traits[, -c(2:8)], by = "ID")

#rename site variable
for (i in c(1:nrow(BB22.metrics.traits))) {
  BB22.metrics.traits$site[i] <- paste(BB22.metrics.traits$location[i], BB22.metrics.traits$replicate[i], sep = "_")
}
BB22.metrics.traits$site <- as.factor(BB22.metrics.traits$site)

# import data on sites
BB22.sites.meta <- read_csv("BB22_sites_2016.csv")%>%
  mutate(site = as_factor(site))

# add coordinates to the dataframe (in LV95)
BB22.metrics.traits <- merge(BB22.metrics.traits,BB22.sites.meta[, c(1,3,4)], by  = "site", all.x=TRUE) 

BB.pasc <- BB22.metrics.traits %>% 
  filter(bbspecies == "B.pascuorum") %>% na.omit
BB.lapi <- BB22.metrics.traits %>% 
  filter(bbspecies == "B.lapidarius")


#### B.pascuorum ####
# look at data
hist(BB.pasc$Shannon)
hist(BB.pasc$NrSpecies)

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
pairs.panels(BB.pasc[,c(8:24)], 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
) 

library(ellipse)
plotcorr(cor(BB.pasc[,c(8:24)], use = "complete.obs"))

library(corrplot)
M<-cor(BB.pasc[,c(8:24)], use = "complete.obs")
corrplot::corrplot(M, method="circle", type="lower")

# calculate p values of correlations
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

# matrix of the p-value of the correlation
p.mat <- cor.mtest(BB.pasc[,c(8:24)])
head(p.mat[, 1:5])
corrplot::corrplot(M, type="upper", order="hclust", p.mat = p.mat, sig.level = 0.01)


# glossa, prementum and fore wings lenght highly correlated to proboscis lenght --> removed from intial model

# PCA for finding correlations
library(ggbiplot)

pca <- prcomp(BB.pasc[,c(8:24)], center = TRUE, scale. = TRUE)
summary(pca)
ggbiplot(pca, groups=BB.pasc$location, ellipse = TRUE) + #PC1 and PC2
  ggtitle("PCA")+
  theme_minimal()+ theme(aspect.ratio=1)
#  first 2 axis explain 0.6684 of variance
ggbiplot(pca, groups=BB.pasc$site, ellipse = TRUE) + #PC1 and PC2
  ggtitle("PCA")+
  theme_minimal()+ theme(aspect.ratio=1)
ggbiplot(pca, groups=BB.pasc$landscape, ellipse = TRUE) + #PC1 and PC2
  ggtitle("PCA")+
  theme_minimal()+ theme(aspect.ratio=1)

#  fit GLMM
#load the libraries
library(lme4)
library(nlme)
library(arm)

# Question 1
#first a random intercept model
mod_lme1<-lmer(Shannon~intertegular_distance + proboscis_length + proboscis_ratio 
               + fore_wing_ratio + corbicula_length + corbicula_ratio + (1|site),
               data=BB.pasc) 
BB.pasc.coord=unique(BB.pasc[,c(1,25:26)])
summary(mod_lme1)

plot(mod_lme1)
qqnorm(residuals(mod_lme1)) # Short-Tailed but ok??
qqline(residuals(mod_lme1))
hist(residuals(mod_lme1))

# formal test for spatial correlation
sims <- simulateResiduals(mod_lme1)
BB.pasc$site <- as.factor(BB.pasc$site)
res.rec.mod_lme1 = recalculateResiduals(sims, group = BB.pasc.coord$site)
auto.cor=testSpatialAutocorrelation(res.rec.mod_lme1, x = BB.pasc.coord$LV95_x, y = BB.pasc.coord$LV95_y, plot = FALSE)

plot(sims)

library(car)
vif(mod_lme1) #cut-off of five (???) to check for colinearity among our explanatory variables

# update model
mod_lme1.1<-lmer(Shannon~intertegular_distance + proboscis_length +
                   fore_wing_ratio  + corbicula_length + corbicula_ratio + (1|site),
                 data=BB.pasc) 
summary(mod_lme1.1)
plot(mod_lme1.1)
qqnorm(residuals(mod_lme1.1)) 
qqline(residuals(mod_lme1.1))
hist(residuals(mod_lme1.1))

vif(mod_lme1.1)

# update model
mod_lme1.2<-lmer(Shannon~intertegular_distance + proboscis_length
                 + fore_wing_ratio  + corbicula_length + (1|site),
                 data=BB.pasc)
summary(mod_lme1.2)
plot(mod_lme1.2)
qqnorm(residuals(mod_lme1.2))
qqline(residuals(mod_lme1.2))
hist(residuals(mod_lme1.2))

vif(mod_lme1.2) #looks ok



# Variable selection: Stepwise regression model
library(lmerTest)

step.model <- drop1(mod_lme1, test="Chisq")
summary(step.model)

#then a random slope plus intercept model
mod_lme2<-lme(log(Shannon+1)~prementum,data=BB.pasc,random=~prementum|site)
summary(mod_lme2)

BB.pasc$Shannon + 1

















