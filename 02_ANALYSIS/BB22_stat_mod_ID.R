################################################
# Statistical Models script INDIVIDUALS
# by Simonetta Selva
#
# Created: January, 24th, 2023
# Project: Bumblebee 2022
################################################

rm(list=ls())

#load libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(glmmTMB)
library(DHARMa)

# load function
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

#function to produce model-checking plots for the fixed effects of an lmer model
fix.check <- function(mod){
  par(mfrow = c(1,3))
  plot(fitted(mod),resid(mod),main="Scale-location plot")	#should have no pattern
  abline(h = 0, col="red", lty=2)
  print(anova(lm(fitted(mod)~resid(mod))))	#should be non-significant
  qqnorm(resid(mod), ylab="Residuals")		#should be approximately straight line
  qqline(resid(mod), col="red")
  plot(density(resid(mod)))					#should be roughly normally distributed
  rug(resid(mod))
  par(mfrow = c(1,1))
}


# set working directory to main repository
input <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/01_DATA"
output <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/03_OUTPUT"

# load data
setwd(input)
BB22.bb.traits <- read_csv("BB22_traits.csv")

#### data preparation ####
# replace the names for urban and rural
for (i in 1:nrow(BB22.bb.traits)) {
  if (BB22.bb.traits$landscape[i] == "urban") {
    BB22.bb.traits$landscape[i] = "U"
  } else {
    BB22.bb.traits$landscape[i] = "R"
  }
}

# rename column site to match other dataframes
BB22.bb.traits <- BB22.bb.traits %>%
  mutate(site = paste(location, landscape, replicate, sep = "")) %>%
  select(-NrSpecies, -Shannon)

#only B.pascuroum 
BB22.bb.traits.sp <- BB22.bb.traits[BB22.bb.traits$bbspecies == "B.pascuorum",]

# import data with spatial information on sites
BB22.sites.meta <- read_csv("BB22_sites_2016.csv")

# import data with funtional diversity of plants per site
BB22.fun.ID <- read_csv("./FD/FD_package_B.pascuorum_ID.short.csv")%>% 
  rename_with(.cols = 1, ~"ID")%>%
  rename_with(.cols = 2, ~"sp_richn")%>%
  rename_with(.cols = 3, ~"fric")%>%
  rename_with(.cols = 5, ~"fdiv")%>%
  rename_with(.cols = 4, ~"feve")
  
# filter functional metrics of plant traits to compare with traits per individual
BB22.fun.ID <- subset(BB22.fun.ID, ID %in% BB22.bb.traits.sp$ID)

# add site coordinates to the trait data frame (in LV95)
BB22.ID <- merge(BB22.bb.traits.sp, BB22.fun.ID, by  = "ID", all.x=TRUE)
BB22.ID <- merge(BB22.ID, BB22.sites.meta[, c(1,2,3)], by  = "site", all.x=TRUE) 

#### look at data ####
# FDis, FRic, FDiv, FEve, FSpe
library(Hmisc)
hist.data.frame(BB22.ID[, c(16:19)])

# Boxplots for all the variables we want to look at with Wilcoxon test
resp <- c("sp_richn", "fric", "fdiv", "feve")
library(psych)
library(rstatix)
plot_list  <- list()
palette.landscape <- c("#E69F00", "#56B4E9") #create color palette for landscape

# for (j in resp) {
  j <- "sp_richn"
  gg.data <- data.frame(landscape=BB22.ID$landscape,value=BB22.ID[[j]])
  w.test <- wilcox_test(gg.data,value~landscape)
  p <- ggplot(gg.data, aes(x=landscape, y = value, fill=landscape)) + 
    geom_boxplot(notch = T) + 
    xlab("") +
    scale_x_discrete(labels=c('rural', 'urban'))+
    theme_classic(base_size = 20) +     
    theme(aspect.ratio=1) + 
    guides(alpha = "none") +
    scale_fill_manual(values=palette.landscape, guide = "none") + 
    labs(subtitle = paste("p = ", w.test$p, sep=""))
  if (j == "sp_richn") {
    p <- p + ylim(0, 18) + ylab("species richness")
  } else if (j == "fric") {
    p <- p + ylim(0, 1) + ylab("funtional richness")
  } else if (j == "fdiv") {
    p <- p + ylim(0.6, 0.85) + ylab("funtional divergence")
  } else {
    p <- p + ylim(0.5, 0.75) + ylab("funtional evenness")
  }
  plot_list[[j]] <- p
  # describeBy(gg.data$value, gg.data$landscape)







