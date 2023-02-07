################################################
# R compute Correlation Functional Diversity vs. Bumblebees
# by Simonetta Selva
#
# Created: February 7th, 2022
# Project: Bumblebee 2022
################################################
rm(list=ls())

#load libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)

# set working directory to main repository
input <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/01_DATA"
output <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/03_OUTPUT"

setwd(input)

# load data
BB22.bb.traits <- read_csv("BB22_traits.csv")

# data preparation -----
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
  mutate(site = paste(location, landscape, replicate, sep = ""),
         region = paste(location, landscape, sep = "")) %>%
  dplyr::select(-NrSpecies, -Shannon)
BB22.bb.traits$site <- as.factor(BB22.bb.traits$site)
BB22.bb.traits$ID <- as.factor(BB22.bb.traits$ID)

# B.pascuorum ----
# select only B.pascuorum
bb.traits <- BB22.bb.traits %>%
  filter(bbspecies == "B.pascuorum")%>%
  dplyr::select(ID, site, intertegular_distance, glossa, prementum, proboscis_length, proboscis_ratio, fore_wing_length,
                fore_wing_ratio, corbicula_length, corbicula_ratio)

# import FD of plants
fd.plants <-  read_csv(paste("./FD/FD_package_B.pascuorum_ID.short.csv", sep = "")) %>% 
  rename(ID = ...1,
         nbsp = nbsp.w,
         FRic = FRic.w,
         FEve = FEve.w,
         FDiv = FDiv.w) 

df.all <- merge(bb.traits, fd.plants, by  = "ID", all.x=TRUE)

sitenames <- unique(bb.traits$site)
traits <- colnames(bb.traits[, -c(1,2)])
metrics <- colnames(fd.plants[,-1])
correlation.df <- c()

for (h in sitenames) {
  for (i in traits) {
    for (j in metrics) {
      # compute correlation and store them in a dataframe
      df.all.site <- df.all[df.all$site == h,]
      correlation <- c(h, i, j, round(cor(df.all.site[[i]], df.all.site[[j]], method=c("pearson"), use = "complete.obs"),3))
      correlation.df <- rbind(correlation.df, correlation)
    }
  }
}

correlation.df <- as.data.frame(correlation.df)%>%
  rename(site = V1,
         trait = V2,
         metric = V3,
         correlation = V4)

library(rlist)
library(ggcorrplot)
pcorr.list <- list()

for (i in sitenames) {
  
  correlation.df.site <- filter(correlation.df, site == i)
  correlation.df.site$trait <- as.factor(correlation.df.site$trait)
  matrix.site.species <- split(correlation.df.site, f = correlation.df.site$trait)
  temp <- list.cbind(matrix.site.species)[, seq(4, 36, 4)]
  colnames(temp) <- traits
  rownames(temp) <- metrics
  matrix.site.species <- as.matrix(temp[,match(traits, colnames(temp))])
  class(matrix.site.species) <- "numeric"
  matrix.site.species <- abs(matrix.site.species)
  colnames(matrix.site.species) <- c("int", "glo", "pre", "pro len", "pro rat", "fwi len", "fwi rat", "cor len", "cor rat")

  pcorr <- ggcorrplot(matrix.site.species, hc.order = F,
             outline.col = "white")+ 
    labs(title =i, fill = "", x = "", y = "") +
    theme_classic(base_size = 20) +
    theme(axis.text.x = element_text(angle = 90)) +
    theme(aspect.ratio=1) +
    scale_fill_gradient2(limit = c(0,1), low = "white", high =  "#07575B")
  
  pcorr.list[[i]] <- pcorr
  
}

# arrange them into one file to export
setwd(output)
plot <- ggarrange(pcorr.list[[1]],pcorr.list[[2]],
                  pcorr.list[[3]],pcorr.list[[4]],
                  pcorr.list[[5]],pcorr.list[[6]],
                  pcorr.list[[7]],pcorr.list[[8]],
                  pcorr.list[[9]],pcorr.list[[10]],
                  pcorr.list[[11]],pcorr.list[[12]],
                  pcorr.list[[13]],pcorr.list[[14]],
                  pcorr.list[[15]],pcorr.list[[16]],
                  ncol = 4, nrow = 4,
                  labels = LETTERS[1:9],
                  common.legend = T)
annotate_figure(plot, top = text_grob("B.pascuorum: correlation bumblebee traits vs. plant FDs", 
                                      face = "bold", size = 22))
ggsave("./functional diversity/correlations/FD_B.pascuorum_bb_plants.png", width = 20, height = 20)
setwd(input)

# B.lapidarius ----
# select only B.lapidarius
bb.traits <- BB22.bb.traits %>%
  filter(bbspecies == "B.lapidarius")%>%
  dplyr::select(ID, site, intertegular_distance, glossa, prementum, proboscis_length, proboscis_ratio, fore_wing_length,
                fore_wing_ratio, corbicula_length, corbicula_ratio)

# import FD of plants
fd.plants <-  read_csv(paste("./FD/FD_package_B.lapidarius_ID.short.csv", sep = "")) %>% 
  rename(ID = ...1,
         nbsp = nbsp.w,
         FRic = FRic.w,
         FEve = FEve.w,
         FDiv = FDiv.w) 

df.all <- merge(bb.traits, fd.plants, by  = "ID", all.x=TRUE)

sitenames <- unique(bb.traits$site)
traits <- colnames(bb.traits[, -c(1,2)])
metrics <- colnames(fd.plants[,-1])
correlation.df <- c()

for (h in sitenames) {
  for (i in traits) {
    for (j in metrics) {
      # compute correlation and store them in a dataframe
      df.all.site <- df.all[df.all$site == h,]
      correlation <- c(h, i, j, round(cor(df.all.site[[i]], df.all.site[[j]], method=c("pearson"), use = "complete.obs"),3))
      correlation.df <- rbind(correlation.df, correlation)
    }
  }
}

correlation.df <- as.data.frame(correlation.df)%>%
  rename(site = V1,
         trait = V2,
         metric = V3,
         correlation = V4)

library(rlist)
library(ggcorrplot)
pcorr.list <- list()

for (i in sitenames) {
  
  correlation.df.site <- filter(correlation.df, site == i)
  correlation.df.site$trait <- as.factor(correlation.df.site$trait)
  matrix.site.species <- split(correlation.df.site, f = correlation.df.site$trait)
  temp <- list.cbind(matrix.site.species)[, seq(4, 36, 4)]
  colnames(temp) <- traits
  rownames(temp) <- metrics
  matrix.site.species <- as.matrix(temp[,match(traits, colnames(temp))])
  class(matrix.site.species) <- "numeric"
  matrix.site.species <- abs(matrix.site.species)
  colnames(matrix.site.species) <- c("int", "glo", "pre", "pro len", "pro rat", "fwi len", "fwi rat", "cor len", "cor rat")
  
  pcorr <- ggcorrplot(matrix.site.species, hc.order = F,
                      outline.col = "white")+ 
    labs(title =i, fill = "", x = "", y = "") +
    theme_classic(base_size = 20) +
    theme(axis.text.x = element_text(angle = 90)) +
    theme(aspect.ratio=1) +
    scale_fill_gradient2(limit = c(0,1), low = "white", high =  "#07575B")
  
  pcorr.list[[i]] <- pcorr
  
}

# arrange them into one file to export
setwd(output)
plot <- ggarrange(pcorr.list[[1]],pcorr.list[[2]],
                  pcorr.list[[3]],pcorr.list[[4]],
                  pcorr.list[[5]],pcorr.list[[6]],
                  pcorr.list[[7]],pcorr.list[[8]],
                  pcorr.list[[9]],pcorr.list[[10]],
                  pcorr.list[[11]],pcorr.list[[12]],
                  pcorr.list[[13]],pcorr.list[[14]],
                  pcorr.list[[15]],pcorr.list[[16]],
                  ncol = 4, nrow = 4,
                  labels = LETTERS[1:9],
                  common.legend = T)
annotate_figure(plot, top = text_grob("B.lapidarius: correlation bumblebee traits vs. plant FDs", 
                                      face = "bold", size = 22))
ggsave("./functional diversity/correlations/FD_B.lapidarius_bb_plants.png", width = 20, height = 20)
setwd(input)





