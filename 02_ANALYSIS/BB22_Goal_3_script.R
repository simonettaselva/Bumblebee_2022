################################################
# GOAL 3: Diet Consistency
# by Simonetta Selva
#
# Created: January 18, 2023
# Project: Bumblebee 2022
################################################

# Aim: Examine the diet consistency between and within urban and rural landscapes for both species.


# xxx ----
## preparation ----
# clear environment
rm(list=ls())

# load required library
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rlist)

# set working directory to main repository
input <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/01_DATA"
output <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/03_OUTPUT"

# set working directory
setwd(input) 

# load data
# import species lists on bee pollen (file BB22_sites_plants_GBIF.R)
BB22.full <- read_csv("BB22_full.csv") %>%
  mutate(site = as_factor(site),
         species = as_factor(species),
         region = substr(site, 1, 3)) 

## B.lapidarius ----
BB22.full.species <- BB22.full[BB22.full$bbspecies == "B.lapidarius",]

# convert into binary data frame
BB22.full.table <- t(as.data.frame.matrix(+(table(BB22.full.species[, c("site", "plant.species")]) > 0)))
BB22.full.table <- as.data.frame(BB22.full.table)
BB22.full.table$species=rownames(BB22.full.table)

# create species lists per site
sitenames <- c("ZHUA", "ZHUB", "ZHUC", "ZHRD", "ZHRE", "ZHRF", "BEUA", "BEUB", "BEUC", "BERD", "BSUA", "BSUB", "BSUC", "BSRD", "BSRE", "BSRF")
site.list.bb <- list()

for (i in sitenames) {
  site.list.bb[[i]]  <- unique(BB22.full.species$plant.species[BB22.full.species$site == i])
}

# import species lists form GBIF and InfoFlora (file BB22_sites_plants_GBIF.R)
site.list.occ <- readRDS("sp_list_gbif_infoflora.RData")

# create list with intersection of species lists for all regions
list.intersections <- list()
list.table <- list()
correlation.df.lapi <- data.frame(species = NA, site1 = NA, site2 = NA, Correlation = NA)

for (i in sitenames) {
  for (j in sitenames) {
    # find common species of two sites
    list.intersections[[paste(i,"/",j, sep="")]] <- intersect(site.list.occ[[i]], site.list.occ[[j]])
    
    # filter binary dataframe for the two sites
    list.table[[paste(i,"/",j, sep="")]] <- select(BB22.full.table, i, j)
    
    # compare bumblebee pollen with pollen lists of landscape
    comparison.table <- BB22.full.table[BB22.full.table$species %in% list.intersections[[paste(i,"/",j, sep="")]], c(i,j)]
    
    # compute correlation and store them in a dataframe
    correlation <- c("B.lapidarius", i, j, round(cor(comparison.table[,1], comparison.table[,2], method=c("pearson")),3))
    correlation.df.lapi <- rbind(correlation.df.lapi, correlation)
    
  } # end loop j
} # end loop i


## B.pascuroum ----
BB22.full.species <- BB22.full[BB22.full$bbspecies == "B.pascuorum",]

# convert into binary data frame
BB22.full.table <- t(as.data.frame.matrix(+(table(BB22.full.species[, c("site", "plant.species")]) > 0)))
BB22.full.table <- as.data.frame(BB22.full.table)
BB22.full.table$species=rownames(BB22.full.table)

# create species lists per site
sitenames <- c("ZHUA", "ZHUB", "ZHUC", "ZHRD", "ZHRE", "ZHRF", "BEUA", "BEUB", "BEUC", "BERD", "BSUA", "BSUB", "BSUC", "BSRD", "BSRE", "BSRF")
site.list.bb <- list()

for (i in sitenames) {
  site.list.bb[[i]]  <- unique(BB22.full.species$plant.species[BB22.full.species$site == i])
}

# import species lists form GBIF and InfoFlora (file BB22_sites_plants_GBIF.R)
site.list.occ <- readRDS("sp_list_gbif_infoflora.RData")

# create list with intersection of species lists for all regions
list.intersections <- list()
list.table <- list()
correlation.df.pasc <- data.frame(species = NA, site1 = NA, site2 = NA, Correlation = NA)

for (i in sitenames) {
  for (j in sitenames) {
    # find common species of two sites
    list.intersections[[paste(i,"/",j, sep="")]] <- intersect(site.list.occ[[i]], site.list.occ[[j]])
    
    # filter binary dataframe for the two sites
    list.table[[paste(i,"/",j, sep="")]] <- select(BB22.full.table, i, j)
    
    # compare bumblebee pollen with pollen lists of landscape
    comparison.table <- BB22.full.table[BB22.full.table$species %in% list.intersections[[paste(i,"/",j, sep="")]], c(i,j)]
    
    # compute correlation and store them in a data frame
    correlation <- c("B.pascuorum", i, j, round(cor(comparison.table[,1], comparison.table[,2], method=c("pearson")),3))
    correlation.df.pasc <- rbind(correlation.df.pasc, correlation)
  
    # correlation.matrix.pasc[k, l] <- round(cor(comparison.table[,1], comparison.table[,2], method=c("pearson")),3)
    
  } # end loop j
} # end loop i

correlation.df.lapi <- correlation.df.lapi[-1, ]
correlation.df.pasc <- correlation.df.pasc[-1, ]
correlation.df <- rbind(correlation.df.lapi, correlation.df.pasc)

# bring the data frames into correlation matrix format
# B.lapidarius
correlation.df.lapi$site1 <- as.factor(correlation.df.lapi$site1)
correlation.matrix.lapi <- split(correlation.df.lapi, f = correlation.df.lapi$site1)
temp <- list.cbind(correlation.matrix.lapi)[, seq(4, 67, 4)]
colnames(temp) <- levels(correlation.df.lapi$site1)
rownames(temp) <- c("ZHUA", "ZHUB", "ZHUC", "ZHRD", "ZHRE", "ZHRF", "BEUA", "BEUB", "BEUC", "BERD", "BSUA", "BSUB", "BSUC", "BSRD", "BSRE", "BSRF")
correlation.matrix.lapi <- as.matrix(temp[,match(sitenames, colnames(temp))])
class(correlation.matrix.lapi) <- "numeric"

# B.pascuroum
correlation.df.pasc$site1 <- as.factor(correlation.df.pasc$site1)
correlation.matrix.pasc <- split(correlation.df.pasc, f = correlation.df.pasc$site1)
temp <- list.cbind(correlation.matrix.pasc)[, seq(4, 67, 4)]
colnames(temp) <- levels(correlation.df.pasc$site1)
rownames(temp) <- c("ZHUA", "ZHUB", "ZHUC", "ZHRD", "ZHRE", "ZHRF", "BEUA", "BEUB", "BEUC", "BERD", "BSUA", "BSUB", "BSUC", "BSRD", "BSRE", "BSRF")
correlation.matrix.pasc <- as.matrix(temp[,match(sitenames, colnames(temp))])
class(correlation.matrix.pasc) <- "numeric"

# plot the correlations
library(corrplot)
library(RColorBrewer)
library(ggcorrplot)

par(mfrow = c(1,2))
corrplot(correlation.matrix.lapi, type="upper", order="hclust",
         col=brewer.pal(n=8, name="RdYlBu"))
corrplot(correlation.matrix.pasc, type="upper", order="hclust",
         col=brewer.pal(n=8, name="RdYlBu"))
par(mfrow = c(1,1))

a <- ggcorrplot(correlation.matrix.lapi, hc.order = TRUE, type = "lower", method = "circle",
           outline.col = "white")+ 
  labs(title ='B. lapidarius', fill = "", x = "", y = "") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradient2(limit = c(-0.02,1), low = "white", high =  "#07575B")

b <- ggcorrplot(correlation.matrix.pasc, hc.order = TRUE, type = "lower", method = "circle",
           outline.col = "white")+ 
  labs(title ='B.pascuroum', fill = "", x = "", y = "") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradient2(limit = c(-0.02,1), low = "white", high =  "#07575B")

# arrange them into one file to export
setwd(output)
ggarrange(a, b, ncol = 2, nrow = 1,
          labels = c("A", "B"), common.legend = TRUE, legend = "right")
ggsave("./04_Goal_3/corr_plots_species.png", width = 20, height = 10)
setwd(input)






