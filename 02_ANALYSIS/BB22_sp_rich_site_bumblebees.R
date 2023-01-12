################################################
# Statistical Models Species Richness Site and Bumblebee
# by Simonetta Selva
#
# Created: January, 12th, 2023
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

# import species lists form GBIF and InfoFlora (file BB22_sites_plants_GBIF.R)
site.list.occ <- readRDS("sp_list_gbif_infoflora.RData")

# load data on pollen (file BB22_sites_plants_GBIF.R)
site.list.bb <- readRDS("site_list_bb.RData")

# combine occurance data with species found in pollen per site
site.list <- list()
sitenames <- c("ZHUA", "ZHUB", "ZHUC", "ZHRD", "ZHRE", "ZHRF", "BEUA", "BEUB", "BEUC", "BERD", "BSUA", "BSUB", "BSUC", "BSRD", "BSRE", "BSRF")
for (i in sitenames) {
  x <- as.factor(site.list.occ [[i]])
  y <- as.factor(site.list.bb [[i]])
  site.list[[i]] <- unique(c(x,y))%>% 
    droplevels()
}

# compile species richness per site in dataframe
rich.occ <- c()
for (i in sitenames) {
  x <- nlevels(site.list.occ[[i]])
  rich.occ <- c(rich.occ, x)
}

rich.bb <- c()
for (i in sitenames) {
  x <- nlevels(site.list.bb[[i]])
  rich.bb <- c(rich.bb, x)
}

df.site <- data_frame(site = sitenames,
                      landscape = rep_len(c(rep("U", 3), rep("R", 3)), length(sitenames)),
                      rich.occ = rich.occ,
                      rich.bb = rich.bb)

# plot the relationship of species richness of data bases and species collected by the bumblebees
ggplot(df.site, aes(x = rich.occ, y = rich.bb)) + 
  geom_point()

# test the relationship with a model
fit <- lm(rich.bb~rich.occ, df.site)
summary(fit)

# plot it with information on the model
palette.landscape <- c("#E69F00", "#56B4E9") #create color palette for landscape
ggplot(df.site, aes(x = rich.occ, y = rich.bb)) + 
  geom_point(aes(color = landscape), size = 3) + 
  labs(x = "species occurances data bases", y = "species richness bumblebees") +
  theme_classic(base_size = 20) + 
  theme(aspect.ratio=1) + 
  geom_smooth(method="lm", se = FALSE, col = "black") +
  scale_color_manual(values = palette.landscape) +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
    label.x = 3,
    size = 7)

# save the figure
setwd(output)
ggsave("./GBIF and InfoFlora/sp_rich_occ_bb.png", width = 10, height = 10)
setwd(input)



