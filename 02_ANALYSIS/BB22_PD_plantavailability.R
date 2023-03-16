################################################
# Statistical Models script PD INDIVIDUALS
# by Simonetta Selva
#
# Created: March 16, 2023
# Project: Bumblebee 2022
################################################

rm(list=ls())
# preparation and load data ----

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
  site.list[[i]] <- unique(c(x,y)) %>% 
    droplevels()
}

# import data with phylogenetic diversity of plants per ID
BB22.fun.ID.pasc <- read_csv("./PD/PD_B.pascuorum_ID.short.csv")%>% 
  dplyr::select(-1)%>%
  rename_with(.cols = 1, ~"ID")%>%
  rename_with(.cols = 2, ~"sp_richn")%>%
  rename_with(.cols = 3, ~"pvar")%>%
  rename_with(.cols = 4, ~"pric")%>%
  rename_with(.cols = 6, ~"pdiv")%>%
  rename_with(.cols = 5, ~"pclu")
BB22.fun.ID.pasc <- BB22.fun.ID.pasc%>%
  mutate(species = rep("B.pascuorum", length(BB22.fun.ID.pasc$ID)))

BB22.fun.ID.lapi <- read_csv("./PD/PD_B.lapidarius_ID.short.csv")%>% 
  dplyr::select(-1)%>%
  rename_with(.cols = 1, ~"ID")%>%
  rename_with(.cols = 2, ~"sp_richn")%>%
  rename_with(.cols = 3, ~"pvar")%>%
  rename_with(.cols = 4, ~"pric")%>%
  rename_with(.cols = 6, ~"pdiv")%>%
  rename_with(.cols = 5, ~"pclu")
BB22.fun.ID.lapi <- BB22.fun.ID.lapi%>%
  mutate(species = rep("B.lapidarius", length(BB22.fun.ID.lapi$ID)))

# combine both species
BB22.fun.ID <- rbind(BB22.fun.ID.pasc[-c(257, 258, 316, 317, 413, 414, 469, 470),], 
                     BB22.fun.ID.lapi[-c(160, 161),]) %>% 
  mutate(site = substr(ID, nchar(ID)-5, nchar(ID)-2))
# BB22.fun.ID$site <- as.factor(BB22.fun.ID$site)
# levels(BB22.fun.ID$site)
# # 
# # sitenames <- c("ZHUA", "ZHUB", "ZHUC", "ZHRD", "ZHRE", "ZHRF", "BEUA", "BEUB", 
# #                "BEUC", "BERD", "BSUA", "BSUB", "BSUC", "BSRD", "BSRE", "BSRF")
# # BB22.fun.ID <- subset(BB22.fun.ID, BB22.fun.ID$site %in% sitenames)

BB22.fun.ID.site <- BB22.fun.ID %>%
  group_by(site, species) %>%
  summarise(pvar = mean(pvar, na.rm = TRUE),
            pric = mean(pric, na.rm = TRUE),
            pdiv = mean(pdiv, na.rm = TRUE),
            pclu = mean(pclu, na.rm = TRUE))

rich.occ <- c()
for (i in sitenames) {
  x <- nlevels(site.list.occ[[i]])
  rich.occ <- c(rich.occ, x)
}
# 
# pric.list <- list()
# for (i in sitenames) {
#   pric.list[[i]] <- BB22.fun.ID.site$pric[BB22.fun.ID.site$site == i]
# }

plant.phy.site <- data.frame(plant.ric = rich.occ,
                             species = as.factor(BB22.fun.ID.site$species),
                             pvar = BB22.fun.ID.site$pvar,
                             pric = BB22.fun.ID.site$pric,
                             pdiv = BB22.fun.ID.site$pdiv,
                             pclu = BB22.fun.ID.site$pclu)

# plot it with information on the model
metrics <- c("pvar", "pric", "pdiv", "pclu")
plot.list <- list()
j <- "pvar"
for (j in metrics) {
  f <- formula(paste("plant.ric ~", j))
  fit <- lme(f, random = ~1|species, data = plant.phy.site, na.action=na.omit)
  plot.list[[j]] <-
    ggplot(plant.phy.site, aes_string(j, "plant.ric", colour = "species")) + 
    geom_point() + 
    theme_classic(base_size = 20) + 
    theme(aspect.ratio=1) + 
    ylim(0, 2500) +
    ylab("plant richness in landscape") +
    geom_smooth(method="lm", se = FALSE) +
    scale_color_manual(values=c("#291600","#e0b802")) + 
    stat_cor(aes(color = species), size = 5)
}

setwd(output)
ggarrange(plot.list[[1]],plot.list[[2]],
          plot.list[[3]],plot.list[[4]],
          ncol = 2, nrow = 2,
          labels = c("A", "B", "C", "D"))
# annotate_figure(plot, top = text_grob("species richness and phylogenetic diversity between species", 
#                                       face = "bold", size = 22))
setwd(output)  
ggsave("./phylogenetic diversity/plant richness landscape/PD_plant_richness_species.png", width = 20, height = 20)
setwd(input)  



