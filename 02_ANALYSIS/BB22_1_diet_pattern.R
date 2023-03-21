################################################
# Foraging patterns of the two bumblebee species
# by Simonetta Selva
#
# Created: May 2nd, 2023
# Project: Bumblebee 2022
################################################

# AIM: Characterize the diet composition and structural (taxonomic, functional, and phylogenetic diversity) 
# and chemical properties of the two bumblebee species in both urban and rural landscapes.

# information: 
# 1) every subsection works in itself
# 2) all the plots not used in the main work are in the script but # are used to not print them


# SPECIES ABBUNDANCES IN POLLEN AND PHYLOGENETIC TREE ----
## preparation ----
# clear environment
rm(list=ls())

# load required library
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(pals)
# library(V.PhyloMaker) old version
library(V.PhyloMaker2)
library(ape)

# set working directory to main repository
input <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/01_DATA"
output <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/03_OUTPUT"

# load data
setwd(input) # set working directory
BB22.full <- read_csv("BB22_full.csv") # import data set with 

# initialize region as column
BB22.full$region <- c() 

# add site and region as columns
for (i in 1:nrow(BB22.full)) {
  BB22.full$site[i] <-paste(BB22.full$location[i], BB22.full$landscape[i], BB22.full$replicate[i], sep="")
  BB22.full$region[i] <- paste(BB22.full$location[i], BB22.full$landscape[i], sep="")
}

## body and corbicula pollen ----

# 1. create a data frame with taxa (species, genus, family)
phylo <- data.frame(species = BB22.full$plant.species, genus = BB22.full$genus, family = BB22.full$family)

# 2. phylogentic tree
# phylogenetic hypotheses under three scenarios based on a backbone phylogeny 
# tree.result <- phylo.maker(phylo, scenarios=c("S1","S2","S3")) # from old version

tree.result <- phylo.maker(phylo) # new version of package
write.tree(tree.result$scenario.3, "BB22_plant_tree.tre") # save for later purpose

# plot the phylogenies with node ages displayed with sceanario 3
tree <- plot.phylo(tree.result$scenario.3, cex = 0.5, main = "Phylogenetic tree of species in pollen"); tree

# get order of species in tree
phylo.order <- data.frame(sps=tree.result$scenario.3$tip.label)
phylo.order$order <- seq(1, length(phylo.order$sps))
colnames(phylo.order) <- c("plant.species", " order")
phylo.order$plant.species <- sub("_", " ", phylo.order$plant.species)

# create new data frame with needed variables
BB22.full.bubble <- BB22.full%>%
  dplyr::group_by(site, plant.species, bbspecies)%>%
  dplyr::summarise(Abundance = sum(Abundance),
                   region = region,
                   landscape = landscape)

# find how many taa were visited per region
BB22.full.taxa <- BB22.full.bubble %>%
  dplyr::group_by(region, bbspecies) %>%
  dplyr::summarise(plant.species = plant.species)%>%
  distinct()

BB22.full.taxa <- BB22.full.taxa%>%
  dplyr::summarise(Nr.taxa = n())

write_csv(BB22.full.taxa, "BB22_NrTaxa_region.csv")


# 3. plotting

# merge two data frames and order along plant species in phylo-tree
BB22.full.bubble_ordered <- merge(x = BB22.full.bubble, y = phylo.order, by.x = "plant.species")%>%
  mutate(plant.species = as_factor(plant.species))
BB22.full.bubble_ordered <- BB22.full.bubble_ordered[order(BB22.full.bubble_ordered$` order`),]
BB22.full.bubble_ordered$plant.species <- factor(BB22.full.bubble_ordered$plant.species, 
                                                 levels = unique(BB22.full.bubble_ordered$plant.species[order(BB22.full.bubble_ordered$` order`)]))

# plot bubble plot with relative abundances along order of species in tree
setwd(output)

# on site level
# palette.site <- kelly(18)[3:18] #create color palette for sites
# ggplot(BB22.full.bubble_ordered, aes(x = site, y =BB22.full.bubble_ordered$plant.species, color = site)) + 
#   geom_point(aes(size = Abundance, fill = site, alpha=0.5)) + 
#   facet_wrap(~bbspecies) +
#   labs(y = "plant species") +
#   theme(axis.text.x = element_text(angle = 90))+
#   theme_classic(base_size = 20) + guides(alpha = "none") +
#   scale_color_manual(values = palette.site, guide = "none")+ #no legend
#   scale_fill_manual(values = palette.site, guide = "none") + #no legend
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# ggsave(paste("./01_Goal 0/Phylo_Bubble_Site.png", sep = ""), width = 16, height = 16)

### FIGURE 1 ####
# region level
palette.landscape <- c("#E69F00", "#56B4E9") #create color palette for landscape
ggplot(BB22.full.bubble_ordered, aes(x = region, y =BB22.full.bubble_ordered$plant.species)) + 
  geom_point(aes(size = Abundance, color = landscape), alpha=0.5) + 
  facet_wrap(~bbspecies) +
  labs(y = "plant species")+    
  theme_classic(base_size = 20) + 
  guides(alpha = "none") +
  scale_color_manual(values = palette.landscape, labels = c("rural", "urban"), name = "Landscape") +
  guides(color = guide_legend(override.aes=list(alpha = 1)))
ggsave(paste("./01_Goal 0/Phylo_Bubble_Region.png", sep = ""), width = 16, height = 16)

# landscape level
# palette.landscape <- c("#E69F00", "#56B4E9") #create color palette for landscape
# ggplot(BB22.full.bubble_ordered, aes(x = landscape, y =BB22.full.bubble_ordered$plant.species, color = landscape)) + 
#   geom_point(aes(size = Abundance, fill = landscape, alpha=0.5)) + 
#   facet_wrap(~bbspecies) +
#   labs(y = "plant species")+ 
#   scale_x_discrete(labels=c('rural', 'urban'))+
#   theme_classic(base_size = 20) + guides(alpha = "none") +
#   scale_color_manual(values = palette.landscape, guide = "none") +
#   scale_fill_manual(values = palette.landscape, guide = "none")
# ggsave(paste("./01_Goal 0/Phylo_Bubble_Landscape.png", sep = ""), width = 16, height = 16)
# setwd(input)


# POLLEN COMPOSITION ----
## preparation ----

# reset environment
rm(list=ls())

#load libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)

# set working directory to main repository
input <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/01_DATA"
output <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/03_OUTPUT"

# load data
setwd(input) # set working directory
BB22.full <- read_csv("BB22_full.csv") # import data set 

# initialize region as column
BB22.full$region <- c() 

# add site and region as columns
for (i in 1:nrow(BB22.full)) {
  BB22.full$site[i] <-paste(BB22.full$location[i], BB22.full$landscape[i], BB22.full$replicate[i], sep="")
  BB22.full$region[i] <- paste(BB22.full$location[i], BB22.full$landscape[i], sep="")
}

## analysis ----
# divide data set into leg and body pollen
BB22.full.leg <- BB22.full[BB22.full$bborgan=="L",]
BB22.full.body <- BB22.full[BB22.full$bborgan=="B",]

## body and corbicula pollen ----
### plant families per site, region and species ----

# 1. find the 30 most abundant families
families.overview <- BB22.full%>%
  dplyr::group_by(family)%>%
  dplyr::summarise(cum.abund = sum(Abundance))
tot <- sum(families.overview$cum.abund)
families.overview <- families.overview %>%
  dplyr::summarise(family = family,
                   cum.abund = cum.abund,
                   rel.abund = cum.abund/tot)
# plot family abundance
ggplot(families.overview, aes(y=reorder(family, -cum.abund), x=cum.abund)) + 
  geom_bar(stat="identity")+ theme_classic()

# use cumulative abundance to select 30 most abundant species
rare.families <- families.overview %>% top_n(nrow(families.overview)-30, -cum.abund)
BB22.full$family.agg <- BB22.full$family
# group all rare families into "other families"
for(h in rare.families$family){
  for(i in 1:nrow(BB22.full)){
    if(BB22.full$family.agg[i] == h){
      BB22.full$family.agg[i] <- "Other families"
    } 
  } # end loop i
} # end loop h

# 2. produce data frame with families' abundances per site with leg pollen
families.site <- BB22.full %>%
  dplyr::group_by(family.agg, site, bbspecies) %>%
  dplyr::summarise(cum.abund = sum(Abundance))
families.site$family.agg <- as.factor(families.site$family.agg)

# 3. plot families per site
# color palette for families
palette.fams=c("#07575B", "#5D535E", "#C4DFE6", "#336B87" , "#FAAF08", "#DFE166", "#1995AD", 
                        "#4897D8", "#80BD9E", "#FA812F", "#66A5AD", "#F34A4A", "#375E97", "#258039",  "#73605B", "#DDBC95",
                        "#A3A599", "#A1D6E2", "#88A550", "#75B1A9", "#D9B44A", "#4F6457", "#ACD0C0", "#0F1B07", 
                        "#F7EFE2", "#D09683", "#F62A00", "#A1BE95", "#20948B", "#9B4F0F", "#CB0000")
                        # filled bar plot per site with abundance of families
# setwd(output)
# ggplot(families.site, aes(fill=family.agg, y=cum.abund, x=site)) + 
#   geom_bar(position="fill", stat="identity")+ 
#   theme_classic(base_size = 20) +
#   facet_wrap(~bbspecies)+ 
#   ggtitle("Plant families per site") +
#   labs(fill='Plant families', x = "sites", y = "relative abundance") +
#   theme(axis.text.x = element_text(angle = 90)) +
#   scale_fill_manual(values=palette.fams, name = "Plant families")
# ggsave(paste("./01_Goal 0/PlantFamilies_per_Site.png", sep = ""), width = 16, height = 8, device = "png")

# 4. produce data frame with families' abundances per region
families.region <- BB22.full %>%
  dplyr::group_by(family.agg, region, bbspecies) %>%
  dplyr::summarise(cum.abund = sum(Abundance))
families.region$family.agg <- as.factor(families.region$family.agg)

### FIGURE S XXX ####
# 5. plot families per region
ggplot(families.region, aes(fill=family.agg, y=cum.abund, x=region)) + 
  geom_bar(position="fill", stat="identity")+ 
  theme_classic(base_size = 20) +
  facet_wrap(~bbspecies)+ 
  ggtitle("Plant families per region") +
  labs(fill='Plant families', x = "regions", y = "relative abundance") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values=palette.fams, name = "Plant families")
ggsave(paste("./01_Goal 0/PlantFamilies_per_Region.png", sep = ""), width = 16, height = 8, device = "png", )
setwd(input) 

### origin status per site, region and species ----
# 1. produce data frame with exotic/native per site
ex.nat.site <- BB22.full %>%
  dplyr::group_by(native_exotic, site, bbspecies) %>%
  dplyr::summarise(abund = sum(binom.abund))
ex.nat.site$native_exotic <- as.factor(ex.nat.site$native_exotic)

# 2. plot native/exotic per site
# create color palette
palette.ex =c("#518A45", "#BDA509")

# setwd(output)
# ggplot(ex.nat.site, aes(fill=native_exotic, y=abund, x=site)) + 
#   geom_bar(position="fill", stat="identity")+ 
#   theme_classic(base_size = 20) +
#   facet_wrap(~bbspecies)+ 
#   ggtitle("Origin status per site LEG") + ylab("Proportion of species in the pollen")+
#   labs(fill='Origin status') +
#   theme(axis.text.x = element_text(angle = 90)) +
#   scale_fill_manual(values=palette.ex, labels=c('Exotic', 'Native'), name = "")
# ggsave(paste("./01_Goal 0/OriginStatus_per_Site.png", sep = ""), width = 16, height = 8)

# 3. produce data frame with exotic/native per region
ex.nat.region <- BB22.full %>%
  dplyr::group_by(native_exotic, region, bbspecies) %>%
  dplyr::summarise(abund = sum(binom.abund))
ex.nat.region$native_exotic <- as.factor(ex.nat.region$native_exotic)

### FIGURE S XXX ####
# 4. plot native/exotic per region    
ggplot(ex.nat.region, aes(fill=native_exotic, y=abund, x=region)) + 
  geom_bar(position="fill", stat="identity")+ 
  theme_classic(base_size = 20) +
  facet_wrap(~bbspecies)+ 
  ggtitle("Origin status per region LEG") + ylab("Proportion of species in the pollen")+
  labs(fill='Origin status') +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values=palette.ex, labels=c('Exotic', 'Native'), name="")
ggsave(paste("./01_Goal 0/OriginStatus_per_Region.png", sep = ""), width = 16, height = 8)
setwd(input)

### growth form per site, region and species ----
# 1. produce data frame with growth form per site
growth.site <- BB22.full %>%
  dplyr::group_by(growth_form_category, site, bbspecies) %>%
  dplyr::summarise(abund = sum(binom.abund))
growth.site$growth_form_category <- as.factor(growth.site$growth_form_category)

# 2. plot growth form per site
# create color palette
palette.growth =c("#375E97", "#518A45", "#FA812F", "#F34A4A", "#5A99AD", "#A3A3A3")

# setwd(output)
# ggplot(growth.site, aes(fill=growth_form_category, y=abund, x=site)) + 
#   geom_bar(position="fill", stat="identity")+ 
#   theme_classic(base_size = 20) +
#   facet_wrap(~bbspecies)+ 
#   ggtitle("Growth Form LEG") + ylab("Proportion of species in the pollen")+ 
#   labs(fill='Growth form') +
#   theme(axis.text.x = element_text(angle = 90)) +
#   scale_fill_manual(values=palette.growth, name="")
# ggsave(paste("./01_Goal 0/GrowthForm_per_Site.png", sep = ""), width = 16, height = 8)

# 3. produce data frame with growth form per region
growth.region <- BB22.full %>%
  dplyr::group_by(growth_form_category, region, bbspecies) %>%
  dplyr::summarise(abund = sum(binom.abund))
growth.region$growth_form_category <- as.factor(growth.region$growth_form_category)

### FIGURE S XXX ####
# 4. plot growth form per region    
ggplot(growth.region, aes(fill=growth_form_category, y=abund, x=region)) + 
  geom_bar(position="fill", stat="identity")+ 
  theme_classic(base_size = 20) +
  facet_wrap(~bbspecies)+ 
  ggtitle("Growth Form LEG") + ylab("Proportion of species in the pollen")+ 
  labs(fill='Growth form') +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values=palette.growth, name="")
ggsave(paste("./01_Goal 0/GrowthForm_per_Region.png", sep = ""), width = 16, height = 8)
setwd(input)

### blossom class per site, region and species ----
# 1. produce data frame with blossom class per site
blossom.site <- BB22.full %>%
  dplyr::group_by(structural_blossom_class, site, bbspecies) %>%
  dplyr::summarise(abund = sum(binom.abund))
blossom.site$structural_blossom_class <- as.factor(blossom.site$structural_blossom_class)

# 2. plot blossom class per site
# create color palette
palette.bloss=c("#375E97", "#80BD9E", "#FA812F", "#758A30", "#07575B", "#D95F56", "#C4DFE6")

# setwd(output)
# ggplot(blossom.site, aes(fill=structural_blossom_class, y=abund, x=site)) + 
#   geom_bar(position="fill", stat="identity")+ 
#   theme_classic(base_size = 20) +
#   facet_wrap(~bbspecies)+ 
#   ggtitle("Blossom Class per site LEG") + ylab("Proportion of species in the pollen")+
#   labs(fill='Blossom Class') +
#   theme(axis.text.x = element_text(angle = 90))+
#   scale_fill_manual(values=palette.bloss, labels=c('Bell Trumpet', 'Brush', "Dish Bowl", "Flag", "Gullet", "Stalk Disk", "Tube"),
#                     name = "")
# ggsave(paste("./01_Goal 0/BlossomClass_per_Site.png", sep = ""), width = 16, height = 8)

# 3. produce data frame with blossom class per region
blossom.region <- BB22.full %>%
  dplyr::group_by(structural_blossom_class, region, bbspecies) %>%
  dplyr::summarise(abund = sum(binom.abund))
blossom.region$structural_blossom_class <- as.factor(blossom.region$structural_blossom_class)

### FIGURE S XXX ####
# 4. plot blossom class per region
ggplot(blossom.region, aes(fill=structural_blossom_class, y=abund, x=region)) + 
  geom_bar(position="fill", stat="identity")+ 
  theme_classic(base_size = 20) +
  facet_wrap(~bbspecies)+ 
  ggtitle("Blossom Class per site LEG") + ylab("Proportion of species in the pollen")+
  labs(fill='Blossom Class') +
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values=palette.bloss, 
                    labels=c('Bell Trumpet', 'Brush', "Dish Bowl",
                             "Flag", "Gullet", "Stalk Disk", "Tube"),
                    name = "")
ggsave(paste("./01_Goal 0/BlossomClass_per_Region.png", sep = ""), width = 16, height = 8)
setwd(input)

# TEST FAMILY, FLOWER, GRWOTH FORM, EXOTIC/NATIVE PROPORTIONS URBAN vs. RURAL ----
## preparation ----
# clear work environment
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

# import data (added row with information on family)
BB22.full.family <- read_csv("BB22.full.family")

## analysis ----

BB.pasc <- BB22.full.family %>% 
  filter(bbspecies == "B.pascuorum")
BB.lapi <- BB22.full.family %>% 
  filter(bbspecies == "B.lapidarius")

# define vectors to loop over
resolution <- c("landscape", "region", "site")
fun.trait <- c("family.agg", "native_exotic", "growth_form_category", "structural_blossom_class")

### B.pascuorum ----
chisq.summary <- NULL
for (i in resolution) {
  for (j in fun.trait) {
    chisq <- chisq.test(BB.pasc[[i]], BB.pasc[[j]]); chisq
    summary <- c(paste(i, j, sep ="~"), round(chisq$statistic, 3), chisq$parameter, chisq$p.value)
    chisq.summary <- rbind(chisq.summary, summary)
  } # end loop j
} # end loop i

chisq.summary <- as.data.frame(chisq.summary)
colnames(chisq.summary) <- c("formula", "X-squared", "df", "p")
rownames(chisq.summary) <- NULL
chisq.summary <- chisq.summary %>%
  mutate(p = format_p(p, stars = TRUE)) %>%
  format_table()

# export the summary of the tests
write.table(chisq.summary, file = "chisq_Bpascuorum.txt", sep = "\t", row.names = TRUE, col.names = NA)

### B. lapidarius ----
chisq.summary.1 <- NULL
for (i in resolution) {
  for (j in fun.trait) {
    chisq <- chisq.test(BB.lapi[[i]], BB.lapi[[j]]); chisq
    summary <- c(paste(i, j, sep ="~"), round(chisq$statistic, 3),chisq$parameter,chisq$p.value)
    chisq.summary.1 <- rbind(chisq.summary.1, summary)
  } # end loop j
} # end loop i

chisq.summary.1 <- as.data.frame(chisq.summary.1)
colnames(chisq.summary.1) <- c("formula", "X-squared", "df", "p")
rownames(chisq.summary.1) <- NULL
chisq.summary.1 <- chisq.summary.1 %>%
  mutate(p = format_p(p, stars = TRUE)) %>%
  format_table()

# export the summary of the tests
write.table(chisq.summary.1, file = "chisq_Blapidarius.txt", sep = "\t", row.names = TRUE, col.names = NA)

# DOCUMENTED PLANT SPECIES PER SITE AND FOUND IN POLLEN ----

## preparation ----
rm(list=ls()) # clear work environment

#load libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)

# set working directory to main repository
input <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/01_DATA"
output <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/03_OUTPUT"

setwd(input)

## analysis ----

# import species lists form GBIF and InfoFlora (file BB22_sites_plants_GBIF.R)
site.list.occ <- readRDS("sp_list_gbif_infoflora.RData")

# load data on pollen (file BB22_sites_plants_GBIF.R)
site.list.bb <- readRDS("site_list_bb.RData")

# combine occurrence data with species found in pollen per site
site.list <- list()
sitenames <- c("ZHUA", "ZHUB", "ZHUC", "ZHRD", "ZHRE", "ZHRF", "BEUA", "BEUB", 
               "BEUC", "BERD", "BSUA", "BSUB", "BSUC", "BSRD", "BSRE", "BSRF")

for (i in sitenames) {
  x <- as.factor(site.list.occ [[i]])
  y <- as.factor(site.list.bb [[i]])
  site.list[[i]] <- unique(c(x,y)) %>% 
    droplevels()
}

# calculate mean of number of species per site
mean.landsacpe <- c()
for (i in sitenames) {
  temp <- c(substring(i, 3, 3), nlevels(site.list[[i]]))
  mean.landsacpe <- rbind(mean.landsacpe, temp)
}

# adapt dataframe to be used with t-test
colnames(mean.landsacpe) <- c("landscape", "species_richness")
mean.landsacpe <- as.data.frame(mean.landsacpe)
mean.landsacpe$species_richness <- as.numeric(mean.landsacpe$species_richness)
str(mean.landsacpe)

# perform t-test and display it in a boxplot
library(rstatix)
w.test <- wilcox_test(mean.landsacpe, species_richness~landscape, paired = F)
palette.landscape <- c("#E69F00", "#56B4E9") #create color palette for landscape
a <- ggplot(mean.landsacpe, aes(x=landscape, y = species_richness,  fill=landscape)) + 
  geom_boxplot(notch = T) + 
  ylab("plant diversity of sites") +
  xlab("") +
  scale_x_discrete(labels=c('rural', 'urban')) +
  theme_classic(base_size = 20) + 
  theme(aspect.ratio=1) + 
  scale_fill_manual(values=palette.landscape, guide = "none") + 
  labs(subtitle = paste("W = ", w.test$statistic, ", p = ", w.test$p, sep="")); a

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
fit <- lm(rich.occ ~ rich.bb, df.site)
summary(fit)

# plot it with information on the model
b <- ggplot(df.site, aes(x = rich.bb , y = rich.occ)) + 
  geom_point(aes(color = landscape), size = 3) + 
  labs(x ="species richness in pollen" , y = "plant diversity of sites") +
  theme_classic(base_size = 20) + 
  theme(aspect.ratio=1) + 
  geom_smooth(method="lm", se = FALSE, col = "black") +
  scale_color_manual(values = palette.landscape) +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
    label.x = 3,
    size = 7); b

# arrange them into one file to export
### FIGURE 3 ----
setwd(output)
ggarrange(a, b, ncol = 2, nrow = 1,
          labels = c("A", "B"))
ggsave("./01_Goal 0/sp_rich_occ_bb.png", width = 20, height = 10)
setwd(input)



# DIET SIMILARITIES ----

rm(list=ls()) # clear work environment

# preparation ----

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
BB22.full <- read_csv("BB22_full.csv") %>%
  mutate(site = as_factor(site),
         species = as_factor(species),
         region = substr(site, 1, 3)) 

#remove information on leg and body from ID
BB22.full$ID.short = as.factor(substring(BB22.full$ID,1, nchar(BB22.full$ID)-1)) 

# sites' names for looping
sitenames <- c("ZHUA", "ZHUB", "ZHUC", "ZHRD", "ZHRE", "ZHRF", "BEUA", "BEUB", "BEUC", "BERD", "BSUA", "BSUB", "BSUC", "BSRD", "BSRE", "BSRF")

## B. pasuorum  ----
# only use B.pascuorum
BB22.full.species <- BB22.full[BB22.full$bbspecies == "B.pascuorum",]

# prepare list element to store correlation matrices for species, genus and family
pasc.ID <- list()
pasc.ID.mean <- list()

### calculate distance ----
for (i in sitenames) {
  # select a site
  BB22.full.species.site <- BB22.full.species[BB22.full.species$site == i,] %>%
    droplevels()
  
  BB22_full.ab <- BB22.full.species.site%>%
    group_by(ID.short, family)%>%
    summarise(abundance = sum(Abundance))%>% 
    distinct() # remove duplicates
  
  BB22_full.ab.new <- c()
  for (h in unique(BB22_full.ab$ID.short)) {
    temp <- BB22_full.ab[BB22_full.ab$ID.short==h,]
    perc <- sum(temp$abundance)
    for (k in 1:nrow(temp)) {
      temp$ab.new[k] <- 100/perc*temp$abundance[k]/100
    }
    BB22_full.ab.new <- rbind(BB22_full.ab.new, temp[, c(1,2,4)])
  }
  
  # convert into matrix
  library(reshape2)
  BB22.full.table <- dcast(BB22_full.ab.new, ID.short ~ family, value.var="ab.new")
  rownames(BB22.full.table) <- BB22.full.table$ID.short
  BB22.full.table[is.na(BB22.full.table)] <- 0
  BB22.full.table <- BB22.full.table[, -1]
  
  # compute measures of distance (or resemblance) and store them in list
  library(vegan)
  dist <- vegdist(BB22.full.table,method="bray")
  pasc.ID[[i]] <- as.numeric(dist[1:1000])
  pasc.ID.mean[[i]] <- mean(dist)
}


### plotting ----
# transform data
pasc.ID.df <- do.call(cbind, pasc.ID)
pasc.ID.df <- melt(pasc.ID.df)[, -1]
colnames(pasc.ID.df) <- c("site", "dist")
pasc.ID.df$landscape <- substr(pasc.ID.df$site, 3, 3)

palette.landscape <- c("#E69F00", "#56B4E9") #create color palette for landscape

#### FIGURE S XXX ----
ggplot(pasc.ID.df, aes(x=site, y = dist, fill=landscape)) + 
  geom_boxplot(notch = T) + 
  xlab("") + ylab("distance") +
  ggtitle("B.pascuorum: distance between individuals per site") +
  theme_classic(base_size = 20) +     
  theme(aspect.ratio=1) + 
  scale_fill_manual(values=palette.landscape, guide = "none") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

setwd(output)
ggsave("./01_Goal 0/distances/pasc_dist_family.png", width = 10, height = 10)
setwd(input)

# B.lapidarius ----
# only use B.B.lapidarius
BB22.full.species <- BB22.full[BB22.full$bbspecies == "B.lapidarius",]

# prepare list element to store correlation matrices for species, genus and family
lapi.ID <- list()
lapi.ID.mean <- list()

### calculate distance ----
for (i in sitenames) {
  # select a site
  BB22.full.species.site <- BB22.full.species[BB22.full.species$site == i,] %>%
    droplevels()
  
  BB22_full.ab <- BB22.full.species.site%>%
    group_by(ID.short, family)%>%
    summarise(abundance = sum(Abundance))%>% 
    distinct() # remove duplicates
  
  BB22_full.ab.new <- c()
  for (h in unique(BB22_full.ab$ID.short)) {
    temp <- BB22_full.ab[BB22_full.ab$ID.short==h,]
    perc <- sum(temp$abundance)
    for (k in 1:nrow(temp)) {
      temp$ab.new[k] <- 100/perc*temp$abundance[k]/100
    }
    BB22_full.ab.new <- rbind(BB22_full.ab.new, temp[, c(1,2,4)])
  }
  
  # convert into matrix
  library(reshape2)
  BB22.full.table <- dcast(BB22_full.ab.new, ID.short ~ family, value.var="ab.new")
  rownames(BB22.full.table) <- BB22.full.table$ID.short
  BB22.full.table[is.na(BB22.full.table)] <- 0
  BB22.full.table <- BB22.full.table[, -1]
  
  # compute measures of distance (or resemblance) and store them in list
  library(vegan)
  dist <- vegdist(BB22.full.table,method="bray")
  lapi.ID[[i]] <- as.numeric(dist[1:1000])
  lapi.ID.mean[[i]] <- mean(dist)
}

### plotting ----
# transform data
lapi.ID.df <- do.call(cbind, lapi.ID)
lapi.ID.df <- melt(lapi.ID.df)[, -1]
colnames(lapi.ID.df) <- c("site", "dist")
lapi.ID.df$landscape <- substr(lapi.ID.df$site, 3, 3)

palette.landscape <- c("#E69F00", "#56B4E9") #create color palette for landscape

#### FIGURE S XXX ----
ggplot(lapi.ID.df, aes(x=site, y = dist, fill=landscape)) + 
  geom_boxplot(notch = T) + 
  xlab("") + ylab("distance") +
  ggtitle("B.lapidarius: distance between individuals per site") +
  theme_classic(base_size = 20) +     
  theme(aspect.ratio=1) + 
  scale_fill_manual(values=palette.landscape, guide = "none") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

setwd(output)
ggsave("./01_Goal 0/distances/lapi_dist_family.png", width = 10, height = 10)
setwd(input)