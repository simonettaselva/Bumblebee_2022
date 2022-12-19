################################################
# R compute Functional Diversity
# by Simonetta Selva
#
# Created: December 14, 2022
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
BB22_full <- read_csv("BB22_full.csv")

# choose only to numeric data
BB22_full.numeric <- BB22_full %>% 
  summarise(ID = as.factor(ID),
            location = as.factor(location),
            landscape = as.factor(landscape),
            replicate = as.factor(replicate),
            bbspecies = as.factor(bbspecies),
            bborgan = as.factor(bborgan),
            site = as.factor(paste(location, landscape, replicate, sep="")),
            plant.species = plant.species,
            binom.abund = binom.abund,
            Flowering_duration = Flowering_months_duration,
            Flowering_start = start_flowering,
            growth_form_numeric = growth_form_numeric,
            structural_blossom_numeric = structural_blossom_numeric,
            sugar.concentration = sugar.concentration,
            symmetry_numeric = symmetry_numeric,
            plant_height_m = plant_height_m)
BB22_full.numeric$Flowering_duration <- as.numeric(BB22_full.numeric$Flowering_duration)

# for (i in levels(BB22_full.numeric$bbspecies)) { #loop trough bumblebee species
#   for (j in c("ID", "site", "landscape")) {
    
    j <- "site"
    i <- "B.pascuorum"
    

# 1. Impute missing values and reduce data dimensionality using PCA 
require(caret)
require(vegan)
   
  BB22_full.loop <- BB22_full.numeric[BB22_full.numeric$bbspecies == i,]%>%
    filter(plant.species!="Fabaceae sp.")%>% # Remove this uninteresting entry, where no traits are found
    mutate(Flowering_duration = as.numeric(Flowering_duration))%>% 
    droplevels()
  
  # find Individuals with less than 3 plant species in pollen and remove them from data set
  if (j == "ID") {
    BB22.NrSpecies <- BB22_full[BB22_full$bbspecies == i,]%>% 
      group_by(ID) %>%
      summarise(NrSpecies=n_distinct(plant.species))
    filter.list <- BB22.NrSpecies$ID[BB22.NrSpecies$NrSpecies < 3]
  }
  
BB22_full.loop.species <- BB22_full.loop %>% 
  select(plant.species, Flowering_duration, Flowering_start, growth_form_numeric, 
         structural_blossom_numeric, sugar.concentration, symmetry_numeric, plant_height_m) %>% 
  distinct() # remove duplicates

trt.mis.pred <- preProcess(as.data.frame(BB22_full.loop.species[,-c(1)]), "knnImpute")
traits <- predict(trt.mis.pred, BB22_full.loop.species[,-c(1)]); head(trt.mis.pred)
traits <- as.data.frame(traits)
rownames(traits)  <- BB22_full.loop.species$plant.species
# traits$Flowering_duration <- as.numeric(traits$Flowering_duration)


# PCA -> reduce dimensionality
trt.pca <- prcomp(traits, scale. = T, center = T)
cumsum(trt.pca$sdev/sum(trt.pca$sdev))
trt.scaled <- scores(trt.pca)[,1:2] # adjust number of axes for each group

# 2. prepare data rames for mFD-package
# 2.1 bring plants species per site/ID into wide format (for each BB species)
library(reshape2)
wide <- dcast(BB22_full.loop, BB22_full.loop[[j]] ~ plant.species, value.var="binom.abund")
rownames(wide)  <- levels(BB22_full.loop[[j]])

sp.pa <- decostand(wide[,-1], "pa")
# rownames(sp.pa)  = wide[[j]] #re-introduce rownames
sp.pa <- as.matrix(sp.pa) #turn into matrix

# 2.2 Summarize my assemblages
asb_sp_summ <- mFD::asb.sp.summary(asb_sp_w = sp.pa)
asb_sp_occ <- asb_sp_summ$"asb_sp_occ"

# 2.3 reate data frame with information on the variables (numeric, character, factor etc.)
# for details go to "mFD: General Workflow", by Camille Magneville 2022
traits_cat <- data_frame(trait_name = colnames(traits),
                         trait_type = rep("Q", length(colnames(traits))))

# 2.4 Species traits summary:
traits_summ <- mFD::sp.tr.summary(
  tr_cat     = traits_cat,   
  sp_tr      = traits, 
  stop_if_NA = TRUE)

# 3. Computing distances between species based on functional traits
sp_dist <- mFD::funct.dist(
  sp_tr         = traits,
  tr_cat        = traits_cat,
  metric        = "gower",
  scale_euclid  = "scale_center",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = TRUE)
round(sp_dist, 3)     # Output of the function mFD::funct.dist()

# 4. Computing functional spaces & their quality
# 4.1. Compute multdimensional functional spaces and assess their quality
fspaces_quality <- mFD::quality.fspaces(
  sp_dist             = sp_dist,
  maxdim_pcoa         = 10,
  deviation_weighting = "absolute",
  fdist_scaling       = FALSE,
  fdendro             = "average")
round(fspaces_quality$"quality_fspaces", 3)  # Quality metrics of spaces


# 4.2. Illustrating the quality of the selected functional spaces
mFD::quality.fspaces.plot(
  fspaces_quality            = fspaces_quality,
  quality_metric             = "mad",
  fspaces_plot               = c("tree_average", "pcoa_2d", "pcoa_3d", 
                                 "pcoa_4d", "pcoa_5d", "pcoa_6d"),
  name_file                  = NULL,
  range_dist                 = NULL,
  range_dev                  = NULL,
  range_qdev                 = NULL,
  gradient_deviation         = c(neg = "darkblue", nul = "grey80", pos = "darkred"),
  gradient_deviation_quality = c(low = "yellow", high = "red"),
  x_lab                      = "Trait-based distance")

# on the first (top) row, the y-axis represents species functional distances in the multidimensional space. 
# Thus, the closer species are to the 1:1 line, the better distances in the functional space fit trait-based ones.
# on the second row, the y-axis shows the raw deviation of species distances in the functional space 
# compared to trait-based distances. Thus, the raw deviation reflects the distance to the 1:1 line.
# on the third row, the y-axis shows the absolute or squared deviation of the (“scaled”) distance 
# in the functional space. It is the deviation that is taken into account for computing the quality metric.
# read more here: https://cmlmagneville.github.io/mFD/articles/Compute_and_interpret_quality_of_functional_spaces.html

# 5. Test correlation between functional axes and traits
sp_faxes_coord <- fspaces_quality$"details_fspaces"$"sp_pc_coord"
tr_faxes <- mFD::traits.faxes.cor(
  sp_tr          = traits, 
  sp_faxes_coord = sp_faxes_coord[ , c("PC1", "PC2", "PC3", "PC4")], 
  plot           = TRUE)

# Print traits with significant effect:
tr_faxes$"tr_faxes_stat"[which(tr_faxes$"tr_faxes_stat"$"p.value" < 0.05), ]
tr_faxes$"tr_faxes_plot"

# 6. Plot functional space
big_plot <- mFD::funct.space.plot(
  sp_faxes_coord  = sp_faxes_coord[ , c("PC1", "PC2", "PC3", "PC4")],
  faxes           = c("PC1", "PC2", "PC3", "PC4"),
  name_file       = NULL,
  faxes_nm        = NULL,
  range_faxes     = c(NA, NA),
  color_bg        = "grey95",
  color_pool      = "darkgreen",
  fill_pool       = "white",
  shape_pool      = 21,
  size_pool       = 1,
  plot_ch         = TRUE,
  color_ch        = "black",
  fill_ch         = "white",
  alpha_ch        = 0.5,
  plot_vertices   = TRUE,
  color_vert      = "blueviolet",
  fill_vert       = "blueviolet",
  shape_vert      = 23,
  size_vert       = 1,
  plot_sp_nm      = NULL,
  nm_size         = 3,
  nm_color        = "black",
  nm_fontface     = "plain",
  check_input     = TRUE)

# Plot the graph with all pairs of axes:
big_plot$patchwork


# 7. Compute functional diversity indices & plot them
# 7.1. Functional alpha diversity indices in a multidimensional space
sp.pa <- as.matrix(sp.pa)

alpha_fd_indices <-mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord[ , c("PC1", "PC2", "PC3", "PC4")],
  asb_sp_w         = sp.pa,
  ind_vect         = c("fdis", "fmpd", "fnnd", "feve", "fric", "fdiv", "fori", 
                       "fspe", "fide"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

fd_ind_values <- alpha_fd_indices$"functional_diversity_indices"
assign(paste("fd",i,j, sep = "_"), fd_ind_values) #return FD data frame
write.csv(assign(paste("fd",i,j, sep = "_"), fd_ind_values), file = paste("./FD/fd",i,j, ".csv", sep = "_"))

# FDis Functional Dispersion: the biomass weighted deviation of species traits values from the center of the functional space filled by the assemblage 
#   i.e. the biomass-weighted mean distance to the biomass-weighted mean trait values of the assemblage.
# FRic Functional Richness: the proportion of functional space filled by species of the studied assemblage, 
#   i.e. the volume inside the convex-hull shaping species. To compute FRic the number of species must be at least higher than the number of functional axis + 1.
# FDiv Functional Divergence: the proportion of the biomass supported by the species with the most extreme functional traits 
#   i.e. the ones located close to the edge of the convex-hull filled by the assemblage.
# FEve Functional Evenness: the regularity of biomass distribution in the functional space using the Minimum Spanning Tree 
#   linking all species present in the assemblage.
# FSpe Functional Specialization: the biomass weighted mean distance to the mean position of species from the global pool 
#   (present in all assemblages).
# FMPD Functional Mean Pairwise Distance: the mean weighted distance between all species pairs.
# FNND Functional Mean Nearest Neighbour Distance: the weighted distance to the nearest neighbor within the assemblage.
# FIde Functional Identity: the mean traits values for the assemblage. FIde is always computed when FDis is computed.
# FOri Functional Originality: the weighted mean distance to the nearest species from the global species pool.


# 7.2. Functional beta diversity indices based on multidimensional space
beta_fd_indices <- mFD::beta.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord[ , c("PC1", "PC2", "PC3", "PC4")],
  asb_sp_occ       = asb_sp_occ,
  check_input      = TRUE,
  beta_family      = c("Jaccard"),
  details_returned = TRUE)

beta_fd_indices$pairasb_fbd_indices

} #end of loop j
} #end of loop i







##### LEG POLLEN trial ####


BB22_full.loop <- BB22_full.numeric[BB22_full.numeric$bbspecies == i,]%>%
  filter(plant.species!="Fabaceae sp.", # Remove this uninteresting entry, where no traits are found
         bborgan =="L")%>% 
  mutate(Flowering_duration = as.numeric(Flowering_duration))%>% 
  droplevels()

BB22_full.ID <- BB22_full.loop%>% 
  group_by(ID) %>%
  summarise(NrSpecies=n_distinct(plant.species)); length(BB22_full.ID$ID)
filter.list <- BB22_full.ID$ID[BB22_full.ID$NrSpecies < 3]; length(filter.list)

BB22_full.ID.ftr <- BB22_full.ID[ ! BB22_full.ID$ID %in% filter.list, ]; length(BB22_full.ID.ftr$ID)






