################################################
# Statistical Models script chemical analysis
# by Simonetta Selva
#
# Created: February, 28th, 2023
# Project: Bumblebee 2022
################################################

rm(list=ls())

# preparation and load data ----

#load libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(glmmTMB)
library(DHARMa)
library(vegan)

# load function
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

BB22.bb.traits.input <- read_csv("BB22_traits.csv")
BB22.chemical.input <- read.csv("BB22_chemical.csv")

# summarize chemical data on site level
BB22.chemical <- BB22.chemical.input %>% 
  group_by(site, bbspecies) %>%
  summarise(location = as.factor(location),
            landscape = as.factor(landscape),
            replicate = as.factor(replicate),
            bbspecies = as.factor(bbspecies),
            site = as.factor(paste(location, landscape, replicate, sep="")),
            AS.mean = mean(AS)) %>%
  mutate(replicate = fct_relevel(replicate,"A", "B", "C", "D", "E", "F")) %>%
  distinct()

# summarize bb traits data on site level
BB22.bb.traits <- BB22.bb.traits.input %>% 
  mutate(landscape = as.factor(landscape))
BB22.bb.traits <- BB22.bb.traits %>%
  mutate(landscape = recode_factor(landscape, rural = "R", urban = "U"),
        site = as.factor(paste(location, landscape, replicate, sep="")))
BB22.bb.traits.mean <- BB22.bb.traits %>%
  group_by(site, bbspecies) %>%
  summarise(location = as.factor(location),
            landscape = as.factor(landscape),
            replicate = as.factor(replicate),
            bbspecies = as.factor(bbspecies),
            glossa.mean = mean(glossa),
            prementum.mean = mean(prementum),
            proboscis_length.mean = mean(proboscis_length),
            proboscis_ratio.mean = mean(proboscis_ratio),
            fore_wing_length.mean = mean(fore_wing_length),
            fore_wing_ratio.mean = mean(fore_wing_ratio),
            corbicula_length.mean = mean(corbicula_length),
            corbicula_ratio.mean = mean(corbicula_ratio))%>%
  distinct()

# combine traits and cwm in new dataframes
BB22.chem.site <- merge(BB22.cwm, BB22.bb.traits.mean[, c(1,9:17)], by  = "ID.short", all.x=TRUE)

# B.PASCUORUM ----------------------------------------------------------------------------------------
#only B.pascuroum 
BB22.bb.sp <- BB22.bb.traits.mean[BB22.bb.traits.mean$bbspecies == "B.pascuorum",]
BB22.chemical.sp <- BB22.chemical[BB22.chemical$bbspecies == "B.pascuorum",]

# combine traits and cwm in new dataframes
BB22.chem.site <- merge(BB22.chemical.sp, BB22.bb.sp[, c(1,6:13)], by  = "site", all.x=TRUE) %>%
  mutate(region = paste(location, landscape, sep=""))

## look at data ----------------------------------------------------------------------------------------

# Boxplots for all the variables we want to look at with Wilcoxon test (landscape level)
library(psych)
library(rstatix)
palette.landscape <- c("#E69F00", "#56B4E9") #create color palette for landscape

w.test <- wilcox_test(BB22.chem.site,AS.mean~landscape)

setwd(output)
ggplot(BB22.chem.site, aes(x=landscape, y = AS.mean, fill=landscape)) + 
  geom_boxplot(notch = T) + 
  xlab("") +
  scale_x_discrete(labels=c('rural', 'urban'))+
  theme_classic(base_size = 20) +     
  theme(aspect.ratio=1) + 
  guides(alpha = "none") +
  scale_fill_manual(values=palette.landscape, guide = "none") + 
  labs(subtitle = paste("p = ", w.test$p, sep=""))
ggsave(paste("./chemical/pasc_site/chemical_pasc_AS_landscapes.png", sep = ""), width = 8, height = 8)
setwd(input)


BB22.chem.site.comp <- BB22.chem.site %>% 
  filter(complete.cases(.)) %>% 
  droplevels()

# plot the relationship of plants traits of one site and bumblebee traits of one site
traits.bb <- colnames(BB22.chem.site.comp[, 7:14]) # bumblebee traits to look at; prepare for loop
plot_list <- list()

for (i in traits.bb) {
  plot_list[[i]] <- 
    ggplot(BB22.chem.site.comp, aes_string(i, "AS.mean", colour = "landscape")) + 
    geom_point() + 
    theme_classic(base_size = 20) + 
    theme(aspect.ratio=1) + 
    geom_smooth(method="lm", se = FALSE) +
    scale_color_manual(values=palette.landscape, labels=c("rural", "urban")) + 
    stat_cor(aes(color = landscape), size = 5)
  
} # end loop i

setwd(output)
plot <- ggarrange(plot_list[[1]],plot_list[[2]],plot_list[[3]],
                  plot_list[[4]],plot_list[[5]],plot_list[[6]],
                  plot_list[[7]],plot_list[[8]],
                  ncol = 4, nrow = 2, 
                  labels = c(LETTERS[1:8]),   
                  common.legend = TRUE)

annotate_figure(plot, top = text_grob(paste("B.pascuorum: comparison of AS and traits across landscapes", sep = ""),
                                       face = "bold", size = 22))
ggsave(paste("./chemical/pasc_site/chemical_pasc_AS_BBtraits_landscapes.png", sep = ""), width = 16, height = 8)
setwd(input)

