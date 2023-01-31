################################################
# Statistical Models script INDIVIDUALS CWM
# by Simonetta Selva
#
# Created: January, 31th, 2023
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

BB22.bb.traits <- read_csv("BB22_traits.csv")
BB22_full <- read.csv("BB22_full.csv")

# compute community weighted means
BB22_full$ID.short = as.factor(substring(BB22_full$ID,1, nchar(BB22_full$ID)-1)) #remove information on leg and body from ID
BB22.cwm <- BB22_full %>% 
  group_by(ID.short) %>%
  summarise(location = as.factor(location),
            landscape = as.factor(landscape),
            replicate = as.factor(replicate),
            bbspecies = as.factor(bbspecies),
            site = as.factor(paste(location, landscape, replicate, sep="")),
            Shannon = diversity(Abundance),
            NrSpecies=n_distinct(species),
            Flowering_duration_cwm = weighted.mean(Flowering_months_duration, Abundance, na.rm=T),
            Flowering_start_cwm = weighted.mean(start_flowering, Abundance, na.rm=T),
            growth_form_cwm = weighted.mean(growth_form_numeric, Abundance, na.rm=T),
            structural_blossom_cwm = weighted.mean(structural_blossom_numeric, Abundance, na.rm=T),
            sugar_concentration_cwm = weighted.mean(sugar.concentration, Abundance, na.rm=T),
            plant_height_cwm = weighted.mean(plant_height_m, Abundance, na.rm=T)
  ) %>%
  mutate(replicate = fct_relevel(replicate,"A", "B", "C", "D", "E", "F"),
         landscape = fct_relevel(landscape, "R", "U")) %>%
  distinct()


# combine trait and cwm in new dataframes
BB22.bb.traits <- BB22.bb.traits %>% 
  rename(ID.short = ID)
BB22.ID.cwm <- merge(BB22.cwm, BB22.bb.traits[, c(1,9:17)], by  = "ID.short", all.x=TRUE)

## B.PASCUORUM ---- 
#only B.pascuroum 
BB22.ID.cwm.sp <- BB22.ID.cwm[BB22.ID.cwm$bbspecies == "B.pascuorum",]


## look at data ----------------------------------------------------------------------------------------

# Boxplots for all the variables we want to look at with Wilcoxon test
traits.pl <- colnames(BB22.ID.cwm.sp[, 7:14])
library(psych)
library(rstatix)
plot_list  <- list()
palette.landscape <- c("#E69F00", "#56B4E9") #create color palette for landscape

for (i in traits.pl) {
  gg.data <- data.frame(landscape=BB22.ID.cwm.sp$landscape,value=BB22.ID.cwm.sp[[i]])
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
  if (i == "Shannon") {
    p <- p + ylim(0, 3) + ylab("shannon diversity index")
  } else if (i == "NrSpecies") {
    p <- p + ylim(0, 20) + ylab("species richness")
  } else if (i == "Flowering_duration_cwm") {
    p <- p + ylim(0, 7) + ylab("flowering duration")
  } else if (i == "Flowering_start_cwm") {
    p <- p + ylim(0, 8) + ylab("flowering start")
  } else if (i == "growth_form_cwm") {
    p <- p + ylim(0, 3) + ylab("growth form (num)")
  } else if (i == "structural_blossom_cwm") {
    p <- p + ylim(0, 6) + ylab("structural blossom (num)")
  } else if (i == "sugar_concentration_cwm") {
    p <- p + ylim(0, 1500) + ylab("sugar concentration")
  }else {
    p <- p + ylim(0, 7) + ylab("plant height")
  }
  plot_list[[i]] <- p
}

# arrange them into one file to export
setwd(output)
plot <- ggarrange(plot_list[[1]],plot_list[[2]],
                  plot_list[[3]],plot_list[[4]],
                  plot_list[[5]],plot_list[[6]],
                  plot_list[[7]],plot_list[[8]],
                  ncol = 4, nrow = 2,
                  labels = LETTERS[1:8])
annotate_figure(plot, top = text_grob("B.pascuorum: CWM across landscapes", 
                                      face = "bold", size = 22))
ggsave("./community weighted means/pasc_ID/CWM_B.pascuorum_ID.png", width = 24, height = 12)
setwd(input)

# plot the relationship of plants traits of one site and bumblebee traits of one site
traits.bb <- colnames(BB22.ID.cwm.sp[, 15:23]) # bumblebee traits to look at; prepare for loop

library(nlme)

# perform loop to output plots per relationship summarized per FD
for (i in traits.pl) {
  x <- 1 # for naming the plots
  for (j in traits.bb) {
    f <- formula(paste(i,"~", j))
    fit <- lme(f, random=~1|landscape, data = BB22.ID.cwm.sp, na.action=na.omit)
    assign(paste("a", x, sep=""), # assign the ggplot to plot name
           # define the ggplot
           ggplot(BB22.ID.cwm.sp, aes_string(j, i, colour = "landscape")) + 
             geom_point() + 
             theme_classic(base_size = 20) + 
             theme(aspect.ratio=1) + 
             geom_smooth(method="lm", se = FALSE) +
             scale_color_manual(values=palette.landscape, labels=c("rural", "urban")) + 
             stat_cor(aes(color = landscape), size = 5))
    x <- x+1
  } # end loop j
  setwd(output)
  plot4 <- ggarrange(a1,a2,a3,a4,a5,a6,a7,a8,a9, # arrange to plots nicely and export them 
                     ncol = 5, nrow = 2, 
                     labels = c(LETTERS[1:9]),   
                     common.legend = TRUE)
  annotate_figure(plot4, top = text_grob(paste("B.pascuorum: comparison of ", i, " and traits across landscapes", sep = ""),
                                         face = "bold", size = 22))
  ggsave(paste("./community weighted means/pasc_ID/CWM_pasc_corr_", i, "_BBtraits_landscapes.png", sep = ""), width = 16, height = 8)
  setwd(input)
} # end loop i

## modelling ----------------------------------------------------------------------------------------
# center and scale all the variables
BB22.ID.cwm.sp[, c(7:23)] <- scale(BB22.ID.cwm.sp[, c(7:23)],center=TRUE,scale=TRUE)

# correlation analysis
# look at the correlation between the explanatory variables
library(corrplot)
M<-cor(BB22.ID.cwm.sp[, 15:23], use = "complete.obs") # subset data; only explanatory variables

# matrix of the p-value of the correlation
p.mat <- cor.mtest(BB22.ID.cwm.sp[, 15:23])
head(p.mat)

#correlation plot
corrplot::corrplot(M, type="upper", order="hclust", 
                   p.mat = p.mat$p, sig.level = 0.01, tl.col = "black",
                   col = COL2('RdBu', 10)) # plot correlation with p-values

# !!! a lot are multicollinear. in the model only use: proboscis_ratio, fore_wing_ratio and corbicula_ratio









