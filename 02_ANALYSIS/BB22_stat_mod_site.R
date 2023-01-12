################################################
# Statistical Models script SITE
# by Simonetta Selva
#
# Created: January, 6th, 2023
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
BB22.bb.traits <- read_csv("BB22.bb.traits.csv")

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
  mutate(site = paste(location, landscape, replicate, sep = ""))

#only B.pascuroum 
BB22.bb.traits.sp <- BB22.bb.traits[BB22.bb.traits$bbspecies == "B.pascuorum",]

# import data with spatial information on sites
BB22.sites.meta <- read_csv("BB22_sites_2016.csv")

# import data with funtional diversity of plants per site
BB22.fun.site <- read_csv("./FD/fd_B.pascuorum_site.csv")%>% 
  rename_with(.cols = 1, ~"site")

# compute means of bumblebee traits to compare with functional metrics of plant traits on site level
BB22.bb.traits.site <- BB22.bb.traits.sp %>%
  group_by(site) %>%
  summarise(site = site,
            location = location,
            landscape = landscape,
            intertegular_distance = mean(intertegular_distance),
            glossa = mean(glossa),
            prementum = mean(prementum),
            proboscis_length = mean(proboscis_length),
            proboscis_ratio = mean(proboscis_ratio),
            fore_wing_length = mean(fore_wing_length),
            fore_wing_ratio = mean(fore_wing_ratio),
            corbicula_length = mean(corbicula_length),
            corbicula_ratio = mean(corbicula_ratio)) %>%
  distinct()

# add site coordinates to the trait data frame (in LV95)
BB22.sites <- merge(BB22.bb.traits.site, BB22.sites.meta[, c(1,2,3)], by  = "site", all.x=TRUE) 
BB22.sites <- merge(BB22.sites, BB22.fun.site, by  = "site", all.x=TRUE)

#### look at data ####
# FDis, FRic, FDiv, FEve, FSpe
library(Hmisc)
hist.data.frame(BB22.sites[, -c(1:14)])

# Boxplots for all the variables we want to look at with Wilcoxon test
resp <- c("sp_richn","fdis", "fric", "fdiv", "feve", "fspe")
library(psych)
library(rstatix)
plot_list  <- list()
palette.landscape <- c("#E69F00", "#56B4E9") #create color palette for landscape
for (j in resp) {
  gg.data <- data.frame(landscape=BB22.sites$landscape,value=BB22.sites[[j]])
  w.test <- wilcox_test(gg.data,value~landscape)
  p <- ggplot(gg.data, aes(x=landscape, y = value, fill=landscape)) + 
    geom_boxplot(notch = T) + 
    ylab(j) + xlab("") +
    scale_x_discrete(labels=c('rural', 'urban'))+
    theme_classic(base_size = 20) + guides(alpha = "none") +
    scale_fill_manual(values=palette.landscape, guide = "none") + 
    labs(subtitle = paste("W = ", w.test$statistic, ", p = ", w.test$p, sep=""))
  plot_list[[j]] <- p
  describeBy(gg.data$value, gg.data$landscape)
}

# arrange them into one file to export
setwd(output)
plot <- ggarrange(plot_list[[1]],plot_list[[2]],
                   plot_list[[3]],plot_list[[4]],
                   plot_list[[5]],plot_list[[6]], 
                  ncol = 1, nrow = 6,
                  labels = c("A", "B", "C", "D", "E", "F"))
annotate_figure(plot, top = text_grob("B.pascuorum: species richness and funtional diversity across landscapes", 
                                       face = "bold", size = 22))
# ggsave("./functional diversity/pasc_site/FD_B.pascuorum.png", width = 4, height = 24)
setwd(input)

# plot the relationship of plants traits of one site and bumblebee traits of one site
traits <- colnames(BB22.sites[, 4:12]) # bumblebee traits to look at; prepare for loop
metrics <- colnames(BB22.sites[, c(15, 16, 19, 20, 21, 23)]) # plant FD to look at (see file BB22_compute_FD); prepare for loop

library(nlme)

# perform loop to output plots per relationship summarized per FD
for (i in metrics) {
  x <- 1 # for naming the plots
  for (j in traits) {
    f <- formula(paste(i,"~", j))
    fit <- lme(f, random=~1|landscape, data = BB22.sites)
    assign(paste("a", x, sep=""), # assign the ggplot to plot name
           # define the ggplot
           ggplot(BB22.sites, aes_string(j, i, colour = "landscape")) + 
             geom_point() + 
             theme_classic(base_size = 20) + 
             theme(aspect.ratio=1) + 
             geom_smooth(method="lm", se = FALSE) +
             scale_color_manual(values=palette.landscape, labels=c("rural", "urban")) + 
             stat_cor(aes(color = landscape), size = 5))
    x <- x+1
  } # end loop j
  setwd(output)
  plot4 <- ggarrange(a1,a2,a3,a4,a5,a6,a7,a8, # arrange to plots nicely and export them 
                     ncol = 4, nrow=2, 
                     labels = c(LETTERS[1:8]),   
                     common.legend = TRUE)
  annotate_figure(plot4, top = text_grob(paste("B.pascuorum: comparison of ", i, " and traits across landscapes", sep = ""),
                                         face = "bold", size = 22))
  # ggsave(paste("./functional diversity/pasc_site/FD_pasc_corr_", i, "_BBtraits_landscapes.png", sep = ""), width = 16, height = 8)
  setwd(input)
} # end loop i

# MODELLING
# center and scale all the variables
BB22.sites[, c(4:12, 15, 16, 19, 20, 21, 23)] <- scale(BB22.sites[, c(4:12, 15, 16, 19, 20, 21, 23)],center=TRUE,scale=TRUE)

# correlation analysis
# look at the correlation between the explanatory variables
library(corrplot)
M<-cor(BB22.sites[, 4:12], use = "complete.obs") # subset data; only explanatory variables

# matrix of the p-value of the correlation
p.mat <- cor.mtest(BB22.sites[, 4:12])
head(p.mat)
#correlation plot
corrplot::corrplot(M, type="upper", order="hclust", 
                   p.mat = p.mat$p, sig.level = 0.01, tl.col = "black",
                   col = COL2('RdBu', 10)) # plot correlation with p-values

# !!! a lot are multicollinear. in the model only use: proboscis_ratio, fore_wing_ratio and corbicula_ratio

#  fit GLMM
#load the libraries
library(lme4)
library(car)
library(MuMIn)
library(arm)


# ----------------------------------------------------- Species Richness ---------------------------------------------------

# built an initial full model based on collinearity 
M1.full <- lmer(sp_richn ~ proboscis_ratio + fore_wing_ratio + corbicula_ratio + 
                 (1|landscape),
             data=BB22.sites)
fix.check(M1.full) # looks ok
vif(M1.full) # looks good


# formal test for spatial correlation
sims <- simulateResiduals(M1.full)
BB22.sites$site <- as.factor(BB22.sites$site)
simulationOutput = recalculateResiduals(sims, group = BB22.sites$site)
testSpatialAutocorrelation(sims, x = BB22.sites.meta$LV95_x, y = BB22.sites.meta$LV95_y, plot = FALSE)
# there is no spatial autocorrelation

# update model: add intertegular distance to the model
M1.full.1 <- lmer(sp_richn ~ intertegular_distance + proboscis_ratio + fore_wing_ratio + 
                   corbicula_ratio + (1|landscape),
               data=BB22.sites) # boundary (singular) fit: see help('isSingular')
fix.check(M1.full.1)
summary(M1.full.1)
vif(M1.full.1) # looks good
anova(M1.full.1, M1.full) # with intertegular_distance better fit than M1.full

# update model: remove random effect as its variance is estimated very near zero
M1.full.lm <- lm(sp_richn ~ intertegular_distance + proboscis_ratio + fore_wing_ratio + corbicula_ratio,
                  data=BB22.sites)
summary(M1.full.lm)
fix.check(M1.full.lm)
vif(M1.full.lm)
anova(M1.full.1, M1.full.lm) 

# dredging
# Things to take into account when dredging:
# 1) epends on the models included in the candidate set. You can’t identify a model as being the 
# “best” fit to the data if you didn’t include the model to begin with!
# 2) The parameter estimates and predictions arising from the “best” model or set of best models 
# should be biologically meaningful.

options(na.action = "na.fail") # Required for dredge to run
std.model <- MuMIn::dredge(M1.full.1)
options(na.action = "na.omit") # set back to default
# Get the top best models
top.mod <- get.models(std.model, subset = delta < 6) ## Delta-AICc < 6 (Burnham et al., 2011)
# Model averaging
avg.model <- model.avg(top.mod,revised.var = TRUE)
summary(avg.model)
avg.model$sw

# the only variable having a significant effect in any model is intertegular distance


# ----------------------------------------------------- Functional Richness ---------------------------------------------------
# built an initial full model based on collinearity 
M2.full <- lmer(fric ~ proboscis_ratio + fore_wing_ratio + corbicula_ratio + 
                  (1|landscape),
                data=BB22.sites)
fix.check(M2.full) # looks ok
vif(M2.full) # looks good
summary(M2.full)

# formal test for spatial correlation
sims <- simulateResiduals(M2.full)
BB22.sites$site <- as.factor(BB22.sites$site)
simulationOutput = recalculateResiduals(sims, group = BB22.sites$site)
testSpatialAutocorrelation(sims, x = BB22.sites.meta$LV95_x, y = BB22.sites.meta$LV95_y, plot = FALSE)
# there is no spatial autocorrelation

# add intertegular distance to the model
M2.full.1 <- lmer(fric ~ intertegular_distance + proboscis_ratio + fore_wing_ratio + 
                    corbicula_ratio + (1|landscape),
                  data=BB22.sites) # boundary (singular) fit: see help('isSingular')
fix.check(M2.full.1)
summary(M2.full.1)
vif(M2.full.1) # looks good
anova(M2.full.1, M2.full) # with intertegular_distance better fit than M1.full

# remove random effect as its variance is estimated very near zero
M2.full.lm <- lm(fric ~ intertegular_distance + proboscis_ratio + fore_wing_ratio + corbicula_ratio,
                 data=BB22.sites)
summary(M2.full.lm)
fix.check(M2.full.lm)
vif(M2.full.lm) # looks good
anova(M2.full.1, M2.full.lm) 
# nothing explains the differences in fric

# dredging
options(na.action = "na.fail") # Required for dredge to run
std.model <- MuMIn::dredge(M2.full.1)
options(na.action = "na.omit") # set back to default
# Get the top best models
top.mod <- get.models(std.model, subset = delta < 6) ## Delta-AICc < 6 (Burnham et al., 2011)
# Model averaging
avg.model <- model.avg(top.mod,revised.var = TRUE)
summary(avg.model)
avg.model$sw

# ----------------------------------------------------- Functional Divergence ---------------------------------------------------
# built an initial full model based on collinearity 
M3.full <- lmer(fdiv ~ proboscis_ratio + fore_wing_ratio + corbicula_ratio + (1|landscape),
                data=BB22.sites)
fix.check(M3.full) # looks ok
vif(M3.full) # looks good
summary(M3.full)

# formal test for spatial correlation
sims <- simulateResiduals(M3.full)
BB22.sites$site <- as.factor(BB22.sites$site)
simulationOutput = recalculateResiduals(sims, group = BB22.sites$site)
testSpatialAutocorrelation(sims, x = BB22.sites.meta$LV95_x, y = BB22.sites.meta$LV95_y, plot = FALSE)
# there is no spatial autocorrelation

# add intertegular distance to the model
M3.full.1 <- lmer(fdiv ~ intertegular_distance + proboscis_ratio + fore_wing_ratio + 
                    corbicula_ratio + (1|landscape),
                  data=BB22.sites)
fix.check(M3.full.1)
summary(M3.full.1)
vif(M3.full.1) # looks good
anova(M3.full.1, M3.full) # with intertegular_distance not a better fit than M1.full

# dredging
options(na.action = "na.fail") # Required for dredge to run
std.model <- MuMIn::dredge(M3.full.1)
options(na.action = "na.omit") # set back to default
# Get the top best models
top.mod <- get.models(std.model, subset = delta < 6) ## Delta-AICc < 6 (Burnham et al., 2011)
# Model averaging
avg.model <- model.avg(top.mod,revised.var = TRUE)
summary(avg.model)
avg.model$sw

# ----------------------------------------------------- Functional Evenness ---------------------------------------------------
# built an initial full model based on collinearity 
M4.full <- lmer(feve ~ proboscis_ratio + fore_wing_ratio + corbicula_ratio + (1|landscape),
                data=BB22.sites) #boundary (singular) fit: see help('isSingular')
fix.check(M4.full) # looks ok
vif(M4.full) # looks good
summary(M4.full)

# formal test for spatial correlation
sims <- simulateResiduals(M4.full)
BB22.sites$site <- as.factor(BB22.sites$site)
simulationOutput = recalculateResiduals(sims, group = BB22.sites$site)
testSpatialAutocorrelation(sims, x = BB22.sites.meta$LV95_x, y = BB22.sites.meta$LV95_y, plot = FALSE)
# there is no spatial autocorrelation

# add intertegular distance to the model
M4.full.1 <- lmer(feve ~ intertegular_distance + proboscis_ratio + fore_wing_ratio + 
                    corbicula_ratio + (1|landscape),
                  data=BB22.sites) # boundary (singular) fit: see help('isSingular')
fix.check(M4.full.1)
summary(M4.full.1)
vif(M4.full.1) # looks good
anova(M4.full.1, M4.full) # with intertegular_distance not a better fit than M1.full

# remove random effect as its variance is estimated very near zero
M4.full.lm <- lm(feve ~ intertegular_distance + proboscis_ratio + fore_wing_ratio + corbicula_ratio,
                 data=BB22.sites)
summary(M4.full.lm)
fix.check(M4.full.lm)
vif(M4.full.lm) # looks good
anova(M4.full.1, M4.full.lm) 
# nothing explains the differences in fric

# dredging
options(na.action = "na.fail") # Required for dredge to run
std.model <- MuMIn::dredge(M4.full.1)
options(na.action = "na.omit") # set back to default
# Get the top best models
top.mod <- get.models(std.model, subset = delta < 6) ## Delta-AICc < 6 (Burnham et al., 2011)
# Model averaging
avg.model <- model.avg(top.mod,revised.var = TRUE)
summary(avg.model)
avg.model$sw















#### only B.lapidarius ####

#### Data preparation ####
BB22.bb.traits.sp <- BB22.bb.traits[BB22.bb.traits$bbspecies == "B.lapidarius",]

# import data with spatial information on sites
BB22.sites.meta <- read_csv("BB22_sites_2016.csv")

# import data with funtional diversity of plants per site
BB22.fun.site <- read_csv("./FD/fd_B.lapidarius_site_.csv")%>% 
  rename_with(.cols = 1, ~"site")

# compute means of bumblebee traits to compare with functional metrics of plant traits on site level
BB22.bb.traits.site <- BB22.bb.traits.sp %>%
  group_by(site) %>%
  summarise(site = site,
            location = location,
            landscape = landscape,
            intertegular_distance = mean(intertegular_distance),
            glossa = mean(glossa),
            prementum = mean(prementum),
            proboscis_length = mean(proboscis_length),
            proboscis_ratio = mean(proboscis_ratio),
            fore_wing_length = mean(fore_wing_length),
            fore_wing_ratio = mean(fore_wing_ratio),
            corbicula_length = mean(corbicula_length),
            corbicula_ratio = mean(corbicula_ratio)) %>%
  distinct()

# add site coordinates to the trait data frame (in LV95)
BB22.sites <- merge(BB22.bb.traits.site, BB22.sites.meta[, c(1,2,3)], by  = "site", all.x=TRUE) 
BB22.sites <- merge(BB22.sites, BB22.fun.site, by  = "site", all.x=TRUE)

# look at data
# FDis, FRic, FDiv, FEve, FSpe
library(Hmisc)
hist.data.frame(BB22.sites[, -c(1:14)])

# Boxplots for all the variables we want to look at
resp <- c("sp_richn","fdis", "fric", "fdiv", "feve", "fspe")
library(psych)
plot_list  <- list()
palette.landscape <- c("#E69F00", "#56B4E9") #create color palette for landscape

for (j in resp) {
  gg.data <- data.frame(landscape=BB22.sites$landscape,value=BB22.sites[[j]])
  w.test <- wilcox_test(gg.data,value~landscape)
  p <- ggplot(gg.data, aes(x=landscape, y = value, fill=landscape)) + 
    geom_boxplot(notch = T) + 
    ylab(j) + xlab("") +
    scale_x_discrete(labels=c('rural', 'urban'))+
    theme_classic(base_size = 20) + guides(alpha = "none") +
    scale_fill_manual(values=palette.landscape, guide = "none") + 
    labs(subtitle = paste("W = ", w.test$statistic, ", p = ", w.test$p, sep=""))
  plot_list[[j]] <- p
  describeBy(gg.data$value, gg.data$landscape)
}

# arrange them into one file to export
setwd(output)
plot <- ggarrange(plot_list[[1]],plot_list[[2]],
                  plot_list[[3]],plot_list[[4]],
                  plot_list[[5]],plot_list[[6]], 
                  ncol = 1, nrow = 6,
                  labels = c("A", "B", "C", "D", "E", "F"))
annotate_figure(plot, top = text_grob("B.lapidarius: species richness and funtional diversity across landscapes", 
                                      face = "bold", size = 22))
# ggsave("./functional diversity/lapi_site/FD_B.lapidarius.png", width = 4, height = 24)
setwd(input)

# plot the relationship of plants traits of one site and bumblebee traits of one site
traits <- colnames(BB22.sites[, 4:12]) # bumblebee traits to look at; prepare for loop
metrics <- colnames(BB22.sites[, c(15, 16, 19, 20, 21, 23)]) # plant FD to look at (see file BB22_compute_FD); prepare for loop

# perform loop to output plots per relationship summarized per FD
for (i in metrics) {
  x <- 1 # for naming the plots
  for (j in traits) {
    assign(paste("a", x, sep=""), # assign the ggplot to plot name
           # define the ggplot
           ggplot(BB22.sites, aes_string(j, i, colour = "landscape")) + 
             geom_point() + 
             theme_classic(base_size = 20) + 
             theme(aspect.ratio=1) + 
             geom_smooth(method="lm", se = FALSE)+
             scale_color_manual(values=palette.landscape, labels=c("rural", "urban")) + 
             stat_cor(aes(color = landscape), size = 5))
    x <- x+1
  } # end loop j
  setwd(output)
  plot4 <- ggarrange(a1,a2,a3,a4,a5,a6,a7,a8, # arrange to plots nicely and export them 
                     ncol = 4, nrow=2, 
                     labels = c(LETTERS[1:8]),   
                     common.legend = TRUE)
  annotate_figure(plot4, top = text_grob(paste("B.lapidarius: comparison of ", i, " and traits across landscapes", sep = ""),
                                         face = "bold", size = 22))
  ggsave(paste("./functional diversity/lapi_site/FD_lapi_corr_", i, "_BBtraits_landscapes.png", sep = ""), width = 16, height = 8)
  setwd(input)
} # end loop i


# MODELLING
#load the libraries
library(lme4)
library(car)
library(MuMIn)
library(arm)
library(corrplot)

# center and scale all the variables
BB22.sites[, c(4:12, 15, 16, 19, 20, 21, 23)] <- scale(BB22.sites[, c(4:12, 15, 16, 19, 20, 21, 23)],center=TRUE,scale=TRUE)

# correlation analysis
# look at the correlation between the explanatory variables
M<-cor(BB22.sites[, 4:12], use = "complete.obs") # subset data; only explanatory variables

# matrix of the p-value of the correlation
p.mat <- cor.mtest(BB22.sites[, 4:12])
head(p.mat)
#correlation plot
corrplot::corrplot(M, type="upper", order="hclust", 
                   p.mat = p.mat$p, sig.level = 0.01, tl.col = "black",
                   col = COL2('RdBu', 10)) # plot correlation with p-values

# !!! a lot are multicollinear. in the model only use: proboscis_ratio, fore_wing_ratio and corbicula_ratio

#  fit GLMM


# ----------------------------------------------------- Species Richness ---------------------------------------------------

# built an initial full model based on collinearity 
M5.full <- lmer(sp_richn ~ proboscis_ratio + fore_wing_ratio + corbicula_ratio + 
                  (1|landscape),
                data=BB22.sites)
fix.check(M5.full) # looks ok
vif(M5.full) # looks good
summary(M5.full)

# formal test for spatial correlation
sims <- simulateResiduals(M1.full)
BB22.sites$site <- as.factor(BB22.sites$site)
simulationOutput = recalculateResiduals(sims, group = BB22.sites$site)
testSpatialAutocorrelation(sims, x = BB22.sites.meta$LV95_x, y = BB22.sites.meta$LV95_y, plot = FALSE)
# there is no spatial autocorrelation

# update model: add intertegular distance to the model
M5.full.1 <- lmer(sp_richn ~ intertegular_distance + proboscis_ratio + fore_wing_ratio + 
                    corbicula_ratio + (1|landscape),
                  data=BB22.sites) # boundary (singular) fit: see help('isSingular')
fix.check(M5.full.1)
summary(M5.full.1)
vif(M5.full.1) # looks not really good
anova(M5.full.1, M5.full) # with intertegular_distance better fit than M1.full

# update model: remove random effect as its variance is estimated very near zero
M5.full.lm <- lm(sp_richn ~ intertegular_distance + proboscis_ratio + fore_wing_ratio + corbicula_ratio,
                 data=BB22.sites)
summary(M5.full.lm)
fix.check(M5.full.lm)
vif(M5.full.lm)
anova(M5.full.1, M5.full.lm) 

# dredging
options(na.action = "na.fail") # Required for dredge to run
std.model <- MuMIn::dredge(M5.full.1)
options(na.action = "na.omit") # set back to default
# Get the top best models
top.mod <- get.models(std.model, subset = delta < 6) ## Delta-AICc < 6 (Burnham et al., 2011)
# Model averaging
avg.model <- model.avg(top.mod,revised.var = TRUE)
summary(avg.model)
avg.model$sw


# ----------------------------------------------------- Functional Richness ---------------------------------------------------
# built an initial full model based on collinearity 
M6.full <- lmer(fric ~ proboscis_ratio + fore_wing_ratio + corbicula_ratio + 
                  (1|landscape),
                data=BB22.sites)
fix.check(M6.full) # looks ok
vif(M6.full) # looks good
summary(M6.full)

# formal test for spatial correlation
sims <- simulateResiduals(M6.full)
BB22.sites$site <- as.factor(BB22.sites$site)
simulationOutput = recalculateResiduals(sims, group = BB22.sites$site)
testSpatialAutocorrelation(sims, x = BB22.sites.meta$LV95_x, y = BB22.sites.meta$LV95_y, plot = FALSE)
# there is no spatial autocorrelation

# add intertegular distance to the model
M6.full.1 <- lmer(fric ~ intertegular_distance + proboscis_ratio + fore_wing_ratio + 
                    corbicula_ratio + (1|landscape),
                  data=BB22.sites)
fix.check(M6.full.1)
summary(M6.full.1)
vif(M6.full.1) # looks good
anova(M6.full.1, M6.full) # with M1.full is a better model

# dredging
options(na.action = "na.fail") # Required for dredge to run
std.model <- MuMIn::dredge(M6.full.1)
options(na.action = "na.omit") # set back to default
# Get the top best models
top.mod <- get.models(std.model, subset = delta < 6) ## Delta-AICc < 6 (Burnham et al., 2011)
# Model averaging
avg.model <- model.avg(top.mod,revised.var = TRUE)
summary(avg.model)
avg.model$sw

# ----------------------------------------------------- Functional Divergence ---------------------------------------------------
# built an initial full model based on collinearity 
M6.full <- lmer(fdiv ~ proboscis_ratio + fore_wing_ratio + corbicula_ratio + (1|landscape),
                data=BB22.sites)
fix.check(M6.full) # looks ok
vif(M6.full) # looks good
summary(M6.full)

# formal test for spatial correlation
sims <- simulateResiduals(M6.full)
BB22.sites$site <- as.factor(BB22.sites$site)
simulationOutput = recalculateResiduals(sims, group = BB22.sites$site)
testSpatialAutocorrelation(sims, x = BB22.sites.meta$LV95_x, y = BB22.sites.meta$LV95_y, plot = FALSE)
# there is no spatial autocorrelation

# add intertegular distance to the model
M6.full.1 <- lmer(fdiv ~ intertegular_distance + proboscis_ratio + fore_wing_ratio + 
                    corbicula_ratio + (1|landscape),
                  data=BB22.sites)
fix.check(M6.full.1)
summary(M6.full.1)
vif(M6.full.1) # looks good
anova(M6.full.1, M6.full) # M6.full model is better

# dredging
options(na.action = "na.fail") # Required for dredge to run
std.model <- MuMIn::dredge(M6.full.1)
options(na.action = "na.omit") # set back to default
# Get the top best models
top.mod <- get.models(std.model, subset = delta < 6) ## Delta-AICc < 6 (Burnham et al., 2011)
# Model averaging
avg.model <- model.avg(top.mod,revised.var = TRUE)
summary(avg.model)
avg.model$sw

# ----------------------------------------------------- Functional Evenness ---------------------------------------------------
# built an initial full model based on collinearity 
M8.full <- lmer(feve ~ proboscis_ratio + fore_wing_ratio + corbicula_ratio + (1|landscape),
                data=BB22.sites) #boundary (singular) fit: see help('isSingular')
fix.check(M8.full) # looks ok
vif(M8.full) # looks good
summary(M8.full)

# formal test for spatial correlation
sims <- simulateResiduals(M8.full)
BB22.sites$site <- as.factor(BB22.sites$site)
simulationOutput = recalculateResiduals(sims, group = BB22.sites$site)
testSpatialAutocorrelation(sims, x = BB22.sites.meta$LV95_x, y = BB22.sites.meta$LV95_y, plot = FALSE)
# there is no spatial autocorrelation

# add intertegular distance to the model
M8.full.1 <- lmer(feve ~ intertegular_distance + proboscis_ratio + fore_wing_ratio + 
                    corbicula_ratio + (1|landscape),
                  data=BB22.sites) # boundary (singular) fit: see help('isSingular')
fix.check(M8.full.1)
summary(M8.full.1)
vif(M8.full.1) # looks good
anova(M8.full.1, M8.full) # with intertegular_distance not a better fit than M1.full

# remove random effect as its variance is estimated very near zero
M8.full.lm <- lm(feve ~ intertegular_distance + proboscis_ratio + fore_wing_ratio + corbicula_ratio,
                 data=BB22.sites)
summary(M8.full.lm)
fix.check(M8.full.lm)
vif(M8.full.lm) # looks good
anova(M8.full.1, M8.full.lm) 
# nothing explains the differences in fric

# dredging
options(na.action = "na.fail") # Required for dredge to run
std.model <- MuMIn::dredge(M8.full.1)
options(na.action = "na.omit") # set back to default
# Get the top best models
top.mod <- get.models(std.model, subset = delta < 6) ## Delta-AICc < 6 (Burnham et al., 2011)
# Model averaging
avg.model <- model.avg(top.mod,revised.var = TRUE)
summary(avg.model)
avg.model$sw
