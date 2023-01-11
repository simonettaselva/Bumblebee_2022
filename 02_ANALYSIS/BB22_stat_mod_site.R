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
#function to produce model-checking plots for the fixed effects of an lmer model
fix.check <- function(mod){
  par(mfrow = c(1,3))
  plot(fitted(mod),resid(mod),main="Scale-location plot")	#should have no pattern
  abline(h = 0, col="red", lty=2)
  print(anova(lm(fitted(mod)~resid(mod))))	#should be non-significant
  qqnorm(resid(mod), ylab="Residuals")		#should be approximately straight line
  qqline(resid(mod), col="red")
  plot(density(resid(mod)))					#should be roughly normally distributed
  rug(resid(mod))}


# set working directory to main repository
input <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/01_DATA"
output <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/03_OUTPUT"

# load data
setwd(input)
BB22.bb.traits <- read_csv("BB22.bb.traits.csv")

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

# look at data
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

# center and scale all the variables
BB22.sites[, c(4:12, 15, 16, 19, 20, 21, 23)] <- scale(BB22.sites[, c(4:12, 15, 16, 19, 20, 21, 23)],center=TRUE,scale=TRUE)

# correlation analysis
# look at the correlation between the explanatory variables
library(corrplot)
M<-cor(BB22.sites[, 4:12], use = "complete.obs") # subset data; only explanatory variables
corrplot::corrplot(M, method="circle", type="lower") # first look

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

# matrix of the p-value of the correlation
p.mat <- cor.mtest(BB22.sites[, 4:12])
head(p.mat)
corrplot::corrplot(M, type="upper", order="hclust", p.mat = p.mat, sig.level = 0.01) # plot correlation with p-values

# !!! a lot are multicollinear. in the model only use: proboscis_ratio, fore_wing_ratio and corbicula_ratio

#  fit GLMM
#load the libraries
library(lme4)


# Wilcoxon-Test for all FD

qnorm(w.test$p/2) # z score = -9.000839
w.test$p # p value = 2.24e-19

ggplot(BB22.sites, aes(x=landscape, y=sp_richn, fill = landscape)) + 
  geom_boxplot(notch = T) + 
  theme_classic(base_size = 20) +              
  scale_fill_manual(values=palette.landscape, labels=c("rural", "urban")) + 
  labs(subtitle = get_test_label(w.test, detailed = TRUE))




# ----------------------------------------------------- Species Richness ---------------------------------------------------









#define formula for the full model
form <- formula(sp_richn ~ intertegular_distance + glossa + prementum + proboscis_length+proboscis_ratio + 
                  fore_wing_length + fore_wing_ratio + corbicula_length + corbicula_ratio) 
M1.Full <- lmer(sp_richn ~ intertegular_distance + glossa + prementum + proboscis_length + proboscis_ratio + 
                  fore_wing_length + fore_wing_ratio + corbicula_length + corbicula_ratio + (1|landscape),
                data = BB22.sites)
fix.check(M1.Full) # check model assumptions

# perform variable selection
### ROUND 1
M1.A <- update(M1.Full, .~. -intertegular_distance)
M1.B <- update(M1.Full, .~. -glossa)
M1.C <- update(M1.Full, .~. -prementum)
M1.D <- update(M1.Full, .~. -proboscis_length)
M1.E <- update(M1.Full, .~. -proboscis_ratio)
M1.F <- update(M1.Full, .~. -fore_wing_length)
M1.G <- update(M1.Full, .~. -fore_wing_ratio)
M1.H <- update(M1.Full, .~. -corbicula_length)
M1.I <- update(M1.Full, .~. -corbicula_ratio)

anova(M1.Full, M1.A)
anova(M1.Full, M1.B)
anova(M1.Full, M1.C)
anova(M1.Full, M1.D)
anova(M1.Full, M1.E)
anova(M1.Full, M1.F)
anova(M1.Full, M1.G)
anova(M1.Full, M1.H)
anova(M1.Full, M1.I)

# variable selection based on AIC or BIC???
# here I'd remove proboscis lenght (makes note really sense thinking ecologically)

# reduced model based on collinearity 
M1.red <- lmer(sp_richn ~ proboscis_ratio + fore_wing_ratio + corbicula_ratio + 
                 (1|landscape),
             data=BB22.sites) 
anova(M1.Full, M1.red) # better fit than full model
fix.check(M1.red) # looks also better than full model

M1.red.1 <- lmer(sp_richn ~ intertegular_distance + proboscis_ratio + fore_wing_ratio + 
                   corbicula_ratio + (1|landscape),
               data=BB22.sites) 
anova(M1.red.1, M1.red) # with intertegular_distance better fit than M1.red
fix.check(M1.red.1)
summary(M1.red.1)

# reduced model based on personal (ecological?) opinion
M1.red.eco <- lmer(sp_richn ~ intertegular_distance + proboscis_length +
                     fore_wing_length + corbicula_length + 
                 (1|landscape),
               data=BB22.sites) 
anova(M1.Full, M1.red.eco) # better fit than full model
fix.check(M1.red.eco) # looks also better than full model

# looking at only tongue length
M1.tongue <- lmer(sp_richn ~ proboscis_length + (1|landscape),
               data=BB22.sites) 
fix.check(M1.tongue) # looks also better than full model
summary(M1.tongue)








#### only B.lapidarius ####
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

