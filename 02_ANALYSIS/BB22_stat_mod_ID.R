################################################
# Statistical Models script INDIVIDUALS
# by Simonetta Selva
#
# Created: January, 24th, 2023
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

## data preparation ####
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

# B.PASCUORUM ---- 
#only B.pascuroum 
BB22.bb.traits.sp <- BB22.bb.traits[BB22.bb.traits$bbspecies == "B.pascuorum",]

# import data with spatial information on sites
BB22.sites.meta <- read_csv("BB22_sites_2016.csv")

# import data with functional diversity of plants per site
BB22.fun.ID <- read_csv("./FD/FD_package_B.pascuorum_ID.short.csv")%>% 
  rename_with(.cols = 1, ~"ID")%>%
  rename_with(.cols = 2, ~"sp_richn")%>%
  rename_with(.cols = 3, ~"fric")%>%
  rename_with(.cols = 5, ~"fdiv")%>%
  rename_with(.cols = 4, ~"feve")
  
# filter functional metrics of plant traits to compare with traits per individual
BB22.fun.ID <- subset(BB22.fun.ID, ID %in% BB22.bb.traits.sp$ID)

# add site coordinates to the trait data frame (in LV95)
BB22.ID <- merge(BB22.bb.traits.sp, BB22.fun.ID, by  = "ID", all.x=TRUE)
BB22.ID <- merge(BB22.ID, BB22.sites.meta[, c(1,2,3)], by  = "site", all.x=TRUE) 

## look at data ----------------------------------------------------------------------------------------
# FDis, FRic, FDiv, FEve, FSpe
library(Hmisc)
# hist.data.frame(BB22.ID[, c(17:20)])

# Boxplots for all the variables we want to look at with Wilcoxon test
resp <- c("sp_richn", "fric", "fdiv", "feve")
library(psych)
library(rstatix)
plot_list  <- list()
palette.landscape <- c("#E69F00", "#56B4E9") #create color palette for landscape

for (j in resp) {
  # j <- "sp_richn"
  gg.data <- data.frame(landscape=BB22.ID$landscape,value=BB22.ID[[j]])
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
  if (j == "sp_richn") {
    p <- p + ylim(0, 18) + ylab("species richness")
  } else if (j == "fric") {
    p <- p + ylim(0, 1) + ylab("funtional richness")
  } else if (j == "fdiv") {
    p <- p + ylim(0.6, 0.85) + ylab("funtional divergence")
  } else {
    p <- p + ylim(0.5, 0.75) + ylab("funtional evenness")
  }
  plot_list[[j]] <- p
  # describeBy(gg.data$value, gg.data$landscape)
}

# arrange them into one file to export
setwd(output)
plot <- ggarrange(plot_list[[1]],plot_list[[2]],
                  plot_list[[3]],plot_list[[4]],
                  ncol = 1, nrow = 6,
                  labels = c("A", "B", "C", "D"))
annotate_figure(plot, top = text_grob("B.pascuorum: species richness and funtional diversity across landscapes", 
                                      face = "bold", size = 22))
# ggsave("./functional diversity/pasc_ID/FD_B.pascuorum_ID.png", width = 6, height = 24)
setwd(input)



# plot the relationship of plants traits of one site and bumblebee traits of one site
traits <- colnames(BB22.ID[, 7:15]) # bumblebee traits to look at; prepare for loop
metrics <- colnames(BB22.ID[, c(17:20)]) # plant FD to look at (see file BB22_compute_FD); prepare for loop

library(nlme)

# perform loop to output plots per relationship summarized per FD
for (i in metrics) {
  x <- 1 # for naming the plots
  for (j in traits) {
    f <- formula(paste(i,"~", j))
    fit <- lme(f, random=~1|landscape, data = BB22.ID, na.action=na.omit)
    assign(paste("a", x, sep=""), # assign the ggplot to plot name
           # define the ggplot
           ggplot(BB22.ID, aes_string(j, i, colour = "landscape")) + 
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
  ggsave(paste("./functional diversity/pasc_ID/FD_pasc_corr_", i, "_BBtraits_landscapes.png", sep = ""), width = 16, height = 8)
  setwd(input)
} # end loop i


## modelling ----------------------------------------------------------------------------------------
# center and scale all the variables
BB22.ID[, c(7:15, 17:20)] <- scale(BB22.ID[, c(7:15, 17:20)],center=TRUE,scale=TRUE)

# correlation analysis
# look at the correlation between the explanatory variables
library(corrplot)
M<-cor(BB22.ID[, 7:15], use = "complete.obs") # subset data; only explanatory variables

# matrix of the p-value of the correlation
p.mat <- cor.mtest(BB22.ID[, 7:15])
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

# create empty list for best models for model testing
models <- list()

### Species Richness ----------------------------------------------------------------------------------------

# built an initial full model based on collinearity 
M1.full <- lmer(sp_richn ~ proboscis_ratio + fore_wing_ratio + corbicula_ratio + 
                  (1|landscape),
                data=BB22.ID)
fix.check(M1.full) # looks ok
vif(M1.full) # looks good


# formal test for spatial correlation
sims <- simulateResiduals(M1.full)
BB22.ID$site <- as.factor(BB22.ID$site)
simulationOutput <- recalculateResiduals(sims, group = BB22.ID$site)
testSpatialAutocorrelation(simulationOutput, x = BB22.sites.meta$LV95_x, y = BB22.sites.meta$LV95_y, plot = FALSE)
# there is no spatial autocorrelation

# update model: add intertegular distance to the model
M1.full.1 <- lmer(sp_richn ~ intertegular_distance + proboscis_ratio + fore_wing_ratio + 
                    corbicula_ratio + (1|landscape),
                  data=BB22.ID)
fix.check(M1.full.1)
summary(M1.full.1)
vif(M1.full.1) # looks good
anova(M1.full.1, M1.full) # with intertegular_distance better fit than M1.full

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

### Functional Richness ----------------------------------------------------------------------------------------
# built an initial full model based on colinearity 
M2.full <- lmer(fric ~ proboscis_ratio + fore_wing_ratio + corbicula_ratio + 
                  (1|landscape),
                data=BB22.ID)
fix.check(M2.full) # looks ok
vif(M2.full) # looks good
summary(M2.full)

# formal test for spatial correlation
sims <- simulateResiduals(M2.full)
BB22.ID$site <- as.factor(BB22.ID$site)
simulationOutput <- recalculateResiduals(sims, group = BB22.ID$site)
testSpatialAutocorrelation(simulationOutput, x = BB22.sites.meta$LV95_x, y = BB22.sites.meta$LV95_y, plot = FALSE)
# there is no spatial autocorrelation

# add intertegular distance to the model
M2.full.1 <- lmer(fric ~ intertegular_distance + proboscis_ratio + fore_wing_ratio + 
                    corbicula_ratio + (1|landscape),
                  data=BB22.ID) # boundary (singular) fit: see help('isSingular')
fix.check(M2.full.1)
summary(M2.full.1)
vif(M2.full.1) # looks good
anova(M2.full.1, M2.full) # with intertegular_distance better fit than M1.full

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

### Functional Divergence ----------------------------------------------------------------------------------------
# built an initial full model based on collinearity 
M3.full <- lmer(fdiv ~ proboscis_ratio + fore_wing_ratio + corbicula_ratio + (1|landscape),
                data=BB22.ID)
fix.check(M3.full) # looks ok
vif(M3.full) # looks good
summary(M3.full)

# formal test for spatial correlation
sims <- simulateResiduals(M3.full)
BB22.ID$site <- as.factor(BB22.ID$site)
simulationOutput <- recalculateResiduals(sims, group = BB22.ID$site)
testSpatialAutocorrelation(simulationOutput, x = BB22.sites.meta$LV95_x, y = BB22.sites.meta$LV95_y, plot = FALSE)
# there is no spatial autocorrelation

# add intertegular distance to the model
M3.full.1 <- lmer(fdiv ~ intertegular_distance + proboscis_ratio + fore_wing_ratio + 
                    corbicula_ratio + (1|landscape),
                  data=BB22.ID)
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

### Functional Evenness ----------------------------------------------------------------------------------------
# built an initial full model based on collinearity 
M4.full <- lmer(feve ~ proboscis_ratio + fore_wing_ratio + corbicula_ratio + (1|landscape),
                data=BB22.ID) #boundary (singular) fit: see help('isSingular')
fix.check(M4.full) # looks ok
vif(M4.full) # looks good
summary(M4.full)

# formal test for spatial correlation
sims <- simulateResiduals(M4.full)
BB22.ID$site <- as.factor(BB22.ID$site)
simulationOutput = recalculateResiduals(sims, group = BB22.ID$site)
testSpatialAutocorrelation(sims, x = BB22.sites.meta$LV95_x, y = BB22.sites.meta$LV95_y, plot = FALSE)
# there is no spatial autocorrelation

# add intertegular distance to the model
M4.full.1 <- lmer(feve ~ intertegular_distance + proboscis_ratio + fore_wing_ratio + 
                    corbicula_ratio + (1|landscape),
                  data=BB22.ID) # boundary (singular) fit: see help('isSingular')
fix.check(M4.full.1)
summary(M4.full.1)
vif(M4.full.1) # looks good
anova(M4.full.1, M4.full) # with intertegular_distance not a better fit than M1.full

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

## Model testing ----------------------------------------------------------------------------------------

### Data splitting  ----------------------------------------------------------------------------------------
library(rsample)
set.seed(123) # for reproducibility
BB22.ID.complete <- BB22.ID %>%
  filter(complete.cases(.))
split <- initial_split(BB22.ID.complete, prop = 0.8)
ddf_train <- training(split)
ddf_test <- testing(split) 

### Model training and predicting ----------------------------------------------------------------------------------------
library(caret) 
library(yardstick)

#### Linear Model ----------------------------------------------------------------------------------------

for (i in metrics) {
  f <- as.formula(paste(i,  "~ intertegular_distance + proboscis_ratio + fore_wing_ratio + corbicula_ratio + landscape", sep = ""))
  train <- train(form = f,
                 data = ddf_train, 
                 method = "lm",
                 trControl = trainControl(method = "cv", number = 10))
  
  ## Get variable importance, and turn into a data frame
  var_imp <- varImp(train, scale=FALSE)$importance
  var_imp <- data.frame(variables=row.names(var_imp), importance=var_imp$Overall)
  
  ## Create a plot of variable importance
  var_imp %>%
    
    ## Sort the data by importance
    arrange(importance) %>%
    
    ## Create a ggplot object for aesthetic
    ggplot(aes(x=reorder(variables, importance), y=importance)) + 
    
    ## Plot the bar graph
    geom_bar(stat='identity') + 
    
    ## Flip the graph to make a horizontal bar plot
    coord_flip() + 
    
    ## Add x-axis label
    xlab('Variables') +
    
    ## Add a title
    labs(title='Linear model variable importance') + 
    
    ## Some layout for the plot
    theme_minimal() + 
    theme(axis.text = element_text(size = 10), 
          axis.title = element_text(size = 15), 
          plot.title = element_text(size = 20), 
    )
  setwd(output)
  ggsave(paste("./functional diversity/pasc_ID/training and predicting/linear model/", i, "_lm_var_imp_pasc.png", sep = ""), width = 8, height = 5)
  setwd(input)
  
  predict_train <- predict(
    ## lm object
    object=train, 
    ## Data to use for predictions
    newdata=ddf_train)
  
  predict_test <- predict(
    ## lm object
    object=train, 
    ## Data to use for predictions
    newdata=ddf_test)
  
  data_metrics_train <- data.frame(truth=ddf_train$feve, pred=predict_train)
  metrics_train <- metrics(data_metrics_train, truth, pred)
  
  data_metrics_test <- data.frame(truth=ddf_test$feve, pred=predict_test)
  metrics_test <- metrics(data_metrics_test, truth, pred)
  
  # plotting prediction
  library(patchwork)
  
  gg_test <- ggplot(data_metrics_test, aes(x=truth, y=pred))+
    geom_point()+
    labs(title = "Density of Data Point for Training Data",
         subtitle = paste("Linear Model: Rsq = ", round(metrics_test$.estimate[2], 3), sep = ""), 
         x = paste("predicted values for ", i, sep = ""),
         y = paste("observed values for ", i, sep = "")) +
    theme(aspect.ratio=1)+
    geom_abline(intercept = 0, slope = 1, color = "#fc5e03")
  
  gg_train <- ggplot(data_metrics_train, aes(x=truth, y=pred))+
    geom_point()+
    labs(title = "Density of Data Point for Training Data",
         subtitle = paste("Linear Model: Rsq = ", round(metrics_train$.estimate[2], 3), sep = ""),
         x = paste("predicted values for ", i, sep = ""),
         y = paste("observed values for ", i, sep = "")) +
    theme(aspect.ratio=1)+
    geom_abline(intercept = 0, slope = 1, color = "#fc5e03")
  
  # arrange them into one file to export
  setwd(output)
  ggarrange(gg_test, gg_train, ncol = 2, nrow = 1,
            labels = c("A", "B"), common.legend = TRUE, legend = "right")
  ggsave(paste("./functional diversity/pasc_ID/training and predicting/linear model/", i, "_lm_prediction_pasc.png", sep = ""), width = 8, height = 5)
  setwd(input)
} # end i loop

#### Random Forest ----------------------------------------------------------------------------------------

for (i in metrics) {
  f <- as.formula(paste(i,  "~ intertegular_distance + proboscis_ratio + fore_wing_ratio + corbicula_ratio + landscape", sep = ""))
  train <- train(form = f,
                 data = ddf_train, 
                 method = "rf",
                 trControl = trainControl(method = "cv", number = 10))
  
  ## Get variable importance, and turn into a data frame
  var_imp <- varImp(train, scale=FALSE)$importance
  var_imp <- data.frame(variables=row.names(var_imp), importance=var_imp$Overall)
  
  ## Create a plot of variable importance
  var_imp %>%
    
    ## Sort the data by importance
    arrange(importance) %>%
    
    ## Create a ggplot object for aesthetic
    ggplot(aes(x=reorder(variables, importance), y=importance)) + 
    
    ## Plot the bar graph
    geom_bar(stat='identity') + 
    
    ## Flip the graph to make a horizontal bar plot
    coord_flip() + 
    
    ## Add x-axis label
    xlab('Variables') +
    
    ## Add a title
    labs(title='Random forest variable importance') + 
    
    ## Some layout for the plot
    theme_minimal() + 
    theme(axis.text = element_text(size = 10), 
          axis.title = element_text(size = 15), 
          plot.title = element_text(size = 20), 
    )
  setwd(output)
  ggsave(paste("./functional diversity/pasc_ID/training and predicting/random forest/", i, "_rf_var_imp_pasc.png", sep = ""), width = 8, height = 8)
  setwd(input)
  
  predict_train <- predict(
    ## lm object
    object=train, 
    ## Data to use for predictions
    newdata=ddf_train)
  
  predict_test <- predict(
    ## lm object
    object=train, 
    ## Data to use for predictions
    newdata=ddf_test)
  
  data_metrics_train <- data.frame(truth=ddf_train$feve, pred=predict_train)
  metrics_train <- metrics(data_metrics_train, truth, pred)
  
  data_metrics_test <- data.frame(truth=ddf_test$feve, pred=predict_test)
  metrics_test <- metrics(data_metrics_test, truth, pred)
  
  # plotting prediction
  library(patchwork)
  
  gg_test <- ggplot(data_metrics_test, aes(x=truth, y=pred))+
    geom_point()+
    labs(title = "Density of Data Point for Training Data",
         subtitle = paste("Random Forest: Rsq = ", round(metrics_test$.estimate[2], 3), sep = ""),
         x = paste("predicted values for ", i, sep = ""),
         y = paste("observed values for ", i, sep = "")) +
    theme(aspect.ratio=1)+
    geom_abline(intercept = 0, slope = 1, color = "#fc5e03")
  
  gg_train <- ggplot(data_metrics_train, aes(x=truth, y=pred))+
    geom_point()+ 
    labs(title = "Density of Data Point for Training Data",
         subtitle = paste("Random Forest: Rsq", round(metrics_train$.estimate[2], 3), sep = ""),
         x = paste("predicted values for ", i, sep = ""),
         y = paste("observed values for ", i, sep = "")) +
    theme(aspect.ratio=1)+
    geom_abline(intercept = 0, slope = 1, color = "#fc5e03")
  
  # arrange them into one file to export
  setwd(output)
  ggarrange(gg_test, gg_train, ncol = 2, nrow = 1,
            labels = c("A", "B"), common.legend = TRUE, legend = "right")
  ggsave(paste("./functional diversity/pasc_ID/training and predicting/random forest/", i, "_rf_prediction_pasc.png", sep = ""), width = 8, height = 5)
  setwd(input)
  
} # end i loop































# B.LAPIDARIUS ---- 
# only B.lapidarius 
BB22.bb.traits.sp <- BB22.bb.traits[BB22.bb.traits$bbspecies == "B.lapidarius",]

# import data with spatial information on sites
BB22.sites.meta <- read_csv("BB22_sites_2016.csv")

# import data with functional diversity of plants per site
BB22.fun.ID <- read_csv("./FD/FD_package_B.lapidarius_ID.short.csv")%>% 
  rename_with(.cols = 1, ~"ID")%>%
  rename_with(.cols = 2, ~"sp_richn")%>%
  rename_with(.cols = 3, ~"fric")%>%
  rename_with(.cols = 5, ~"fdiv")%>%
  rename_with(.cols = 4, ~"feve")

# filter functional metrics of plant traits to compare with traits per individual
BB22.fun.ID <- subset(BB22.fun.ID, ID %in% BB22.bb.traits.sp$ID)

# add site coordinates to the trait data frame (in LV95)
BB22.ID <- merge(BB22.bb.traits.sp, BB22.fun.ID, by  = "ID", all.x=TRUE)
BB22.ID <- merge(BB22.ID, BB22.sites.meta[, c(1,2,3)], by  = "site", all.x=TRUE) 

## look at data ----------------------------------------------------------------------------------------
# FDis, FRic, FDiv, FEve, FSpe
library(Hmisc)
hist.data.frame(BB22.ID[, c(17:20)])

# Boxplots for all the variables we want to look at with Wilcoxon test
resp <- c("sp_richn", "fric", "fdiv", "feve")
library(psych)
library(rstatix)
plot_list  <- list()
palette.landscape <- c("#E69F00", "#56B4E9") #create color palette for landscape

for (j in resp) {
  # j <- "sp_richn"
  gg.data <- data.frame(landscape=BB22.ID$landscape,value=BB22.ID[[j]])
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
  if (j == "sp_richn") {
    p <- p + ylim(0, 18) + ylab("species richness")
  } else if (j == "fric") {
    p <- p + ylim(0, 1) + ylab("funtional richness")
  } else if (j == "fdiv") {
    p <- p + ylim(0.6, 0.85) + ylab("funtional divergence")
  } else {
    p <- p + ylim(0.5, 0.75) + ylab("funtional evenness")
  }
  plot_list[[j]] <- p
  # describeBy(gg.data$value, gg.data$landscape)
}

# arrange them into one file to export
setwd(output)
plot <- ggarrange(plot_list[[1]],plot_list[[2]],
                  plot_list[[3]],plot_list[[4]],
                  ncol = 1, nrow = 6,
                  labels = c("A", "B", "C", "D"))
annotate_figure(plot, top = text_grob("B.lapidarius: species richness and funtional diversity across landscapes", 
                                      face = "bold", size = 22))
# ggsave("./functional diversity/lapi_ID/FD_B.lapidarius_ID.png", width = 6, height = 24)
setwd(input)



# plot the relationship of plants traits of one site and bumblebee traits of one site
traits <- colnames(BB22.ID[, 7:15]) # bumblebee traits to look at; prepare for loop
metrics <- colnames(BB22.ID[, c(17:20)]) # plant FD to look at (see file BB22_compute_FD); prepare for loop

library(nlme)

# perform loop to output plots per relationship summarized per FD
for (i in metrics) {
  x <- 1 # for naming the plots
  for (j in traits) {
    
    i <- "fric"
    j <- "corbicula_ratio"
    f <- formula(paste(i,"~", j))
    fit <- lme(f, random=~1|landscape, data = BB22.ID, na.action=na.omit)
    assign(paste("a", x, sep=""), # assign the ggplot to plot name
           # define the ggplot
           ggplot(BB22.ID, aes_string(j, i, colour = "landscape")) + 
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
  annotate_figure(plot4, top = text_grob(paste("B.lapidarius: comparison of ", i, " and traits across landscapes", sep = ""),
                                         face = "bold", size = 22))
  # ggsave(paste("./functional diversity/lapi_ID/FD_pasc_corr_", i, "_BBtraits_landscapes.png", sep = ""), width = 16, height = 8)
  setwd(input)
} # end loop i



## modelling ----------------------------------------------------------------------------------------
# center and scale all the variables
BB22.ID[, c(7:15, 17:20)] <- scale(BB22.ID[, c(7:15, 17:20)],center=TRUE,scale=TRUE)

# correlation analysis
# look at the correlation between the explanatory variables
library(corrplot)
M<-cor(BB22.ID[, 7:15], use = "complete.obs") # subset data; only explanatory variables

# matrix of the p-value of the correlation
p.mat <- cor.mtest(BB22.ID[, 7:15])

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


### Species Richness ----------------------------------------------------------------------------------------

# built an initial full model based on collinearity 
M1.full <- lmer(sp_richn ~ proboscis_ratio + fore_wing_ratio + corbicula_ratio + 
                  (1|landscape),
                data=BB22.ID)
fix.check(M1.full) # looks ok
vif(M1.full) # looks good


# formal test for spatial correlation
sims <- simulateResiduals(M1.full)
BB22.ID$site <- as.factor(BB22.ID$site)
simulationOutput <- recalculateResiduals(sims, group = BB22.ID$site)
testSpatialAutocorrelation(simulationOutput, x = BB22.sites.meta$LV95_x, y = BB22.sites.meta$LV95_y, plot = FALSE)
# there is no spatial autocorrelation

# update model: add intertegular distance to the model
M1.full.1 <- lmer(sp_richn ~ intertegular_distance + proboscis_ratio + fore_wing_ratio + 
                    corbicula_ratio + (1|landscape),
                  data=BB22.ID)
fix.check(M1.full.1)
summary(M1.full.1)
vif(M1.full.1) # looks good
anova(M1.full.1, M1.full) # with intertegular_distance better fit than M1.full

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


### Functional Richness ----------------------------------------------------------------------------------------
# built an initial full model based on colinearity 
M2.full <- lmer(fric ~ proboscis_ratio + fore_wing_ratio + corbicula_ratio + 
                  (1|landscape),
                data=BB22.ID)
fix.check(M2.full) # looks ok
vif(M2.full) # looks good
summary(M2.full)

# formal test for spatial correlation
sims <- simulateResiduals(M2.full)
BB22.ID$site <- as.factor(BB22.ID$site)
simulationOutput <- recalculateResiduals(sims, group = BB22.ID$site)
testSpatialAutocorrelation(simulationOutput, x = BB22.sites.meta$LV95_x, y = BB22.sites.meta$LV95_y, plot = FALSE)
# there is no spatial autocorrelation

# add intertegular distance to the model
M2.full.1 <- lmer(fric ~ intertegular_distance + proboscis_ratio + fore_wing_ratio + 
                    corbicula_ratio + (1|landscape),
                  data=BB22.ID) # boundary (singular) fit: see help('isSingular')
fix.check(M2.full.1)
summary(M2.full.1)
vif(M2.full.1) # looks good
anova(M2.full.1, M2.full) # with intertegular_distance better fit than M1.full

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

### Functional Divergence ----------------------------------------------------------------------------------------
# built an initial full model based on collinearity 
M3.full <- lmer(fdiv ~ proboscis_ratio + fore_wing_ratio + corbicula_ratio + (1|landscape),
                data=BB22.ID)
fix.check(M3.full) # looks ok
vif(M3.full) # looks good
summary(M3.full)

# formal test for spatial correlation
sims <- simulateResiduals(M3.full)
BB22.ID$site <- as.factor(BB22.ID$site)
simulationOutput <- recalculateResiduals(sims, group = BB22.ID$site)
testSpatialAutocorrelation(simulationOutput, x = BB22.sites.meta$LV95_x, y = BB22.sites.meta$LV95_y, plot = FALSE)
# there is no spatial autocorrelation

# add intertegular distance to the model
M3.full.1 <- lmer(fdiv ~ intertegular_distance + proboscis_ratio + fore_wing_ratio + 
                    corbicula_ratio + (1|landscape),
                  data=BB22.ID)
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

### Functional Evenness ----------------------------------------------------------------------------------------
# built an initial full model based on collinearity 
M4.full <- lmer(feve ~ proboscis_ratio + fore_wing_ratio + corbicula_ratio + (1|landscape),
                data=BB22.ID) #boundary (singular) fit: see help('isSingular')
fix.check(M4.full) # looks ok
vif(M4.full) # looks good
summary(M4.full)

# formal test for spatial correlation
sims <- simulateResiduals(M4.full)
BB22.ID$site <- as.factor(BB22.ID$site)
simulationOutput = recalculateResiduals(sims, group = BB22.ID$site)
testSpatialAutocorrelation(sims, x = BB22.sites.meta$LV95_x, y = BB22.sites.meta$LV95_y, plot = FALSE)
# there is no spatial autocorrelation

# add intertegular distance to the model
M4.full.1 <- lmer(feve ~ intertegular_distance + proboscis_ratio + fore_wing_ratio + 
                    corbicula_ratio + (1|landscape),
                  data=BB22.ID) # boundary (singular) fit: see help('isSingular')
fix.check(M4.full.1)
summary(M4.full.1)
vif(M4.full.1) # looks good
anova(M4.full.1, M4.full) # with intertegular_distance not a better fit than M1.full

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




