#####################################################################################
### Hahs and Fournier et al. 2022 - Compute functional diversity indices ------------
### Script to calculate various functional diversity metrics
### Contact: Amy Hahs <amy.hahs@unimelb.edu.au>
#####################################################################################


rm(list = ls())


### load useful function ("Toolkit") -------------------------------------------
source("ADAPT_PATH/Function_ToolKit.R")


### Load data related to a taxa ------------------------------------------------
setwd("ADAPT_PATH/Amphibian")
sp <- read.delim(file="sp.txt", row.names = 1, header=T)
site <- read.delim(file="site.txt")
trt <- read.delim(file="trt.txt", row.names = 1)

# remove empty sites
site <- site[rowSums(sp)>0,]
sp <- sp[rowSums(sp)>0,]

# remove species with no occurrence
trt = trt[colSums(sp)>0, ]
sp = sp[ ,colSums(sp)>0]

dim(sp); dim(site); dim(trt)

# compute presence/absence
require(vegan)
sp.pa <- decostand(sp, "pa")



### Calculate species richness and Shannon diversity ------------------------------------------------

# Species richness
rich = specnumber(sp.pa)

# Shannon
sha <- diversity(sp.pa)



### calculate CWMs ------------------------

# CWM
CWM = CWM_calc(trt, site, sp.pa)

for(i in 1:ncol(CWM)){
  print(i)
  CWM[,i][is.na(CWM[,i])]=mean(na.omit(CWM[,i]))
}


### Imput trait mnissing values and compute PCA ----------------------------------

# Imput missing values and reduce data dimensionality using PCA 
# -> enable to calculate FD metrics for sites with low diversity
require(caret)
require(vegan)
trt.mis.model = preProcess(trt, "knnImpute")
trt.mis.pred = predict(trt.mis.model, trt); head(trt.mis.model)

# PCA -> reduce dimensionality
trt.pca <- prcomp(trt.mis.pred, scale. = T, center = T)
cumsum(trt.pca$sdev/sum(trt.pca$sdev))
trt.scaled <- scores(trt.pca)[,1:2] # adjust number of axes for each group


### calculate TOP and TED from Fontana et al. 2016 ------------------------

### TOP
TOP.mat <- TOP_over(trt = trt.scaled, sp=sp.pa)
TOP = TOP.mat[,2]
TOP[is.na(TOP)] = 0



### TED

# List of possible REFERENCES
# Define maximum number of points (max1) and number of trts under consideration (dim1)
# Alternatively, it is possible to manually define max1 and dim1!!
max1 <- max(specnumber(sp))
dim1 <- ncol(trt)

# create reference matrix
ref.matrix<-matrix(ncol=2,nrow=max1)
if (dim1 == 1) {
  i=0.9 } else { i=1.9 }
n <- 0
rows1<-0

while(rows1<max1){
  i=i+0.1
  n=n+1
  trts.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(trts.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
}

k <- i+1
while(i<k){
  i=i+0.1
  n=n+1
  trts.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(trts.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
}

ref.matrix<-na.omit(ref.matrix)

# compute TED
TED.mat <- TED_over(trt = trt.scaled, sp=sp)
TED = TED.mat[,1]
TED[is.na(TED)] = 1


### FD incices from mFD package of Magneville et al. 2022 ------------------------

### mFD
FD_mFD <- calc_mFD(trt=trt.mis.pred, sp=sp.pa)




### FDis from package FD and Rao functional dispersion ------------------------

# FD: Rao and FDis
require(cluster)
require(mgcv)
require(FD)

trt.dist <- fundist_spe(trt=trt.mis.pred, sp=sp.pa)
dis <- as.matrix(trt.dist)

# Rao functional dispersion - de Bello et al. 2010
dat.rao <- Rao(sample=t(sp.pa), dfunc=dis, dphyl=NULL, weight=FALSE, Jost=TRUE, structure=NULL)
FDrao <- dat.rao$FD$Alpha

# FDis - Laliberté and Legendre (2010)
FDis <- fdisp(trt.dist, as.matrix(sp.pa))
FDis <- FDis$FDis


### FEver from Ricotta et al. 2014 ------------------------
FEve.Ricotta <- FeveR(sp.pa, trt.dist)



### FD indices from Villegier et al. 2008 --------------------------------------
require(FD)
ex3 <- dbFD(trt.scaled, sp.pa, calc.FRic=F)
FEve.Villegier = ex3$FEve
FDiv.Villegier = ex3$FDiv # not used in manuscript
FDis.Villegier = ex3$FDis # not used in manuscript
RaoQ.Villegier = ex3$RaoQ # not used in manuscript



