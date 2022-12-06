#####################################################################################
### Hahs and Fournier et al. 2022 - Urban functional diversity ----------------------
### Collection of utility functions used in FD metric computation 
### and random forests analyses
### Contact: Amy Hahs <amy.hahs@unimelb.edu.au>
#####################################################################################


# Utility functions based on the package "mFD" of Magneville et al. 2022 ----------

# functional distance among species
fundist_spe <- function(trt, sp){
  require("mFD")
  
  trt.cat = data.frame(trait_name=colnames(trt), 
                       trait_type = rep("Q", ncol(trt)),                       
                       trait_type = rep(1, ncol(trt)),
                       fuzzy_name =NA
  )
  
  
  sp_dist_trt <- mFD::funct.dist(
    sp_tr         = trt,
    tr_cat        = trt.cat,
    metric        = "gower",
    scale_euclid  = "scale_center",
    ordinal_var   = "classic",
    weight_type   = "equal",
    stop_if_NA    = TRUE)
  
  return(sp_dist_trt)  
}


# FD indice from mFD
calc_mFD <- function(trt, sp){
  
  require("mFD")
  
  trt.cat = data.frame(trait_name=colnames(trt), 
                       trait_type = rep("Q", ncol(trt)),                       
                       trait_type = rep(1, ncol(trt)),
                       fuzzy_name =NA
  )
  
  
  sp_dist_trt <- mFD::funct.dist(
    sp_tr         = trt,
    tr_cat        = trt.cat,
    metric        = "gower",
    scale_euclid  = "scale_center",
    ordinal_var   = "classic",
    weight_type   = "equal",
    stop_if_NA    = TRUE)
  
  
  FD2max <- mFD::alpha.fd.hill(
    asb_sp_w = as.matrix(sp), 
    sp_dist  = sp_dist_trt, 
    tau      = "max", 
    q        = 2)
  
  ### Compute multimensional functional spaces and assess their quality
  fspaces <- mFD::quality.fspaces(
    sp_dist             = sp_dist_trt,
    maxdim_pcoa         = 10,
    deviation_weighting = "absolute",
    fdist_scaling       = FALSE,
    fdendro             = "average")
  
  sel.axis = which(fspaces$"quality_fspaces" == min(fspaces$"quality_fspaces"))
  sel.axis=2
  fspace.sel <- fspaces$details_fspaces$sp_pc_coord[,1:sel.axis]
  
  
  alpha_fd <- mFD::alpha.fd.multidim(
    sp_faxes_coord   = fspace.sel,
    asb_sp_w         = as.matrix(sp[rich>sel.axis,]),
    ind_vect         = c("fdis", "feve", "fric"),
    scaling          = TRUE,
    check_input      = TRUE,
    details_returned = TRUE)
  
  dim(alpha_fd$functional_diversity_indices)
  dim(sp)
  
  FD_mat_Magneville <- as.data.frame(matrix(NA, 
                                            ncol=ncol(alpha_fd$functional_diversity_indices), 
                                            nrow=nrow(sp)
  ))
  colnames(FD_mat_Magneville) = colnames(alpha_fd$functional_diversity_indices)
  rownames(FD_mat_Magneville) = rownames(sp)
  
  for(i in 1:nrow(sp)){
    if(rownames(FD_mat_Magneville)[i] %in% rownames(alpha_fd$functional_diversity_indices)){
      row.sel = which(rownames(alpha_fd$functional_diversity_indices) == rownames(FD_mat_Magneville)[i])
      FD_mat_Magneville[i,] = alpha_fd$functional_diversity_indices[row.sel,]  
    }else FD_mat_Magneville[i,] = rep(NA, ncol(FD_mat_Magneville))
  }
  
  return(data.frame(FD_mat_Magneville, Hill_2_max = FD2max$asb_FD_Hill))
}





# Utility function to calculate the functional evenness "FEver" based on --------

# R function 'FeveR' for calculating the functional evenness of a species' assemblage. 
FeveR<-function(abundances, distances)
{
  rel_abundance<-sweep(abundances, 1, rowSums(abundances), "/")
  rel_abu_matrix<-as.matrix(rel_abundance)
  n_plot<-nrow(rel_abu_matrix)
  n_species<-ncol(rel_abu_matrix)
  dista_matrix<-as.matrix(distances)
  ###BULLA1
  index_array<-rep(NA, n_plot)
  for(i in 1:n_plot)
  {mat_develop<-rep(rel_abu_matrix[i,],n_species)
  matrix_two<-(t(matrix(mat_develop, nrow=n_species,ncol=n_species)))
  matrix_non<-matrix(mat_develop, nrow=n_species,ncol=n_species)
  uni<-matrix(1,nrow=n_species,ncol=n_species)
  subtraction<- uni-matrix_non
  division<-matrix_two/subtraction
  moltiplication<-division*dista_matrix
  row_sums<-rowSums(moltiplication)
  per_abundance<-row_sums*rel_abu_matrix[i,]
  great_sum<-sum(per_abundance)
  divised<-per_abundance/great_sum
  espress<-which(per_abundance>0)
  S<-length(espress)
  bulla<-rep(NA, length(divised))
  for(l in 1:length(divised))
  {
    bulla[l]<-min(divised[l], 1/S)
  }
  indice<-sum(bulla)
  norm_index<-(indice-(1/S))/(1-(1/S))
  index_array[i]<-norm_index
  }
  return(index_array)
}





# Utility functions to calculate the TOP, TED, and FDis indices based on Fontana et al. 2016 --------

# TOP index 
library(geometry)


TOP.index <- function(traitdat){
  
  # TOP
  
  dim1 <- ncol(traitdat)
  
  #definitions: index i, area as empty vector
  
  i=0
  area<-matrix(ncol=2,nrow=nrow(traitdat))
  
  while(nrow(traitdat)>dim1){
    i=i+1
    
    # use of convhulln function
    
    # area
    area[i,2] <- convhulln(traitdat, c("FA","QJ"))$area
    
    # identity of vertices
    vert0<-convhulln(traitdat, c("Fx TO 'vert.txt'","QJ"))
    vert1<-scan("vert.txt",quiet=T)
    vert2<-vert1+1
    
    vertices <- vert2[-1]
    
    traitdat <- traitdat[-vertices,]
    
    area[i,1] <- length(vertices)
    
  }
  
  
  area<-na.omit(area)
  
  # Output (2 numbers): Number of points touched by areas; Sum of the areas (TOP index)
  colSums(area)
  
}

# function to compute TOP for all communities
TOP_over <- function(trt, sp){
  TOP.mat <- matrix(NA, nrow(sp), 2)
  colnames(TOP.mat) <- c("Nb_points", "TOP")
  rownames(TOP.mat) <- rownames(sp)
  for (i in 1:nrow(sp)) {
    if(length(colnames(sp)[sp[i,]>0]) > ncol(trt)){
      trtdat <- as.data.frame(trt.scaled[rownames(trt.scaled) %in% colnames(sp)[sp[i,]>0],])
      TOP.mat[i,] <- TOP.index(trtdat)
    }else{
      TOP.mat[i,] <- c(NA,NA)
    }
  }
  return(TOP.mat)
}


# TED index
# comparison between sample and an assemblage of equidistant points (with the same number of individuals)

library(geozoo)
library(flexmix)
library(geometry)

#  TED index calculation 

TED.index <- function(traitdat){
  
  ##########################################################################################
  # Find the best REFERENCE (minimum number of individuals >= individuals in the sample)
  ##########################################################################################
  
  n.sample<-nrow(traitdat)
  
  diff1<-matrix(ncol=2,nrow=length(ref.matrix)/2)
  diff1[,1] <- ref.matrix[,1]
  diff1[,2] <- ref.matrix[,2]-n.sample
  min.diff1<-min(diff1[,2][diff1[,2]>=0])
  select.i<-diff1[diff1[,2]==min.diff1,][1]
  traits.ref <- sphere.solid.grid(p=dim1, n=select.i)
  
  ###################################
  # Transform REFERENCE in data frame
  ###################################
  
  traits.ref <- as.vector(traits.ref$points)
  ind<-length(traits.ref)/dim1
  reference<-matrix(ncol=dim1,nrow=ind)
  for (j in 1:dim1){
    reference[,j] <- traits.ref[((j-1)*ind+1):(j*ind)]
  }
  traits.ref <- as.data.frame(reference)
  
  
  ############################################################################
  # Ev. delete individuals in order to have the same number as in the sample
  ############################################################################
  
  x <- nrow(traits.ref)-nrow(traitdat)
  
  if (x!=0){
    
    # coordinates of the center of gravity of the vertices (Gv)
    baryv<-apply(traits.ref,2,mean)
    
    # euclidian dstances to Gv (dB)
    distbaryv<-rep(0,nrow(traits.ref))
    for (j in 1:nrow(traits.ref))
      distbaryv[j]<-( sum((traits.ref[j,]-baryv)^2) )^0.5
    
    merge1<-data.frame(traits.ref,distbaryv)
    
    #sort by distbaryv (descending)
    sort1 <- merge1[order(-distbaryv),]
    traits.ref<-sort1[-1:-x,-(ncol(sort1))]
    
  }
  
  
  #######################
  # Compare with sample
  #######################
  
  Distance.method <- "euclidean"
  D1 <- dist(traitdat, method=Distance.method)
  density.D <- density(D1)$y
  rm(D1)
  D.ref <- dist(traits.ref, method=Distance.method)
  density.D.ref <- density(D.ref)$y
  rm(D.ref)
  
  results <- KLdiv(cbind(density.D, density.D.ref))
  
  value <- results[1,2]
  
  TED <- 1-log10(value+1)
  TED
  
}


TED_over <- function(trt, sp){
  TED.mat <- matrix(NA, nrow(sp), 1)
  colnames(TED.mat) <- c("TED")
  rownames(TED.mat) <- rownames(sp)
  for (i in 1:nrow(sp)) {
    if(length(colnames(sp)[sp[i,]>0]) > 2){
      trtdat <- as.data.frame(trt.scaled[rownames(trt.scaled) %in% colnames(sp)[sp[i,]>0],])
      TED.mat[i,] <- TED.index(trtdat)
    }else{
      TED.mat[i,] <- c(NA)
    }
  }
  return(TED.mat)
}





# Utility function to calculated CWM based on Peres-Neto et al 2016 ---------------

# Appendix 3: R code for simulations, calculation of correlations and permutation tests.
generate_community <- function(tolerance,E,T,preset_nspecies,preset_ncommunities){
  repeat {
    # one trait, one environmental variable
    h <- runif(preset_nspecies,min=0.3,max=1)
    sigma <- runif(preset_nspecies)*tolerance
    L <- matrix(data=0,nrow=preset_ncommunities,ncol=preset_nspecies)
    for(j in 1:preset_nspecies){
      L[,j] <- rpois(preset_ncommunities,30*h[j]*exp(-(E-T[j])^2/(2*sigma[j]^2)))
    }
    n_species_c <- sum(colSums(L)!=0) # _c for check
    n_communities_c <- sum(rowSums(L)!=0)
    if ((n_species_c == preset_nspecies) & (n_communities_c==preset_ncommunities)){break}
  }
  return(L)
}


library(ade4)

CWM_gen <- function(L,T, Chessel = TRUE){
  
  T<-as.matrix(T)
  L<-as.matrix(L)
  #  centering_mat <- function(X,w){ X - rep(1,length(w))%*%t(w)%*% X }
  standardize_w <- function(X,w){
    ones <- rep(1,length(w))
    Xc <- X - ones %*% t(w)%*% X
    Xc / ones%*%sqrt(t(ones)%*%(Xc*Xc*w)) 
  } 
  
  # check_L()
  rows<-seq_len(nrow(L))
  cols<-seq_len(ncol(L))
  rni <-which(rowSums(L)==0)
  repeat {
    if (length(rni)) {L <- L[-rni,,drop = FALSE]; rows <-rows[-rni]}
    ksi <- which(colSums(L)==0)
    if (length(ksi)) {L <- L[,-ksi, drop = FALSE]; cols <- cols[-ksi]}
    rni <-which(rowSums(L)==0)
    if ( length(rni)==0 & length(ksi)==0){break}
  }
  T <-T[cols,,drop = FALSE]
  # end check_L()
  
  L<-L/sum(L)
  # dimensions
  #S <- ncol(L) # number of species
  #n <- nrow(L) # number of communities
  q <- ncol(T) # number of traits
  
  # setting up matrices
  Wn <- rowSums(L)
  Ws <- colSums(L)
  
  # cor matrices are trait by environment
  CWM <- L%*%T/Wn  # weighted means wrt to T
  return(data.frame(CWM=CWM))
}


CWM_calc <- function(trt, site, sp){
  CWM = matrix(NA, nrow(site), ncol(trt))
  colnames(CWM) = colnames(trt)
  rownames(CWM) = rownames(site)
  for(i in 1:ncol(trt)) CWM[rowSums(sp[,!(is.na(trt[,i]))])>0,][,i] <- CWM_gen(L=sp[,!(is.na(trt[,i]))],T=na.omit(trt[,i]))[,1]
  return(CWM)
}


TraitEnvCor <- function(L,E,T, Chessel = TRUE){
  
  E<-as.matrix(E)
  T<-as.matrix(T)
  L<-as.matrix(L)
  #  centering_mat <- function(X,w){ X - rep(1,length(w))%*%t(w)%*% X }
  standardize_w <- function(X,w){
    ones <- rep(1,length(w))
    Xc <- X - ones %*% t(w)%*% X
    Xc / ones%*%sqrt(t(ones)%*%(Xc*Xc*w)) 
  } 
  
  # check_L()
  rows<-seq_len(nrow(L))
  cols<-seq_len(ncol(L))
  rni <-which(rowSums(L)==0)
  repeat {
    if (length(rni)) {L <- L[-rni,,drop = FALSE]; rows <-rows[-rni]}
    ksi <- which(colSums(L)==0)
    if (length(ksi)) {L <- L[,-ksi, drop = FALSE]; cols <- cols[-ksi]}
    rni <-which(rowSums(L)==0)
    if ( length(rni)==0 & length(ksi)==0){break}
  }
  E <-E[rows,,drop = FALSE]
  T <-T[cols,,drop = FALSE]
  # end check_L()
  
  L<-L/sum(L)
  # dimensions
  #S <- ncol(L) # number of species
  #n <- nrow(L) # number of communities
  p <- ncol(E) # number of environmental predictors
  q <- ncol(T) # number of traits
  
  # setting up matrices
  Wn <- rowSums(L)
  Ws <- colSums(L)
  # cor matrices are trait by environment
  CWM <- L%*%T/Wn  # weighted means wrt to T
  CWM.cor <- cor(CWM,E)
  
  SNC <- t(L)%*%E/Ws  # weighted means wrt to E
  SNC.cor <- cor(T,SNC)
  
  CWMstd_w  <- standardize_w(CWM,Wn)
  Estd_w <- standardize_w(E,Wn)
  wCWM.cor <- t(t(Estd_w)%*%(CWMstd_w*Wn))
  
  SNCstd_w <- standardize_w(SNC,Ws)
  Tstd_w <-  standardize_w(T,Ws)
  wSNC.cor <- t(Tstd_w)%*%(SNCstd_w*Ws)
  
  # Fourth corner calculated as W_n weighted covariance between 
  # CWM and standardized T (trait)
  
  CWM_std_tw <- L%*%Tstd_w/Wn #CWM wrt to standardized T (trait)
  Fourthcorner <- t(CWM_std_tw)%*%(Estd_w*Wn)
  if (Chessel){
    singular_val1 <- sqrt(ade4::dudi.coa(L, scannf = FALSE)$eig[1])
    Chessel.4thcor<-Fourthcorner/ singular_val1
  }else { Chessel.4thcor<-NA;singular_val1 <-1}
  
  
  # variation components
  # Among communities 
  Among.Variation <- sum(diag(t(CWM_std_tw)%*%(CWM_std_tw* Wn)))
  # Within communities 
  Within.Variation <- 1 - Among.Variation
  
  # result specialized to one trait and one environment variables; use array(0, dim(6,k,p)) in the general case
  # array.result<-matrix(c(CWM.cor,wCWM.cor,SNC.cor,wSNC.cor,Fourthcorner,Chessel.4thcor,Mean.Variation),ncol=1)
  array.result<-array(0, dim=c(8,q,p))
  rownames(array.result)<- c("CWM.cor","wCWM.cor","SNC.cor","wSNC.cor","Fourthcorner","Chessel.4thcor","Among Wn-variance (%)", "Within Wn-variance (%)")
  array.result[1,,]<-CWM.cor
  array.result[2,,]<-wCWM.cor
  array.result[3,,]<-SNC.cor
  array.result[4,,]<-wSNC.cor
  array.result[5,,]<-Fourthcorner
  array.result[6,,]<-Chessel.4thcor
  array.result[7,,]<-Among.Variation * 100
  array.result[8,,]<-Within.Variation * 100
  return(array.result[,,])
}


CorPermutationTest <- function(L, E, T, nrepet = 999){
  E<-as.matrix(E)
  T<-as.matrix(T)
  L<-as.matrix(L)
  obs <- TraitEnvCor(L,E,T)[1:5]
  sim.row <- matrix(0, nrow = nrepet, ncol = ncol(E) * 5)
  sim.col <- matrix(0, nrow = nrepet, ncol = ncol(E) * 5)
  for(i in 1:nrepet){
    per.row <- sample(nrow(L))
    per.col <- sample(ncol(L))
    sim.row[i, ] <- c(as.matrix(data.frame(TraitEnvCor(L,E[per.row,,drop= FALSE],T))))[1:5]
    sim.col[i, ] <- c(as.matrix(data.frame(TraitEnvCor(L,E,T[per.col,,drop= FALSE]))))[1:5]
  }
  pval.row <- (rowSums(apply(sim.row^2, 1, function(i) i >= obs^2)) + 1)  / (nrepet + 1)
  pval.col <- (rowSums(apply(sim.col^2, 1, function(i) i >= obs^2)) + 1)  / (nrepet + 1)
  
  result <- cbind(cor = obs, prow = pval.row, pcol = pval.col, pmax = apply(cbind(pval.row, pval.col), 1, max))
  return(result)
}

CorPermutationTest_fast <- function(L, E, T, nrepet = 999){
  E<-as.matrix(E)
  T<-as.matrix(T)
  L<-as.matrix(L)
  obs <- TraitEnvCor(L,E,T)[1:5]
  sim.row <- matrix(0, nrow = nrepet, ncol = ncol(E) * 5)
  for(i in 1:nrepet){
    per.row <- sample(nrow(L))
    sim.row[i, ] <- c(as.matrix(data.frame(TraitEnvCor(L,E[per.row,,drop= FALSE],T))))[1:5]
  }
  pval.row <- (rowSums(apply(sim.row^2, 1, function(i) i >= obs^2)) + 1)  / (nrepet + 1)
  
  result <- cbind(cor = obs, prow = pval.row)
  return(result)
}



# Utility function to calculate Rao functional dispersion based on de Bell --------
# 	The Rao function computes alpha, gamma and beta-components for taxonomic, functional and phylogenetic diversity with the Rao index          
# 	The script integrates two functions: "Qdecomp", by Villeger & Mouillot (J Ecol, 2008) modify by Wilfried Thuiller, and "disc", by S. Pavoine, in the package ade4.
# 	For a regional assemblage of C local communities gamma = mean(alpha) + beta, where:
#  	gamma is the diversity of the regional pool
#  	alpha are the diversities of the local communities
#  	beta is the turn over between local communities
#  	diversity is estimated with the Rao quadratic entropy index (Rao 1982)
#                                      
# INPUTS:                                                                                 
#	- "abundances": matrix of abundances (c x s) of the s species for the c local communities (or samples)           
#	- "dfunct": matrix (s x s) or dist object with pairwise functional trait distances between the s species
#	- "dphyl": as dfunct but for phylogenetic distances
#	- "weight": defining if the correction by Villeger & Mouillot (J Ecol, 2008) is applied or not
#	- "Jost": defining if the correction Jost correction is applied (this paper and Jost 2007)
#	- "structure": a data frame containing the name of the group to which samples belong see                                
#      NA are not allowed in 'functdist'
#      NA are automatically replaced by 0 in 'abundances'
#                                                                                         
# OUTPUTS:                                                                           
#	- The results are organized for Taxonomic diversity ($TD), Functional diversity ($FD) and phylogenetical diversity ($PD). Beta and gamma diversities are calculated for the whole data set and for each pair of samples ("Pairwise_samples") 
#	- "$Richness_per_plot"(number of species per sample)
#	- "$Relative_abundance" (species relative abundances per plot)
#	- "$Pi" (species regional relative abundance)
#	- "$Wc" (weigthing factor),                               
#	- "$Mean_Alpha" (mean aplpha diversity; for taxonomic diversity the Simpson index is calculated)                                   
#	- "$Alpha" (alpha diversity for each sample; for taxonomic diversity the Simpson index is calculated)                                       
#	- "$Gamma" (gamma diversity; for taxonomic diversity the Simpson index is calculated)                                      
#	- "$Beta_add" (Gamma-Mean_Alpha)                                          
#	- "$Beta_prop" (Beta_add*100/Gamma)                                                         
#	- "$Pairwise_samples$Alpha" (mean alpha for each pair of samples)                        
#	- "$Pairwise_samples$Gamma" (gamma for each pair of samples)
#	- "$Pairwise_samples$Beta_add" (beta for each pair of samples as Gamma-Mean_Alpha)  
#	- "$Pairwise_samples$Beta_prop" (beta for each pair of samples as Beta_add*100/Gamma)  

Rao<-function(sample, dfunc, dphyl, weight=F, Jost=F, structure=NULL)   {
  library(ade4)
  
  ####function Qdecomp by by VillÃger & Mouillot (J Ecol, 2008) modify by Wilfried Thuiller #####
  
  Qdecomp = function(functdist,abundances, w=TRUE) {
    
    # number and names of local communities
    c<-dim(abundances)[1] ; namescomm<-row.names(abundances)
    abundances<-as.matrix(abundances)
    
    # if necessary, transformation of functdist into matrix object
    if (is.matrix(functdist)==F) functdist<-as.matrix(functdist)
    
    # checking 'abundances' and 'functdist' dimensions
    if (dim(functdist)[1]!=dim(functdist)[2])  stop("error : 'functdist' has different number of rows and columns")
    if (dim(abundances)[2]!=dim(functdist)[1]) stop("error : different number of species in 'functdist' and 'abundances' ")
    
    # checking NA absence in 'functdist'
    if (length(which(is.na(functdist)==T))!=0)  stop("error : NA in 'functdist'")
    
    # replacement of NA by 0 in abundances
    if (is.na(sum(abundances))==T)  {
      for (i in 1:dim(abundances)[1])
        for (j in 1:dim(abundances)[2] )
        { if(is.na(abundances[i,j])==T) abundances[i,j]<- 0 } # end of i j
    } # end of if
    
    #  species richness and total abundances in local communities
    abloc<-apply(abundances,1,sum)
    nbsploc<-apply(abundances,1,function(x) {length(which(x>0))} )
    
    # relative abundances inside each local community
    locabrel<-abundances/abloc
    
    # alpha diversity
    Qalpha=apply(locabrel, 1, function(x) t(x) %*%  functdist %*% x)
    
    #Wc
    Wc = abloc/sum(abloc)
    
    # abundance-weighted mean alpha
    mQalpha<-as.numeric(Qalpha%*%abloc/sum(abloc) )
    
    #Villeger's correction
    if(w==T) {
      # abundance-weighted mean alpha
      mQalpha<-as.numeric(Qalpha%*%abloc/sum(abloc) )
      totabrel<-apply(abundances,2,sum)/sum(abundances) 
      Qalpha = Qalpha*Wc
    }	
    
    # Rao's original definition: mean of Pi
    else {
      mQalpha<-mean(Qalpha)
      totabrel<-apply(locabrel,2,mean)  
    }
    
    
    # gamma diversity
    Qgamma<-( totabrel %*% functdist %*% totabrel ) [1]
    
    # beta diversity
    Qbeta<-as.numeric( Qgamma-mQalpha )
    
    # standardized beta diversity
    Qbetastd<-as.numeric(Qbeta/Qgamma )
    
    # list of results
    resQ<-list(Richness_per_plot=nbsploc, Relative_abundance= locabrel, Pi=totabrel, Wc=Wc, Species_abundance_per_plot=abloc, Alpha=Qalpha, Mean_alpha=mQalpha, Gamma=Qgamma, Beta=Qbeta, Standardize_Beta =Qbetastd )
    
    return(resQ)
    
  } 
  
  
  ###########function disc originally from S. Pavoine####
  
  disc = function (samples, dis = NULL, structures = NULL, Jost = F)
  {
    if (!inherits(samples, "data.frame"))
      stop("Non convenient samples")
    if (any(samples < 0))
      stop("Negative value in samples")
    if (any(apply(samples, 2, sum) < 1e-16))
      stop("Empty samples")
    if (!is.null(dis)) {
      if (!inherits(dis, "dist"))
        stop("Object of class 'dist' expected for distance")
      # if (!is.euclid(dis))
      #stop("Euclidean property is expected for distance")
      dis <- as.matrix(dis)
      if (nrow(samples) != nrow(dis))
        stop("Non convenient samples")
    }
    if (!is.null(structures)) {
      if (!inherits(structures, "data.frame"))
        stop("Non convenient structures")
      m <- match(apply(structures, 2, function(x) length(x)),
                 ncol(samples), 0)
      if (length(m[m == 1]) != ncol(structures))
        stop("Non convenient structures")
      m <- match(tapply(1:ncol(structures), as.factor(1:ncol(structures)),
                        function(x) is.factor(structures[, x])), TRUE, 0)
      if (length(m[m == 1]) != ncol(structures))
        stop("Non convenient structures")
    }
    Structutil <- function(dp2, Np, unit, Jost) {
      if (!is.null(unit)) {
        modunit <- model.matrix(~-1 + unit)
        sumcol <- apply(Np, 2, sum)
        Ng <- modunit * sumcol
        lesnoms <- levels(unit)
      }
      else {
        Ng <- as.matrix(Np)
        lesnoms <- colnames(Np)
      }
      sumcol <- apply(Ng, 2, sum)
      Lg <- t(t(Ng)/sumcol)
      colnames(Lg) <- lesnoms
      Pg <- as.matrix(apply(Ng, 2, sum)/nbhaplotypes)
      rownames(Pg) <- lesnoms
      deltag <- as.matrix(apply(Lg, 2, function(x) t(x) %*%
                                  dp2 %*% x))
      ug <- matrix(1, ncol(Lg), 1)
      if(Jost) {
        #dp2 <- as.matrix(as.dist(dfunct01))
        deltag <- as.matrix(apply(Lg, 2, function(x) t(x) %*% dp2 %*% x))
        X=t(Lg) %*% dp2 %*% Lg
        alpha=1/2 * (deltag %*% t(ug) + ug %*% t(deltag))
        Gam = (X + alpha)/2
        alpha = 1/(1-alpha) #Jost correction
        Gam = 1/(1-Gam)  #Jost correction
        Beta_add = Gam - alpha
        Beta_mult = 100*(Gam - alpha)/Gam
      }
      else {
        deltag <- as.matrix(apply(Lg, 2, function(x) t(x) %*% dp2 %*% x))
        X=t(Lg) %*% dp2 %*% Lg
        alpha=1/2 * (deltag %*% t(ug) + ug %*% t(deltag))
        Gam = (X + alpha)/2
        Beta_add = Gam - alpha
        Beta_mult = 100*(Gam - alpha)/Gam
      }
      colnames(Beta_add) <- lesnoms
      rownames(Beta_add) <- lesnoms
      return(list(Beta_add = as.dist(Beta_add), Beta_mult = as.dist(Beta_mult),
                  Gamma=as.dist(Gam), Alpha=as.dist(alpha), Ng = Ng, Pg = Pg))
    }
    Diss <- function(dis, nbhaplotypes, samples, structures, Jost) {
      structutil <- list(0)
      structutil[[1]] <- Structutil(dp2 = dis, Np = samples, NULL, Jost)
      diss <- list(structutil[[1]]$Alpha, structutil[[1]]$Gamma, structutil[[1]]$Beta_add, structutil[[1]]$Beta_mult)
      if (!is.null(structures)) {
        for (i in 1:length(structures)) {
          structutil[[i + 1]] <- Structutil(as.matrix(structutil[[1]]$Beta_add), 
                                            structutil[[1]]$Ng, structures[, i], Jost)
        }
        diss <- c(diss, tapply(1:length(structures), factor(1:length(structures)), 
                               function(x) as.dist(structutil[[x + 1]]$Beta_add)))
      }    
      return(diss)
    }
    nbhaplotypes <- sum(samples)
    diss <- Diss(dis, nbhaplotypes, samples, structures, Jost)
    if (!is.null(structures)) {
      names(diss) <- c("Alpha", "Gamma", "Beta_add", "Beta_prop", "Beta_region")
      return(diss)
    }
    names(diss) <- c("Alpha", "Gamma", "Beta_add", "Beta_prop")
    return(diss)
  }
  
  
  
  
  
  
  TD<-FD<-PD<-NULL
  
  #Taxonomic diversity
  dS <- matrix(1, nrow(sample), nrow(sample)) - diag(rep(1, nrow(sample)))
  temp_qdec<- Qdecomp(dS,t(sample), w=weight)   #Call the Qdecomp function for alpha, gamma and beta estimations.
  TD$Richness_per_plot = temp_qdec$Richness_per_plot
  TD$Relative_abundance = temp_qdec$Relative_abundance
  TD$Pi = temp_qdec$Pi
  TD$Wc = temp_qdec$Wc
  if(Jost){
    TD$Mean_Alpha = 1/(1-temp_qdec$Mean_alpha)
    TD$Alpha = 1/(1-temp_qdec$Alpha)
    TD$Gamma = 1/(1-temp_qdec$Gamma)
    TD$Beta_add = (TD$Gamma -TD$Mean_Alpha )
    TD$Beta_prop = 100*TD$Beta_add/TD$Gamma
    #Call the disc function for alpha, gamma and beta estimations for each pair of samples
    TD$Pairwise_samples<- disc(as.data.frame(sample), as.dist(dS), structure=structure, Jost=Jost)
  }
  else {
    TD$Mean_Alpha = temp_qdec$Mean_alpha
    TD$Alpha = temp_qdec$Alpha
    TD$Gamma = temp_qdec$Gamma
    TD$Beta_add = (TD$Gamma -TD$Mean_Alpha )
    TD$Beta_prop = 100*TD$Beta_add/TD$Gamma
    #Call the disc function for alpha, gamma and beta estimations for each pair of samples
    # TD$Pairwise_samples <- disc(as.data.frame(sample), as.dist(dS), structure=structure, Jost=Jost)
  }
  
  #Functional diversity estimation
  if(!is.null(dfunc)){
    FD<-list()
    if(Jost){
      if(max(dfunc)>1) dfunc <- dfunc/max(dfunc)   #Make sure the distance are between 0 and 1 for the Jost correction
      temp_qdec<- Qdecomp(dfunc,t(sample), w=weight)   #Call the Qdecomp function for alpha, gamma and beta estimations.
      #  FD$Alpha = 1/(1-temp_qdec$Alpha)
      #  FD$Mean_Alpha = mean(FD$Alpha)
      FD$Mean_Alpha = 1/(1-temp_qdec$Mean_alpha)
      FD$Alpha = 1/(1-temp_qdec$Alpha)
      FD$Gamma = 1/(1-temp_qdec$Gamma)
      FD$Beta_add = (FD$Gamma -FD$Mean_Alpha )
      FD$Beta_prop = 100*FD$Beta_add/FD$Gamma
      #Call the disc function for alpha, gamma and beta estimations for each pair of samples
      FD$Pairwise_samples<- disc(as.data.frame(sample), as.dist(dfunc), structure=structure, Jost=Jost)
    }
    else {
      temp_qdec<- Qdecomp(dfunc,t(sample), w=weight) #Call the Qdecomp function for alpha, gamma and beta estimations.
      FD$Mean_Alpha = temp_qdec$Mean_alpha
      FD$Alpha = temp_qdec$Alpha
      FD$Gamma = temp_qdec$Gamma
      FD$Beta_add = (FD$Gamma -FD$Mean_Alpha )
      FD$Beta_prop = 100*FD$Beta_add/FD$Gamma
      #FD$Beta =  temp_qdec$Beta#
      #Call the disc function for alpha, gamma and beta estimations for each pair of samples
      # FD$Pairwise_samples <- disc(as.data.frame(sample), as.dist(dfunc), structure=structure, Jost=Jost)
    }
  }
  #Phylogenetic diversity estimation
  if(!is.null(dphyl)){
    PD<-list()
    if(Jost){
      if(max(dphyl)>1) dphyl <- dphyl/max(dphyl)   #Make sure the distance are between 0 and 1 for the Jost correction
      temp_qdec<- Qdecomp(dphyl,t(sample), w=weight)   #Call the Qdecomp function for alpha, gamma and beta estimations.
      PD$Mean_Alpha = 1/(1-temp_qdec$Mean_alpha)
      PD$Alpha = 1/(1-temp_qdec$Alpha)
      PD$Gamma = 1/(1-temp_qdec$Gamma)
      PD$Beta_add = (PD$Gamma -PD$Mean_Alpha )
      PD$Beta_prop = 100*PD$Beta_add/PD$Gamma
      #Call the disc function for alpha, gamma and beta estimations for each pair of samples
      # PD$Pairwise_samples<- disc(as.data.frame(sample), as.dist(dphyl), structure=structure, Jost=Jost)
    }
    else {
      temp_qdec<- Qdecomp(dphyl,t(sample), w=weight)  #Call the Qdecomp function for alpha, gamma and beta estimations.
      PD$Mean_Alpha = temp_qdec$Mean_alpha
      PD$Alpha = temp_qdec$Alpha
      PD$Gamma = temp_qdec$Gamma
      PD$Beta_add = (PD$Gamma -PD$Mean_Alpha )
      PD$Beta_prop = 100*PD$Beta_add/PD$Gamma
      #PD$Beta =  temp_qdec$Beta
      #Call the disc function for alpha, gamma and beta estimations for each pair of samples
      # PD$Pairwise_samples <- disc(as.data.frame(sample), as.dist(dphyl), structure=structure, Jost=Jost)
    }
    
    
    
    
  }
  out <- list(TD, FD, PD)
  names(out) <- c("TD", "FD", "PD")
  return(out)
  
}



# Utility function to create partial dependence plot ------------------------------------------
make_pdp <- function(predictor.variable.names=predictor.variable.names,
                     model){
  
  desired_length <- length(predictor.variable.names)
  pdp.list <- vector(mode = "list", length = desired_length)
  names(pdp.list) = predictor.variable.names
  
  for(k in 1:desired_length){
    var = predictor.variable.names[k]
    
    dat.pdp <- pdp::partial(
      model,
      train = training,
      pred.var = var,
      plot = FALSE,
    )
    
    library(ggplotify)
    
    p_pdp <- ggplot(dat.pdp, aes(x=dat.pdp[,var], y=yhat)) + 
      ggtitle(dependent.variable.name) + 
      xlab(var) + 
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black")) + 
      ylab(dependent.variable.name) + 
      geom_point() +
      geom_smooth(method=loess, span=1, se=T,  col='black') +
      geom_smooth(method = lm, se = T)
    
    pdp.list[[k]] = list(dat.pdp, p_pdp)
  }
  return(pdp.list)
}



# Utility function to deregister parallel backend ------------------------------
unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
unregister_dopar()



