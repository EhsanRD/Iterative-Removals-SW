# install.packages("rsconnect")
# library("rsconnect")

library("shiny")
library("ggplot2")
library("reshape2")
library("plyr")
library("swCRTdesign")
library("matrixcalc")
library("scales")
library("tidyverse")
library("shinythemes")
library("Matrix")
library("plotly")
library("RColorBrewer")
library("gganimate")
library("gifski")
#library("cargo")
library("tidyr")
library("animation")
library("htmlwidgets")
library("plotly")
#library("manipulateWidget")
#library("htmlwidgets")
library("rsconnect")
require("gridExtra")

#setwd("~/Google Drive/Shared drives/Ehsan PhD work/Codes/")
#setwd("G:\\Shared drives\\Ehsan PhD work\\Codes\\Git\\Iterative-Removals-SW\\IterativeInfContent")
source("CRTVarAdj_func.R", local=TRUE)
#source("ICcell_appfunc.R", local=TRUE)

#function for generating design matrices
SWdesmat <- function(Tp) {
  Xsw <- matrix(data=0, ncol = Tp, nrow = (Tp-1))
  for(i in 1:(Tp-1)) {
    Xsw[i,(i+1):Tp] <- 1
  }
  return(Xsw)
}
#######

#assume one cluster is randomised to each sequence   

#function to calculate the information contents utilising centrosymmetric property 
ICPair = function(Xdes,m,rho0,r,type){
  
  Tp <- ncol(Xdes)
  K  <- nrow(Xdes)
  varD=CRTVarGeneralAdj(Xdes,m,rho0,r,type)
  ICmat<-matrix(data=NA, nrow=nrow(Xdes), ncol=ncol(Xdes))
  
    for(i in 1:nrow(Xdes)){
      for (j in 1:ncol(Xdes)){
        if(is.na(Xdes[i,j])==TRUE | is.na(Xdes[K-i+1,Tp-j+1])==TRUE){
          ICmat[i,j] <- NA
          ICmat[K-i+1,Tp-j+1] <- NA
        }
        else if(is.na(Xdes[i,j])==FALSE & is.na(Xdes[K-i+1,Tp-j+1])==FALSE) {
          Xdesij <- Xdes
          Xdesij[i,j] <- NA
          Xdesij[K-i+1,Tp-j+1] <- NA
          ICmat[i,j] <- CRTVarGeneralAdj(Xdesij,m,rho0,r,type)/varD
          ICmat[K-i+1,Tp-j+1] <- ICmat[i,j]
        }
        #avoid replicating loop
        if (i*j==(nrow(Xdes)*ncol(Xdes)/2)){
          break
        }
      }
    }
  return(ICmat)
}


IterRemove = function(Tp,m,rho0,r,type){
  K=Tp-1
  mval <- list()    #minimum lowest information content values
  Xdlist <- list()  #design matrix
  dlist <- list()   #information content matrix
  Xdlist[[1]] <- SWdesmat(Tp)
  dlist[[1]] <- ICPair(Xdlist[[1]],m,rho0,r,type)
  
  varvec<- c()
  varvec[1]<-CRTVarGeneralAdj(Xdlist[[1]],m,rho0,r,type)
  #removal of pairs of cells
  for (i in 2:(Tp*K/2-1)){ #most minimal design 
    mval <- which(dlist[[i-1]]==min(dlist[[i-1]],na.rm = TRUE), arr.ind = TRUE)
    
    # if (is.na(mval[1])) {
    #   dlist[[i-1]]<- NULL
    #   break
    #    }
    
    Xdlist[[i]]=Xdlist[[i-1]]
    #remove a pair of centrosymmetric cells
    #re-order indices by cluster and period
    mval <- mval[order(mval[,1],mval[,2]),]
    cellid <- mval[1,]
    clustid<-cellid[1]
    perid<-cellid[2]
    Xdlist[[i]][clustid,perid]<- NA
    Xdlist[[i]][K+1-clustid,Tp+1-perid]<- NA
    
    varvec[i]<-CRTVarGeneralAdj(Xdlist[[i]],m,rho0,r,type)#modify CRTVarGeneralAdj to give Na if variance cannt be calculated
    
    if (is.na(varvec[i])) {
    Xdlist[[i]]<- NULL
    break
    }
    
    dlist[[i]] = ICPair(Xdlist[[i]],m,rho0,r,type)
    #checking if dlist contains all NAs, if yes calculate the variance one more time
   }
  return(list(dlist,Xdlist, varvec))
}
