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
setwd("G:\\Shared drives\\Ehsan PhD work\\Codes\\Git\\Iterative-Removals-SW\\IterativeInfContent")
source("CRTVarAdj_func.R", local=TRUE)
#source("ICcell_appfunc.R", local=TRUE)
# Functions for generating design matrices
SWdesmat <- function(Tp) {
  Xsw <- matrix(data=0, ncol = Tp, nrow = (Tp-1))
  for(i in 1:(Tp-1)) {
    Xsw[i,(i+1):Tp] <- 1
  }
  return(Xsw)
}
#######
IterRemoval <- function(Tp,m, rho0, r, type,cpnum) {

        K=Tp-1
        Xdes <- SWdesmat(Tp)
        varmatall <- c()
        varmatall<- CRTVarGeneralAdj(Xdes,m,rho0,r,type)
        varmat_excl<-matrix(data=NA, nrow=nrow(Xdes), ncol=ncol(Xdes))
        
        
        #Need to change the way of coding
        IC_func2 = function(Xdes,varmatall){
          for(i in 1:nrow(Xdes)){
            for (j in 1:ncol(Xdes)){
              if(is.na(Xdes[i,j])==TRUE | is.na(Xdes[K-i+1,Tp-j+1])==TRUE){
                varmat_excl[i,j] <- NA
                varmat_excl[K-i+1,Tp-j+1] <- NA
              }
              else if(is.na(Xdes[i,j])==FALSE & is.na(Xdes[K-i+1,Tp-j+1])==FALSE) {
                Xdesij <- Xdes
                Xdesij[i,j] <- NA
                Xdesij[K-i+1,Tp-j+1] <- NA
                varmat_excl[i,j] <- CRTVarGeneralAdj(Xdesij,m,rho0,r,type)/varmatall
                varmat_excl[K-i+1,Tp-j+1] <- varmat_excl[i,j]
              }
            }
          }
          return(varmat_excl)
        }
        mval <- list()
        Xdlist <- list()
        dlist <- list()
        Xdlist[[1]] <- Xdes
        dlist [[1]]<- IC_func2(Xdlist[[1]],varmatall)
        varmatall[1] <- varmatall
        
        #remove all low-infomration content cells and updating.
        for (i in 2:((Tp*K)/1)){
        #lapply(2:((T*K)/2))
          mval[[i-1]] <- tryCatch(which(IC_func2(Xdlist[[i-1]],varmatall[i-1])==min(IC_func2(Xdlist[[i-1]],varmatall[i-1]),na.rm = TRUE), arr.ind = TRUE), warning=function(w) NA)
        
          if (is.na(mval[[i-1]][1])) {
            dlist[[i-1]]<- NULL
            break
          }
          
          Xdlist[[i]]=Xdlist[[i-1]]
          
            if (cpnum==0){
          # Loop is not good here, chenge it in future
          tryCatch(for (j in 1:dim(mval[[i-1]])[1]){
            Xdlist[[i]][mval[[i-1]][[j]],mval[[i-1]][[dim(mval[[i-1]])[1]+j]]]<- NA
            Xdlist[[i]][K+1-mval[[i-1]][[j]],Tp+1-mval[[i-1]][[dim(mval[[i-1]])[1]+j]]]<- NA
          }, error=function(e) NA)
            }
          #remove the smallest cluster and period, and removing the corresponding pair. 
          else if (cpnum==1){
                    mval[[i-1]] <- mval[[i-1]][order(mval[[i-1]][,1],mval[[i-1]][,2]),]
                    Xdlist[[i]][mval[[i-1]][[1]],mval[[i-1]][[dim(mval[[i-1]])[1]+1]]]<- NA
                    Xdlist[[i]][K+1-mval[[i-1]][[1]],Tp+1-mval[[i-1]][[dim(mval[[i-1]])[1]+1]]]<- NA
                  }
          varmatall[i] <-tryCatch(CRTVarGeneralAdj(Xdlist[[i]],m,rho0,r,type),error=function(e) NA)
          if (is.na(varmatall[i])) {
            Xdlist[[i]]<- NULL
            break
          }
          dlist[[i]] = IC_func2(Xdlist[[i]],varmatall[i])
        }
        return(list(dlist,Xdlist,varmatall))

}
