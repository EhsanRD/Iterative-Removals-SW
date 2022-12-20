source("CRTVarAdj_func.R", local=TRUE)
#function for generating design matrices
SWdesmat <- function(Tp) {
  Xsw <- matrix(data=0, ncol = Tp, nrow = (Tp-1))
  for(i in 1:(Tp-1)) {
    Xsw[i,(i+1):Tp] <- 1
  }
  return(Xsw)
}

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
          ICmat[i,j] <- round(CRTVarGeneralAdj(Xdesij,m,rho0,r,type)/varD,10)
          ICmat[K-i+1,Tp-j+1] <- ICmat[i,j]
          
       if (is.na(ICmat[i,j])==TRUE & is.na(ICmat[K-i+1,Tp-j+1]==TRUE)) {
          ICmat[i,j] <-101.101
          ICmat[K-i+1,Tp-j+1] <- ICmat[i,j]
          }
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
    Xdlist[[i]]=Xdlist[[i-1]]
    #remove a pair of centrosymmetric cells
    #re-order indices by cluster and period
    mval <- mval[order(mval[,1],mval[,2]),]
    cellid <- mval[1,]
    clustid<-cellid[1]
    perid<-cellid[2]
    Xdlist[[i]][clustid,perid]<- NA
    Xdlist[[i]][K+1-clustid,Tp+1-perid]<- NA
    
    varvec[i]<-CRTVarGeneralAdj(Xdlist[[i]],m,rho0,r,type)
    
    dlist[[i]] = ICPair(Xdlist[[i]],m,rho0,r,type)

   }
  return(list(dlist,Xdlist, varvec))
}
