# Functions for generating design matrices
SWdesmat <- function(T) {
  Xsw <- matrix(data=0, ncol = T, nrow = (T-1))
  for(i in 1:(T-1)) {
    Xsw[i,(i+1):T] <- 1
  }
  return(Xsw)
}
####################
CRTVarGeneralAdj <- function(Xmat, m, rho0, r, type) {
  totalvar <- 1
  sig2CP <- rho0*totalvar
  sig2E <- totalvar - sig2CP
  
  sig2 <- sig2E/m
  M1=Xmat[,colSums(!is.na(Xmat)) > 0]
  
if (identical(M1,Xmat)==TRUE){
    
    T <- ncol(Xmat)
    K <- nrow(Xmat)
    
    Xvec <- as.vector(t(Xmat))
    stackI <- matrix(rep(diag(1,T)), nrow=K*T, ncol=T, byrow=TRUE)
    Zmat <- cbind(stackI[!is.na(Xvec),], Xvec[!is.na(Xvec)])
    
    #Variance matrix for one cluster, with decay in correlation over time
    #Vi <- diag(sig2,T) + matrix(sig2CP,nrow=T, ncol=T)
    #Constant decay var if type==0
    if(type==0) {
      Vi <-diag(sig2 +(1-r)*sig2CP, T) + matrix(data=sig2CP*r, nrow=T, ncol=T)
    }
    #exponential decay structure
    if(type==1) {
      Vi <- diag(sig2,T) + sig2CP*(r^abs(matrix(1:T,nrow=T, ncol=T, byrow=FALSE) - matrix(1:T,nrow=T, ncol=T, byrow=TRUE)))
    }
}  
  #Need to create the Vi for the *original* design matrix
  #and then chop out the relevant row and column
  #T+length(t) is to allow for this.
  
    else {

    T <- ncol(M1)
    K <- nrow(Xmat)
    Xvec <- as.vector(t(M1))
    
    stackI <- matrix(rep(diag(1,T)), nrow=K*T, ncol=T, byrow=TRUE)
    Zmat <- cbind(stackI[!is.na(Xvec),], Xvec[!is.na(Xvec)])
    
    tp=c()
    tp=which(colSums(is.na(Xmat)) == nrow(Xmat))
    fullT<-T+length(tp)
    #Variance matrix for one cluster, with decay in correlation over time
    #Vi <- diag(sig2,T) + matrix(sig2CP,nrow=T, ncol=T)
    #Constant decay var if type==0
    if(type==0) { 
      Vi <-diag(sig2 +(1-r)*sig2CP, fullT) + matrix(data=sig2CP*r, nrow=fullT, ncol=fullT)
    }
    #exponential decay structure
    if(type==1) {
      Vi <- diag(sig2,fullT) + sig2CP*(r^abs(matrix(1:fullT,nrow=fullT, ncol=fullT, byrow=FALSE) - matrix(1:fullT,nrow=fullT, ncol=fullT, byrow=TRUE)))
      #Vi <- diag(sig2,T) + sig2CP*(r^abs(matrix(1:T,nrow=T, ncol=T, byrow=FALSE) - matrix(1:T,nrow=T, ncol=T, byrow=TRUE)))
    }
    #Variance matrix for all clusters
    #Vi <- Vi[-c(t[1:(length(t))]), -c(t[1:(length(t))])]
    Vi <- Vi[-tp,-tp]
}

  #Variance matrix for all clusters
  Vall <- kronecker(diag(1,K), Vi)
  Vall <- Vall[!is.na(Xvec),!is.na(Xvec)]
  
  #vartheta <- solve((t(Zmat)%*%solve(Vall)%*%Zmat))[ncol(Zmat),ncol(Zmat)]
  
  #xtry <- try(solve((t(Zmat)%*%solve(Vall)%*%Zmat)))
  #if('try-error' %in% class(xtry)) return(NA)
  
  #there will be problems if Zmat is not of full column rank
  if(rankMatrix(Zmat)[1] < ncol(Zmat)) return(NA)
  else return(solve((t(Zmat)%*%solve(Vall)%*%Zmat))[ncol(Zmat),ncol(Zmat)])
  
  #return(vartheta)
}

