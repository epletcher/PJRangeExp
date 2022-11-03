sampleSigma<-function(Nlat,Npred,a,b){  ###This function will return variance. Need to sqrt to get SD
  n<-length(c(Nlat))
  tmp.a<-a+n/2
  tmp.b<-b+(0.5*sum((Nlat-Npred)^2))
  
  s2=rinvgamma(1,tmp.a,tmp.b)
  
}

sampleObS<-function(Nlat,N,a,b){  ###This function will return variance. Need to sqrt to get SD

  n<-length(c(Nlat))
  tmp.a<-a+n/2
  tmp.b<-b+(0.5*sum((N-Nlat)^2,na.rm=T))
  return(rinvgamma(1,tmp.a,tmp.b))
  
}


#####for t+1 Minv%*%diag(1/G)%*%Npred[t+1,]



sampleLatent<-function(Npred,Nlat,N,G,M,Minv,sig.o,sig.p,tmax) {  
  ###This function sequentially adds info from the different data sources to inform latent state
  #zero out
  s2<-sig.p^2
  o1<-sig.o^2
  Vi=0
  v=0
  
  #1. Contribution of Observations
    Vi=1/o1
    v=(N[t,]/o1)
  
  
    #2. Contribution of previous time step. Not included if we are in the first time step
  if(t>1){
    mu1<-M%*%(diag(G[t,])%*%Nlat[t-1,])
    Vi=Vi+1/s2
    v=v+mu1/s2
  }
  
    #2. Contribution of previous time step
    # should fix G to be matrix
  if(t<tmax){
    mu2<-Minv%*%(diag(1/G[t+1,])%*%Nlat[t+1,]) ####Here y is arranged year by pixel and X is year by pixel by covariate 
    Vi=Vi+1/s2
    v=v+mu2/s2
  }
  
 
  V<-1/Vi
  Nlat_t<-rnorm(length(v),c(V*v),sqrt(V))
  return(Nlat_t)
  
}  




###############Sampling Dispersal and 

UpdateBeta<-function(tmax,beta,Nlat,M,p){
Npred<-matrix(NA,tmax,p)
G<-matrix(NA,tmax,p)
for (t in 2:tmax){
G[t,]<-exp(beta[1]+beta[2]*Nlat[t-1,])
Npred[t,]<-M%*%(diag(G[t,])%*%Nlat[t-1,])



  }
return(list(Npred=Npred,G=G))

}


UpdateDispersal<-function(tmax,tau,Nlat,G,p,D){
  Npred<-matrix(NA,tmax,p)
  Mint<-exp(-(D/tau))
  M<-t(Mint/apply(Mint,1,sum))
   for (t in 2:tmax){
    Npred[t,]<-M%*%(diag(G[t,])%*%Nlat[t-1,])
    
  }
  return(list(Npred=Npred,M=M))
}







