#.libPaths("C:/Rpackages/R/win-library/4.1") # (elise setting package library location)

# load workspacefile with prepped data: 
# load("G:/.shortcut-targets-by-id/1FPlPAVacVgAROSPXMiiOGb2Takzm2241/PJ_Photo/cover_spread/Scripts/data_prepped.RData")

library("splus2R")
library('LaplacesDemon')

###Data####
N<-N[1:31,] # observed data, assumed to be a matrix that is year by pixel (remove last 5 years here, to test forecast)
tmax<-dim(N)[1] 
pmax<-dim(N)[2]
D<-Dsq
X<-enviro.var[1:tmax,,1:2] # covariates this is where they go (1 = ppt, 2 = tmean, 3 = vpdmax, 4 = heatload, 5 = elev)
bmax<-2
gmax<-2 # number of climate covariates (minus topographic vars)


####priors####


###Starting Values###
Nlat<-N #Starting values for latent states is the observed data
beta0<-.019 ###Give beta some starting values based on what we know
beta1<--0.001
tau<-.033###Give tau a reasonable starting value. 
gam0<-0 # regression coef for ppt
gam1<-0 # regresstion coef for tmean
sig.p<-1.3##give sig.p reasonable starting values
o1<-sig.o<-1##give sig.o reasonable starting values
ro <- 0.5
qo1 <- (ro/o1)+1


Mint<-exp(-(D/tau))
M<-t(Mint/apply(Mint,1,sum))  ##calculate M starting M given Tau

Npred<-G<-matrix(NA,tmax,pmax)

for (t in 2:tmax){
  G[t,]<-exp((beta0+X[t,,1]*gam0+X[t,,2]*gam1)+beta1*Nlat[t-1,])
  Npred[t,]<-M%*%(diag(G[t,])%*%Nlat[t-1,])
}


Niter<-20000 ###Number of interations. Keep in mind this will need to be more than you needed for stan  
checkpoint=Niter*0.01
###Containers####
tauOut<-matrix(NA,Niter,)
betaOut<-matrix(NA,Niter,bmax)
gammaOut <- matrix(NA,Niter,gmax) 
NlatOut<-array(NA,c(tmax,pmax,Niter/10)) # change to all pixels, but only every 10th iteration
NlatOutLast<-matrix(NA,pmax,Niter)
#rep.pix <- c(115:145, 910:940, 1865:1895) # representative pixels (high,med,low density)
tenIter <- seq(10,20000, by = 10) # vector of every 10th iteration
sig.pOut<-sig.oOut<-matrix(NA,Niter,1)

accept.gam1=accept.gam0=accept.beta1=accept.beta0=accept.tau=0
#beta.tune=diag(c(.000001,.000001))
beta0.tune=.0001
beta1.tune=.0001
tau.tune=.001
gam0.tune=.001
gam1.tune=.001

for (i in 1:Niter){
  
  #beta 0
  beta0.star=rnorm(1,beta0,beta0.tune)
  Out=UpdateBetaClim(tmax=tmax,b0=beta0.star,b1=beta1,X=X,gam0=gam0,gam1=gam1,Nlat=Nlat,M=M,p=p)
  Npred.star<-Out$Npred
  G.star<-Out$G
  now=UpdateBetaClim(tmax=tmax,b0=beta0,b1=beta1,X=X,gam0=gam0,gam1=gam1,Nlat=Nlat,M=M,p=p)
  Npred<-now$Npred
  mh1=sum(dnorm(Nlat[-1,],(Npred.star[-1,]),sig.p,log=TRUE)) #implied uniform prior
  mh2=sum(dnorm(Nlat[-1,],(Npred[-1,]),sig.p,log=TRUE))      #implied uniform prior
  mh=min(exp(mh1-mh2),1)
  if(mh>runif(1)){
    G=G.star
    beta0=beta0.star
    accept.beta0=accept.beta0+1
    
  }
  betaOut[i,1]<-beta0
  #beta 1
  beta1.star=rnorm(1,beta1,beta1.tune)
  Out=UpdateBetaClim(tmax=tmax,b0=beta0,b1=beta1.star,X=X,gam0=gam0,gam1=gam1,Nlat=Nlat,M=M,p=p)
  Npred.star<-Out$Npred
  G.star<-Out$G
  now=UpdateBetaClim(tmax=tmax,b0=beta0,b1=beta1,X=X,gam0=gam0,gam1=gam1,Nlat=Nlat,M=M,p=p)
  Npred<-now$Npred
  mh1=sum(dnorm(Nlat[-1,],(Npred.star[-1,]),sig.p,log=TRUE)) #implied uniform prior
  mh2=sum(dnorm(Nlat[-1,],(Npred[-1,]),sig.p,log=TRUE))      #implied uniform prior
  mh=min(exp(mh1-mh2),1)
  if(mh>runif(1)){
    G=G.star
    beta1=beta1.star
    accept.beta1=accept.beta1+1
    
  }
  betaOut[i,2]<-beta1
  # gamma0
  gam0.star=rnorm(1,gam0,gam0.tune)
  Out=UpdateBetaClim(tmax=tmax,b0=beta0,b1=beta1,X=X,gam0=gam0.star,gam1=gam1,Nlat=Nlat,M=M,p=p)
  Npred.star<-Out$Npred
  G.star<-Out$G
  now=UpdateBetaClim(tmax=tmax,b0=beta0,b1=beta1,X=X,gam0=gam0,gam1=gam1,Nlat=Nlat,M=M,p=p)
  Npred<-now$Npred
  mh1=sum(dnorm(Nlat[-1,],(Npred.star[-1,]),sig.p,log=TRUE)) #implied uniform prior
  mh2=sum(dnorm(Nlat[-1,],(Npred[-1,]),sig.p,log=TRUE))      #implied uniform prior
  mh=min(exp(mh1-mh2),1)
  if(mh>runif(1)){
    G=G.star
    gam0=gam0.star
    accept.gam0=accept.gam0+1
    
  }
  gammaOut[i,1]<-gam0
  #gamma1
  gam1.star=rnorm(1,gam1,gam1.tune)
  Out=UpdateBetaClim(tmax=tmax,b0=beta0,b1=beta1,X=X,gam0=gam0,gam1=gam1.star,Nlat=Nlat,M=M,p=p)
  Npred.star<-Out$Npred
  G.star<-Out$G
  now=UpdateBetaClim(tmax=tmax,b0=beta0,b1=beta1,X=X,gam0=gam0,gam1=gam1,Nlat=Nlat,M=M,p=p)
  Npred<-now$Npred
  mh1=sum(dnorm(Nlat[-1,],(Npred.star[-1,]),sig.p,log=TRUE)) #implied uniform prior
  mh2=sum(dnorm(Nlat[-1,],(Npred[-1,]),sig.p,log=TRUE))      #implied uniform prior
  mh=min(exp(mh1-mh2),1)
  if(mh>runif(1)){
    G=G.star
    gam1=gam1.star
    accept.gam1=accept.gam1+1
    
  }
  gammaOut[i,2]<-gam1
  #tau
  tau.star=rnorm(1,tau,tau.tune)
  Out=UpdateDispersal(tmax=tmax,tau=tau.star,Nlat=Nlat,G=G,p=p,D=D)
  Npred.star<-Out$Npred
  M.star<-Out$M
  now=UpdateDispersal(tmax=tmax,tau=tau,Nlat=Nlat,G=G,p=p,D=D)
  Npred<-now$Npred
  mh1=sum(dnorm(Nlat[-1,],(Npred.star[-1,]),sig.p,log=TRUE)) #implied uniform prior
  mh2=sum(dnorm(Nlat[-1,],(Npred[-1,]),sig.p,log=TRUE))      #implied uniform prior
  mh=min(exp(mh1-mh2),1)
  if(mh>runif(1)){
    M=M.star
    Minv<-chol2inv(M.star)
    tau<-tau.star
    accept.tau=accept.tau+1
    
  }
  tauOut[i,]<-tau
  
  
  #sigmas
  for (t in 2:tmax){
    Npred[t,]<-M%*%(diag(G[t,])%*%Nlat[t-1,])
  }
  sig.p<-sqrt(sampleSigma(Nlat=c(Nlat[-1,]),Npred=c(Npred[-1,]),a=3,b=.5))
  sig.pOut[i,]<-sig.p
  
  #sig.o<-sqrt(sampleObS(Nlat=c(Nlat),N=c(N),a=3,b=.5))
  sig.oOut[i,]<-sig.o
  
  #latent
  for (t in 1:tmax){
    Nlat[t,]<-sampleLatent(Npred,Nlat,N,G,M,Minv,sig.o,sig.p,tmax)
  }
  if(i %in% tenIter) {
    NlatOut[,,i/10] <- Nlat}
  
  NlatOutLast[,i]<-Nlat[tmax,]
  
  print(i)
  #tuning
  if(i%%checkpoint==0){
    if(accept.beta0/i<0.35) beta0.tune=beta0.tune*.9
    if(accept.beta0/i>0.45) beta0.tune=beta0.tune*1.1
    
    if(accept.beta1/i<0.35) beta1.tune=beta1.tune*.9
    if(accept.beta1/i>0.45) beta1.tune=beta1.tune*1.1
    
    if(accept.gam0/i<0.35) gam0.tune=gam0.tune*.9
    if(accept.gam0/i>0.45) gam0.tune=gam0.tune*1.1
    
    if(accept.gam1/i<0.35) gam1.tune=gam1.tune*.9
    if(accept.gam1/i>0.45) gam1.tune=gam1.tune*1.1
    
    if(accept.tau/i<0.35) tau.tune=tau.tune*.9
    if(accept.tau/i>0.45) tau.tune=tau.tune*1.1
  }
  
  
  
}

save.image(file = "R:/Shriver_Lab/PJspread/sampleroutput/sampler_clim_v2.RData")


