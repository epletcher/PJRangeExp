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
X<-enviro.var[1:tmax,,4:5] # heat load and elev; this is where covars go (1 = tmean, 2 = tmean, 3 = vpdmax, 4 = heatload, 5 = elev)
amax<-3
bmax<-3 # number of covariates 


####priors####


###Starting Values###
Nlat<-N #Starting values for latent states is the observed data
alpha0<-.019 ###Give beta some starting values based on what we know
alpha1<-0 # regression coef for heatload
alpha2<-0 # regresstion coef for elev
beta0<--0.001
beta1<-0 # regression coef for heatload
beta2<-0 # regresstion coef for elev
tau<-.033###Give tau a reasonable starting value. 
sig.p<-1.3##give sig.p reasonable starting values
o1<-sig.o<-1##give sig.o reasonable starting values # could try 1.5 or 2
ro <- 0.5
qo1 <- (ro/o1)+1


Mint<-exp(-(D/tau))
M<-t(Mint/apply(Mint,1,sum))  ##calculate M starting M given Tau

Npred<-G<-matrix(NA,tmax,pmax)

for (t in 2:tmax){
  G[t,]<-exp((alpha0+X[t,,1]*alpha1+X[t,,2]*alpha2)+(beta0+X[t,,1]*beta1+X[t,,2]*beta2)*Nlat[t-1,])
  Npred[t,]<-M%*%(diag(G[t,])%*%Nlat[t-1,])
}


Niter<-20000 ###Number of interations. Keep in mind this will need to be more than you needed for stan  
checkpoint=Niter*0.01
###Containers####
tauOut<-matrix(NA,Niter,)
alphaOut <- matrix(NA,Niter,amax) 
betaOut<-matrix(NA,Niter,bmax)
NlatOut<-array(NA,c(tmax,pmax,Niter/10)) # change to all pixels, but only every 10th iteration
NlatOutLast<-matrix(NA,pmax,Niter)
#rep.pix <- c(115:145, 910:940, 1865:1895) # representative pixels (high,med,low density)
tenIter <- seq(10,20000, by = 10) # vector of every 10th iteration
sig.pOut<-sig.oOut<-matrix(NA,Niter,1)

accept.alpha2=accept.alpha1=accept.alpha0=accept.beta2=accept.beta1=accept.beta0=accept.tau=0
#beta.tune=diag(c(.000001,.000001))
alpha0.tune=.0001
alpha1.tune=.001
alpha2.tune=.001
beta0.tune=.0001
beta1.tune=.001
beta2.tune=.001
tau.tune=.001

for (i in 1:Niter){
  
  #alpha 0
  alpha0.star=rnorm(1,alpha0,alpha0.tune)
  Out=UpdateBetaTop(tmax=tmax,a0=alpha0.star,a1=alpha1,a2=alpha2,X=X,b0=beta0,b1=beta1,b2=beta2,Nlat=Nlat,M=M,p=p)
  Npred.star<-Out$Npred
  G.star<-Out$G
  now=UpdateBetaTop(tmax=tmax,a0=alpha0,a1=alpha1,a2=alpha2,X=X,b0=beta0,b1=beta1,b2=beta2,Nlat=Nlat,M=M,p=p)
  Npred<-now$Npred
  mh1=sum(dnorm(Nlat[-1,],(Npred.star[-1,]),sig.p,log=TRUE)) #implied uniform prior
  mh2=sum(dnorm(Nlat[-1,],(Npred[-1,]),sig.p,log=TRUE))      #implied uniform prior
  mh=min(exp(mh1-mh2),1)
  if(mh>runif(1)){
    G=G.star
    alpha0=alpha0.star
    accept.alpha0=accept.alpha0+1
    
  }
  alphaOut[i,1]<-alpha0
  
  #alpha1
  alpha1.star=rnorm(1,alpha1,alpha1.tune)
  Out=UpdateBetaTop(tmax=tmax,a0=alpha0,a1=alpha1.star,a2=alpha2,X=X,b0=beta0,b1=beta1,b2=beta2,Nlat=Nlat,M=M,p=p)
  Npred.star<-Out$Npred
  G.star<-Out$G
  now=UpdateBetaTop(tmax=tmax,a0=alpha0,a1=alpha1,a2=alpha2,X=X,b0=beta0,b1=beta1,b2=beta2,Nlat=Nlat,M=M,p=p)
  Npred<-now$Npred
  mh1=sum(dnorm(Nlat[-1,],(Npred.star[-1,]),sig.p,log=TRUE)) #implied uniform prior
  mh2=sum(dnorm(Nlat[-1,],(Npred[-1,]),sig.p,log=TRUE))      #implied uniform prior
  mh=min(exp(mh1-mh2),1)
  if(mh>runif(1)){
    G=G.star
    alpha1=alpha1.star
    accept.alpha1=accept.alpha1+1
    
  }
  alphaOut[i,2]<-alpha1
  
  #alpha2
  alpha2.star=rnorm(1,alpha2,alpha2.tune)
  Out=UpdateBetaTop(tmax=tmax,a0=alpha0,a1=alpha1,a2=alpha2.star,X=X,b0=beta0,b1=beta1,b2=beta2,Nlat=Nlat,M=M,p=p)
  Npred.star<-Out$Npred
  G.star<-Out$G
  now=UpdateBetaTop(tmax=tmax,a0=alpha0,a1=alpha1,a2=alpha2,X=X,b0=beta0,b1=beta1,b2=beta2,Nlat=Nlat,M=M,p=p)
  Npred<-now$Npred
  mh1=sum(dnorm(Nlat[-1,],(Npred.star[-1,]),sig.p,log=TRUE)) #implied uniform prior
  mh2=sum(dnorm(Nlat[-1,],(Npred[-1,]),sig.p,log=TRUE))      #implied uniform prior
  mh=min(exp(mh1-mh2),1)
  if(mh>runif(1)){
    G=G.star
    alpha2=alpha2.star
    accept.alpha2=accept.alpha2+1
    
  }
  alphaOut[i,3]<-alpha2
  
  #beta 0
  beta0.star=rnorm(1,beta0,beta0.tune)
  Out=UpdateBetaTop(tmax=tmax,a0=alpha0,a1=alpha1,a2=alpha2,X=X,b0=beta0.star,b1=beta1,b2=beta2,Nlat=Nlat,M=M,p=p)
  Npred.star<-Out$Npred
  G.star<-Out$G
  now=UpdateBetaTop(tmax=tmax,a0=alpha0,a1=alpha1,a2=alpha2,X=X,b0=beta0,b1=beta1,b2=beta2,Nlat=Nlat,M=M,p=p)
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
  Out=UpdateBetaTop(tmax=tmax,a0=alpha0,a1=alpha1,a2=alpha2,X=X,b0=beta0,b1=beta1.star,b2=beta2,Nlat=Nlat,M=M,p=p)
  Npred.star<-Out$Npred
  G.star<-Out$G
  now=UpdateBetaTop(tmax=tmax,a0=alpha0,a1=alpha1,a2=alpha2,X=X,b0=beta0,b1=beta1,b2=beta2,Nlat=Nlat,M=M,p=p)
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
  
  #beta 2
  beta2.star=rnorm(1,beta2,beta2.tune)
  Out=UpdateBetaTop(tmax=tmax,a0=alpha0,a1=alpha1,a2=alpha2,X=X,b0=beta0,b1=beta1,b2=beta2.star,Nlat=Nlat,M=M,p=p)
  Npred.star<-Out$Npred
  G.star<-Out$G
  now=UpdateBetaTop(tmax=tmax,a0=alpha0,a1=alpha1,a2=alpha2,X=X,b0=beta0,b1=beta1,b2=beta2,Nlat=Nlat,M=M,p=p)
  Npred<-now$Npred
  mh1=sum(dnorm(Nlat[-1,],(Npred.star[-1,]),sig.p,log=TRUE)) #implied uniform prior
  mh2=sum(dnorm(Nlat[-1,],(Npred[-1,]),sig.p,log=TRUE))      #implied uniform prior
  mh=min(exp(mh1-mh2),1)
  if(mh>runif(1)){
    G=G.star
    beta2=beta2.star
    accept.beta2=accept.beta2+1
    
  }
  betaOut[i,3]<-beta2
  
  # dispersal param
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
  
  
  
  for (t in 2:tmax){
    Npred[t,]<-M%*%(diag(G[t,])%*%Nlat[t-1,])
  }
  sig.p<-sqrt(sampleSigma(Nlat=c(Nlat[-1,]),Npred=c(Npred[-1,]),a=3,b=.5))
  sig.pOut[i,]<-sig.p
  
  #sig.o<-sqrt(sampleObS(Nlat=c(Nlat),N=c(N),a=3,b=.5))
  sig.oOut[i,]<-sig.o
  
  
  for (t in 1:tmax){
    Nlat[t,]<-sampleLatent(Npred,Nlat,N,G,M,Minv,sig.o,sig.p,tmax)
  }
  if(i %in% tenIter) {
    NlatOut[,,i/10] <- Nlat}
  
  NlatOutLast[,i]<-Nlat[tmax,]
  
  print(i)
  
  if(i%%checkpoint==0){
    if(accept.alpha0/i<0.35) alpha0.tune=alpha0.tune*.9
    if(accept.alpha0/i>0.45) alpha0.tune=alpha0.tune*1.1
    
    if(accept.alpha1/i<0.35) alpha1.tune=alpha1.tune*.9
    if(accept.alpha1/i>0.45) alpha1.tune=alpha1.tune*1.1
    
    if(accept.alpha2/i<0.35) alpha2.tune=alpha2.tune*.9
    if(accept.alpha2/i>0.45) alpha2.tune=alpha2.tune*1.1
    
    if(accept.beta0/i<0.35) beta0.tune=beta0.tune*.9
    if(accept.beta0/i>0.45) beta0.tune=beta0.tune*1.1
    
    if(accept.beta1/i<0.35) beta1.tune=beta1.tune*.9
    if(accept.beta1/i>0.45) beta1.tune=beta1.tune*1.1
    
    if(accept.beta2/i<0.35) beta2.tune=beta2.tune*.9
    if(accept.beta2/i>0.45) beta2.tune=beta2.tune*1.1
    
    if(accept.tau/i<0.35) tau.tune=tau.tune*.9
    if(accept.tau/i>0.45) tau.tune=tau.tune*1.1
  }
  
  if(i %in% seq(1000,20000, by = 1000)) {
    save.image(file = "R:/Shriver_Lab/PJspread/sampleroutput/sampler_topo_v1.RData")
  }
  
}

#save.image(file = "R:/Shriver_Lab/PJspread/sampleroutput/sampler_topo_v1.RData")


