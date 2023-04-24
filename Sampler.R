#.libPaths("C:/Rpackages/R/win-library/4.1") # (elise setting package library location)

# load sampler functions script
source("R:/Shriver_Lab/PJspread/PJ_spread_repo/SamplerFunctions.R")

# load 'dataprepped' workspace

library("splus2R")
library('LaplacesDemon')
###Data####
# N # observed data, assumed to be a matrix that is year by pixel
Noos <- N[31:36,] # out of sample data
tmax<-dim(N)[1]-5 # leave off last 5 years so that we can evaluate out of sample predictions
pmax<-dim(N)[2]
D<-Dsq
#X<-if you have covariates this is where they go
bmax<-2 #length(X[1,]) number of covariates

####priors####


###Starting Values###
Nlat<-N[1:31,] #Starting values for latent states is the observed data
alpha0<-0.0197 #Give tarting values based on what we know
beta0<--0.001
tau<-0.0332###Give tau a reasonable starting value. 
sig.p<-1.335##give sig.p reasonable starting values
o1<-sig.o<-1##give sig.o reasonable starting values
ro <- 0.5
qo1 <- (ro/o1)+1

Mint<-exp(-(D/tau))
M<-t(Mint/apply(Mint,1,sum))  ##calculate M starting M given Tau


Npred<-G<-matrix(NA,tmax,pmax)

for (t in 2:tmax){
  G[t,]<-exp(alpha0+beta0*Nlat[t-1,])
  Npred[t,]<-M%*%(diag(G[t,])%*%Nlat[t-1,])
}


Niter<-50000 ###Number of iterations. Base this off initial runs: still moving around until 12,000, upping from 20000 to 50000
burnin<-Niter*0.5
checkpoint=Niter*0.01

###Containers####
tauOut<-betaOut<-alphaOut<-matrix(NA,Niter,)
NlatOut<-array(NA,c(tmax,pmax,Niter/10)) # change to all pixels, but only every 10th iteration
NlatOutLast<-matrix(NA,pmax,Niter)
#rep.pix <- c(115:145, 910:940, 1865:1895) # representative pixels (high,med,low density)
tenIter <- seq(10,Niter, by = 10) # vector of every 10th iteration
sig.pOut<-sig.oOut<-matrix(NA,Niter,1)

# out of sample prediction evaluation
rmseTotOut<-biasOut<-denseOut<-matrix(NA,5,burnin) # evaluation metrics
Npredoos <- matrix(NA,5,pmax) # out of sample predictions

accept.beta0=accept.alpha0=accept.tau=0
#beta.tune=diag(c(.000001,.000001))
alpha0.tune=.000001
beta0.tune=.0001
tau.tune=.001


for (i in 1:Niter){ # edit starting iteration if start/stopping
  
  alpha0.star=rnorm(1,alpha0,alpha0.tune)
  Out=UpdateBeta(tmax=tmax,a0=alpha0.star,b0=beta0,Nlat=Nlat,M=M,p=p)
  Npred.star<-Out$Npred
  G.star<-Out$G
  now=UpdateBeta(tmax=tmax,a0=alpha0,b0=beta0,Nlat=Nlat,M=M,p=p)
  Npred<-now$Npred
  mh1=sum(dnorm(Nlat[-1,],(Npred.star[-1,]),sig.p,log=TRUE)) #implied uniform prior
  mh2=sum(dnorm(Nlat[-1,],(Npred[-1,]),sig.p,log=TRUE))      #implied uniform prior
  mh=min(exp(mh1-mh2),1)
  if(mh>runif(1)){
    G=G.star
    alpha0=alpha0.star
    accept.alpha0=accept.alpha0+1
    
  }
  alphaOut[i,]<-alpha0
  
  beta0.star=rnorm(1,beta0,beta0.tune)
  Out=UpdateBeta(tmax=tmax,a0=alpha0,b0=beta0.star,Nlat=Nlat,M=M,p=p)
  Npred.star<-Out$Npred
  G.star<-Out$G
  now=UpdateBeta(tmax=tmax,a0=alpha0,b0=beta0,Nlat=Nlat,M=M,p=p)
  Npred<-now$Npred
  mh1=sum(dnorm(Nlat[-1,],(Npred.star[-1,]),sig.p,log=TRUE)) #implied uniform prior
  mh2=sum(dnorm(Nlat[-1,],(Npred[-1,]),sig.p,log=TRUE))      #implied uniform prior
  mh=min(exp(mh1-mh2),1)
  if(mh>runif(1)){
    G=G.star
    beta0=beta0.star
    accept.beta0=accept.beta0+1
    
  }
  betaOut[i,]<-beta0
  
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
    
    if(accept.beta0/i<0.35) beta0.tune=beta0.tune*.9
    if(accept.beta0/i>0.45) beta0.tune=beta0.tune*1.1
    
    if(accept.tau/i<0.35) tau.tune=tau.tune*.9
    if(accept.tau/i>0.45) tau.tune=tau.tune*1.1
  }
  
  ## out of sample prediction
  if(i %in% seq(burnin+1,Niter,1)) {
    
    # withheld years
    Nt <- Nlat[31,] # set initial cover value as actual latent value
    
    for (t in 1:5){
      
      Gnew<-exp(alpha0+beta0*Nt)
      
      Nmean <-M%*%(diag(Gnew)%*%Nt)
      
      Nt <- rnorm(pmax, Nmean, sig.p)
      
      Npredoos[t,] <- Nt
      
    }
    
    # prediction evaluation metrics
    for(t in 1:5) {
      rmseTotOut[t,i] <- rmsefunc(pred=Npredoos[t,], obs=Noos[t,]) # cumulative rmse
      biasOut[t,i] <- biasfunc(pred=Npredoos[t,], obs=Noos[t,]) # bias
      denseOut[t,i] <- densefunc(pred=Npredoos[t,], obs=Noos[t,], sig_o=sig.o) # density
      }
    }
    
  # save
  if(i %in% seq(1000,Niter, by = 1000)) {
    save.image(file = "R:/Shriver_Lab/PJspread/sampleroutput/sampler_base_v4_c1.RData")
  } 
}



