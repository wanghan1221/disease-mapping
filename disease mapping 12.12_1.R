#!/usr/bin/env Rscript
#Use ni,pi to generate data
#use two model fitting
#use Ohio dataset
#use real n and y
#tau_true=1
rm(list=ls(all=TRUE))
library(MASS)


dir0 ='C:/Users/lenovo/Desktop/diseasemappingresult/2013/data/'
dir ='C:/Users/lenovo/Desktop/diseasemappingresult/2013/result/'
##################################1.some useful functions
## for CAR model, use Gibbs to generate phi
phigen=function(tau_in,W,n)  
{
 #tau is the precise, W is the proximity matrix,n is the number of iteration
 #begin with all zero
 I=dim(W)[1]
 phi=rep(0,I)
 sn=round(n/3,0)
 res=matrix(0,nrow=sn,ncol=I)
 #using Gibbs for update
 for(t in 1:n)
 {
  for(i in 1:I)
  {
   thesum=sum(W[i,])
   theu=(W[i,]%*%phi)/thesum
   thevar=sqrt(1/(tau_in*thesum))
   phi[i]=rnorm(1,theu,thevar)
  }
  phi=phi-mean(phi)
  if(t>(n-sn))
  {
   res[t-n+sn,]=phi
  }
 }
 phi=apply(res,2,median)
 return(phi)
}



## calculate the posterior of phi using revised
post_phi1=function(thephi,i,W,phi_in,beta_in,tau_in)
{
  #thephi is the current value of phi[i], i the index of phi
  thesum=sum(W[i,])
  thesum
  phibar=((W[i,]%*%phi_in)/thesum)[1,1] #this value isn't influenced by the current value of phi[i]
  phibar
  part1=0.5*tau_in*thesum*(thephi-phibar)^2
  part1
  part2=N[i]*(1+exp((-1)*(x[i,]%*%beta_in+thephi)))^(-1)
  part2
  part3=y[i]*log(1+exp((-1)*(x[i,]%*%beta_in+thephi)))
  part3
  res=part1*(-1)-part2-part3
  return(res)
}




## update phi using MCMC W.mat,phi,beta,tau,sigma_in_phi,Racc_phi,Eused
update_phi=function(W,phi_in,beta_in,tau_in,sigma_in,Racc_phi_in,E_in)
{
  #index is 1 for revised model, 2 for standard model
  #phi_in,beta_in,tau_in are the current value
  #update each phi[i]
  #sigma is the variance for proposal dist
  for(i in 1:I)
  {
   #get proposal
   thecurr=phi_in[i]
   thestar=rnorm(1,thecurr,sigma_in[i]) #############i???????????###############
   logA=post_phi1(thestar,i,W,phi_in,beta_in,tau_in)
   logB=post_phi1(thecurr,i,W,phi_in,beta_in,tau_in)
   logrr=logA-logB
   logu=log(runif(1,0,1))
   if(logrr>logu)
   {
    phi_in[i]=thestar
    Racc_phi_in[i]=Racc_phi_in[i]+1
   }
  }
  #then centering phi
  phi_in=phi_in-mean(phi_in)
  return(list(phi_in,Racc_phi_in))
}



## calculate the posterior of beta using revised#############################################beta2
post_beta1=function(thebeta,i,phi_in,beta_in)
{
 #thebeta is the current value of beta[i], i the index of beta
 newbeta=beta_in
 newbeta[i]=thebeta    ######？？？？？？？？？？？？？？？？？############
 a=N/(1+exp((x%*%newbeta+phi_in)*(-1)))
 part1=sum(a)
 part2=sum(y*log(a))
 res=part2-part1
 return(res)
}



## update beta      phi,beta,sigma_in_beta,Racc_beta,Eused  
update_beta=function(phi_in,beta_in,sigma_in,Racc_beta_in,E_in)
{
  #index is 1 for revised model, 2 for standard model
  #phi_in,beta_in are the current value
  #update each beta[i]
  #sigma is the variance for proposal dist
  for(i in 1:length(beta_in))
  {
   #get proposal
   thecurr=beta_in[i]
   thestar=rnorm(1,thecurr,sigma_in[i])
   logA=post_beta1(thestar,i,phi_in,beta_in)
   logB=post_beta1(thecurr,i,phi_in,beta_in)
   logrr=logA-logB
   logu=log(runif(1,0,1))
   if(logrr>logu)
   {
    beta_in[i]=thestar
    Racc_beta_in[i]=Racc_beta_in[i]+1
   }
  }
  return(list(beta_in,Racc_beta_in))
}



## update tau 
update_tau=function(phi_in,a,b,W)
{
  #calculte(phi_i-phi_j)^2
  philist=rep(phi_in,I)
  p1=matrix(philist,nrow=I,ncol=I,byrow=T)
  p2=matrix(philist,nrow=I,ncol=I,byrow=F)
  pp=(p1-p2)^2
  ppp=pp[lower.tri(pp)]
  www=W[lower.tri(W)]
  newb=b+sum(ppp*www)/2
  newa=a+I/2
  newtau=rgamma(1,newa,newb)
  return(newtau)
}



###################################2.read spatial map, distance and N
W.mat=read.table(paste(dir0,"W.mat.txt",sep=""),header=F)
W.mat=as.matrix(W.mat)
####dist=read.table(paste(dir0,"distance.txt",sep=""),header=F)
#######dist=as.matrix(dist)
Ns=read.csv(paste(dir0,"Ns_1.csv",sep=""),header=T)
N=Ns[,3]
I=dim(W.mat)[1]
Ys=read.csv(paste(dir0,"Ys_1.csv",sep=""),header=T)
y=Ys[,3]
x0=rep(1,I)
X=x0
M=2
x=read.csv(paste(dir0,"x.csv",sep=""),header = F)

x<-as.matrix(x,nrow=I, ncol=M)
#x<-log(x)
## (3)init some parameters
#### some global parameters
a=b=1
n.T=100000
n.burn0=10000
n.burn=10000
n.tune=5000
times=10 #times for tune sigma
step=5 #save posterior every stepth,so the total number is (n.T-n.burn)/step
saveT=(n.T-n.burn)/step

 #### (c) init some saved matrix for posteriors
 ###### (c1) get the final used sigma for MH
 sigma0=0.3 #the inital value for sigma when choosing phi   #####自己试！！！
 sigma1=0.04  #the inital value for sigma when choosing beta    
 Rsigma_phi=rep(0,I)  #save sigma for each parameter
 Rsigma_beta=rep(0,M)
 Racc_phi=rep(0,I)  #count for acceptance
 Racc_beta=rep(0,M)

 ###### (c2)save posterior draws for revised model
 Rphisave=matrix(0,nrow=saveT,ncol=I)
 Rbetasave=matrix(0,nrow=saveT,ncol=M)
 Rtausave=rep(0,saveT)
 ####### (c21)using spatial fitting
 Rprobsave=matrix(0,nrow=saveT,ncol=I) #save for probability p[i]
 RSMRsave=matrix(0,nrow=saveT,ncol=I) #save for SMR
 REsave=matrix(0,nrow=saveT,ncol=I) #save for E

 
 
 #### (d)estimate using revised model
 #### the data we use are N,y,x1
 ###### (d1)initial parameters
 #Eused=rep(0,I) #we do not use Eused here
 beta=rnorm(M,0,1)
 tau=rgamma(1,a,b)
 phi=mvrnorm(I,0,1) 
 phi=phi[,1] 
 sigma_in_phi=rep(sigma0,I)
 sigma_in_beta=rep(sigma1,M)
 beta=as.matrix(beta,nrow=M,ncol=1)
 beta
 ###### (d2)start with burn-in
 for(nn in 1:n.burn0)
 {
  ######## (1)update each phi using M-H  ###### 此处用了Gibbs抽样的想法！依次更新参数代入下个个参数的后验分布中！
  theres=update_phi(W.mat,phi,beta,tau,sigma_in_phi,Racc_phi,Eused)
  phi=theres[[1]]
  Racc_phi=theres[[2]]
  ######## (2)update each beta using M-H
  theres=update_beta(phi,beta,sigma_in_beta,Racc_beta,Eused)
  beta=theres[[1]]
  Racc_beta=theres[[2]]
  ######## (3)update tau
  tau=update_tau(phi,a,b,W.mat)
 }
  
 ###### (d3)tune up parameter
 for(tt in 1:times)
 {
  # init acceptance count
  Racc_phi=rep(0,I)  #count for acceptance
  Racc_beta=rep(0,M)
  for(nn in 1:n.tune)
  {
   ######## (1)update each phi using MH
   theres=update_phi(W.mat,phi,beta,tau,sigma_in_phi,Racc_phi,Eused)
   theres
   phi=theres[[1]]
   Racc_phi=theres[[2]]
   ######## (2)update each beta using M-H
   theres=update_beta(phi,beta,sigma_in_beta,Racc_beta,Eused)
   beta=theres[[1]]
   Racc_beta=theres[[2]]
   ######## (3)update tau
   tau=update_tau(phi,a,b,W.mat)
  }
  # check whether the acc_rate statisfied   ####### 按规定MH算法中参数的接受率要在某个区间！！
  Racc_phi=Racc_phi/n.tune
  Racc_beta=Racc_beta/n.tune
  for(ii in 1:I)
  {
   if(Racc_phi[ii]>0.30)
   {
    sigma_in_phi[ii]=sigma_in_phi[ii]+0.05
   }
   if(Racc_phi[ii]<0.2)
   {
    sigma_in_phi[ii]=sigma_in_phi[ii]-0.05
   }
  }
  for(ii in 1:M)
  {
   if(Racc_beta[ii]>0.35)
   {
    sigma_in_beta[ii]=sigma_in_beta[ii]+0.008
   }
   if(Racc_beta[ii]<0.2)
   {
    sigma_in_beta[ii]=sigma_in_beta[ii]-0.008
   }
  }
 }
 Rsigma_phi=sigma_in_phi
 Rsigma_beta=sigma_in_beta
 
 
 ###### (d4)model fitting
 Racc_phi=rep(0,I)  #count for acceptance
 Racc_beta=rep(0,M)
 for(nn in 1:n.T)
 {
  ######## (1)update each phi using MH
  theres=update_phi(W.mat,phi,beta,tau,Rsigma_phi,Racc_phi,Eused)
  phi=theres[[1]]
  Racc_phi=theres[[2]]
  ######## (2)update each beta using M-H
  theres=update_beta(phi,beta,Rsigma_beta,Racc_beta,Eused)
  beta=theres[[1]]
  Racc_beta=theres[[2]]
  ######## (3)update tau
  tau=update_tau(phi,a,b,W.mat)
  ######## (4)save posteriors after burnin
  if((nn>n.burn)&((nn-n.burn)%%step==0))
  {
   Rphisave[(nn-n.burn)/step,]=phi
   Rbetasave[(nn-n.burn)/step,]=beta
   Rtausave[(nn-n.burn)/step]=tau  
  }
 }
 
 Racc_phi=Racc_phi/n.T
 Racc_beta=Racc_beta/n.T
 
 
 Racc_phi
 Racc_beta
 
 
 ######## (5) calculate E_est for standard model
 for(nn in 1:saveT)
 {
   ## for revised model using spatial
   aa=x%*%Rbetasave[nn,]+Rphisave[nn,] #linear prediction
   aa
   Rbetasave
   Rphisave
   p=exp(aa)/(1+exp(aa))
   p
   bb=sum(N*p)/sum(N)
   bb
   Rprobsave[nn,]=p
   Rprobsave
  
   RSMRsave[nn,]=p/bb
   RSMRsave
   REsave[nn,]=N*bb
   REsave
 }
 E_est=apply(REsave,2,mean)
 ################end estimation using revised model
 
 

 ################ (g)after model fitting, get SMR
 ## (g1)for revised, spatial
 #### (g12)for prob
 prob_revised_each_mean=apply(Rprobsave,2,mean)
 prob_revised_each_sd=apply(Rprobsave,2,sd)
 prob_revised_each_upper=upper=apply(Rprobsave,2,function(x){return(quantile(x,0.95))})
 prob_revised_each_lower=lower=apply(Rprobsave,2,function(x){return(quantile(x,0.05))})
 Rprob=cbind(prob_revised_each_mean,prob_revised_each_sd,
  prob_revised_each_upper,prob_revised_each_lower)

 #### (g13)for SMR
 SMR_revised_each_mean=apply(RSMRsave,2,mean)
 SMR_revised_each_sd=apply(RSMRsave,2,sd)
 SMR_revised_each_upper=upper=apply(RSMRsave,2,function(x){return(quantile(x,0.95))})
 SMR_revised_each_lower=lower=apply(RSMRsave,2,function(x){return(quantile(x,0.05))})
 RSMR=cbind(SMR_revised_each_mean,SMR_revised_each_sd,
  SMR_revised_each_upper,SMR_revised_each_lower)
 


 #### (g15)for beta
 beta_revised_each_mean=apply(Rbetasave,2,mean)
 beta_revised_each_sd=apply(Rbetasave,2,sd)
 beta_revised_each_upper=upper=apply(Rbetasave,2,function(x){return(quantile(x,0.95))})
 beta_revised_each_lower=lower=apply(Rbetasave,2,function(x){return(quantile(x,0.05))})
 Rbeta=cbind(beta_revised_each_mean,beta_revised_each_sd,
  beta_revised_each_upper,beta_revised_each_lower)
 #########################end g1


 


write.csv(Rprob,paste(dir,"Rprob.csv",sep=""))
write.csv(RSMR,paste(dir,"RSMR.csv",sep=""))
###write.csv(RE,paste(dir,"RE.csv",sep=""))
write.csv(Rbeta,paste(dir,"Rbeta.csv",sep=""))
write.csv(Rbetasave,paste(dir,"Rbetasave.csv",sep = ""))
'########write.csv(Rprob_2,paste(dir,"Rprob_2.csv",sep=""))
#########write.csv(RSMR_2,paste(dir,"RSMR_2.csv",sep=""))
#########write.csv(Sprob_2,paste(dir,"Sprob_2.csv",sep=""))
#########write.csv(SSMR_2,paste(dir,"SSMR_2.csv",sep=""))

RSMR=read.csv(paste(dir,"RSMR.csv",sep=""),header=T)
##SSMR=read.csv(paste(dir,"SSMR.csv",sep=""),header=T)
RSMR=RSMR[,-1]
##SSMR=SSMR[,-1]
Rlen=RSMR[,3]-RSMR[,4]
##Slen=SSMR[,3]-SSMR[,4]
mean(Rlen<Slen)

