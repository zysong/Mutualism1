initial.dir<-getwd()
setwd("D:\\ZYS\\Documents\\12Winter\\Codes\\sim01")

library("msm")
library("Matrix")
options(max.print = 10^7)

# set parameters and initial conditions
set.seed(1)
ra<-1 # animal growth rate 
rp<-1 # flowering plant growth rate due to pollination
rs<-0.01 # flowering plant growth rate due to selfing 
ma<-.001 # animal death rate
mp<-.0005 # plant death rate
h<-10 ## scaling parameter of P/A ratio 
k<-50 ## discrimination efficiency
mua<-10^(-2) # animal mutation rate
mup<-10^(-2) # plant mutation rate
Dv<-.5
generation<-10^4 ##Set the length the simulation
Na0<-10
Np0<-1000

##B&b: genotype and phenotype of discrete reward level
B0<-.05
Bmin<-0
Bmax<-1
bmax<-2
bmin<--2
binterval<-0.01
nb<-(bmax-bmin)/binterval+1
nB<-(Bmax-Bmin)/binterval+1
blist<-seq(bmin,bmax,by=binterval)
Blist<-seq(Bmin,Bmax,by=binterval)
sdb<-sqrt(Blist*Dv)

##discrete choosiness
c0<-0.05
cmax<-1
cmin<-0
cinterval<-0.01
nc<-(cmax-cmin)/cinterval+1
clist<-seq(cmin,cmax,by=cinterval)

##the plant saturation ratio
Hratio<-function(np, nac){
  np/(np+h*sum((1-clist)*nac))
}

## probability of chosen reward level (q(b,c))
bchosen<- function(distb, distB, H) {
  1-H+H*exp(k*(clist)%*%t(blist))/
    colSums(exp(k*(blist)%*%t(clist))*
    as.vector(distb%*%distB))
  ##distb%*%distB is the phenotype distribution in the whole populaiton
  ##row is b and column is c
}

##average chosen reward level
bmeanchosen<-function(distb, distB, qbc) {
  as.vector(qbc%*%(blist*distb)%*%distB)
}

##mutation
###continuous mutation
#mut.con<- function(nlist, mu) {
#  forward<-nlist[-1]*mu/2
#  backward<-nlist[-length(nlist)]*mu/2
#  mutant<-c(forward,0)+c(0,backward)-c(0,forward)-c(backward,0)
#  return(mutant)
#}

###discrete mutation
mut.dis<- function(nlist, mu) {
  forward<-rbinom(length(nlist)-1, nlist[-1], mu/2)
  backward<-rbinom(length(nlist)-1, nlist[-length(nlist)], mu/2)
  mutant<-c(forward,0)+c(0,backward)-c(0,forward)-c(backward,0)
  return(mutant)
}

##population growth functions

##animal
agrowth<- function(na, nac, bmean, H) {
  abirth<-ra*H*bmean*(1-clist)
  newborn<-floor(nac*abirth)+rbinom(nc, 1, nac*abirth-floor(nac*abirth))
  newborn<-newborn+mut.dis(nlist=newborn, mu=mua)
  adeath<-floor(na*ma*nac)+rbinom(nc,1,na*ma*nac-floor(na*ma*nac))
  nac.new<-newborn+nac-adeath
  return(nac.new)
}

##plant
pgrowth<- function(np, npB, distb, distc, qbc, H) {
  GB<-(1-H)/(1-sum(distc*clist))*distc%*%((1-clist)*qbc)%*%distb
  pbirth<-(rp*GB+rs*(1-GB))*(1-Blist)
  newborn<-floor(npB*pbirth)+rbinom(nB, 1, npB*pbirth-floor(npB*pbirth))
  newborn<-newborn+mut.dis(nlist=newborn, mu=mup)
  pdeath<-floor(np*mp*npB)+rbinom(nB,1,np*mp*npB-floor(np*mp*npB))
  npB.new<- newborn+npB-pdeath
  return(npB.new)
}

# simulations of population dynamics and evolutionary dynamics  
### b ~ truncated normal distribution N(B,sdb) for each genotype B
bdist<-matrix(rep(0,nb*nB), nrow=nb) ## row is b, column is B
for (i in 2:(nb-1)) {
  bdist[i,]<-ptnorm(rep(bmin + (i-.5)*binterval, nB), Blist, sdb, bmin, bmax) - 
    ptnorm(rep(bmin + (i-1.5)*binterval, nB), Blist, sdb, bmin, bmax)
}
bdist[1,]<-ptnorm(rep(bmin + 0.5*binterval, nB), Blist, sdb, bmin, bmax)
bdist[nb,]<-1 - ptnorm(rep(bmax - 0.5*binterval, nB), Blist, sdb, bmin, bmax)

### c
cdist<-rep(0, nc)
cdist[(c0-cmin)/cinterval+1]<-1
Nac<-round(Na0*cdist)
Nalist<-rep(0, generation)
Nalist[1]<-Na0  
cmat<-Matrix(0, nrow=generation, ncol=nc,
             dimnames = list(Time=1:generation, Choosiness=clist),
             sparse=TRUE)
cmat[1,]<-cdist

### B: genotype that determines the mean reward
Bdist<-rep(0,nB)
Bdist[(B0-Bmin)/binterval+1]<-1 ## the initial genotype frequency: f(B=B0)=1
Bmat<-Matrix(0, nrow=generation, ncol=nB,
             dimnames = list(Time=1:generation, Reward=Blist),
             sparse=TRUE)
Bmat[1,]<-Bdist
NpB<-round(Np0*Bdist)
Nplist<-rep(0, generation)
Nplist[1]<-Np0

### H
Hlist<-rep(0,generation)
Hlist[1]<-Hratio(Nplist[1], Nac)

###eco-evol dynamics step by step
for (i in 1:(generation-1)) {
  q<-bchosen(distb=bdist, distB=Bdist, Hlist[i])
  meanb<-bmeanchosen(distb=bdist, distB=Bdist, qbc=q)
  
  Nac.new<-agrowth(na=Nalist[i], nac=Nac, bmean=meanb, Hlist[i])
  Nalist[i+1]<-sum(Nac.new)
  cdist<-Nac.new/Nalist[i+1]
  
  NpB.new<-pgrowth(np=Nplist[i], npB=NpB, distb=bdist, distc=cdist, qbc=q, Hlist[i])
  Nplist[i+1]<-sum(NpB.new)
  Bdist<-as.vector(NpB.new/Nplist[i+1])
  
  Nac<-Nac.new
  NpB<-NpB.new
  Bmat[i+1,]<-Bdist
  cmat[i+1,]<-cdist
  Hlist[i+1]<-Hratio(Nplist[i+1], Nac)
  
  if (Nalist[i+1]*Nplist[i+1]==0)  break
}

Bsumm<-summary(Bmat)
csumm<-summary(cmat)

Bstat<-data.frame(Time = rownames(Bmat)[Bsumm$i],
                  Reward = colnames(Bmat)[Bsumm$j],
                  Dist = Bsumm$x)

cstat<-data.frame(Time = rownames(cmat)[csumm$i],
                  Choosiness = colnames(cmat)[csumm$j],
                  Dist = csumm$x)

sink(file="Bdist.dat", append=FALSE)
Bstat
sink()

sink(file="cdist.dat", append=FALSE)
cstat
sink()

attach(Bstat)
plot(Time, Reward, pch=".", cex=5, col=gray(1-Dist), xlab="time",ylab="reward")
detach(Bstat)

detach("package:msm")
detach("package:Matrix")
setwd(initial.dir)