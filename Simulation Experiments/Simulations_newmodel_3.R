# Simulation results - Exp III - set1

require(ggplot2)
require(rmutil)
library(gridExtra)
library(snowfall)

#setwd('D:/Trambak/Side Information/simulations')
source('sideinfo_lib.R')

Ns<-seq(500,5000,100)
m<-50
#n<-5000
repts<- 1000

#------ set 1 -------------------------------
risk.u<- matrix(0,length(Ns),1)
risk.uv.1<- risk.u
risk.uv.2<- risk.u
risk.or<- risk.u
risk.ebt<- risk.u
risk.ejs<- risk.u

sfInit(parallel=TRUE,cpu=6)
sfSource('sideinfo_lib.R')
sfLibrary(rmutil)

wrapper<- function(i){
  n<-Ns[i]
  size.1<-n/100
  size.2<- 4*n/100
  size<- size.1+size.2
  mu.mat<- matrix(0,repts,n)
  eta<- matrix(0,repts,n)
  eta.1<- matrix(0,repts,n)
  eta.2<- matrix(0,repts,n)
  theta<-matrix(0,repts,n)
  q<- n^{-0.5}
  
  #Generate sparse mu
  for(reps in 1:repts){
    set.seed(reps)
    mu.mat[reps,]<- c(runif(size.1,6,7),runif(size.2,2,3),rep(0,n-size))
  }
  #Generate sparse perturbation
  for(reps in 1:repts){
    set.seed(reps^2)
    eta[reps,]<- (runif(n)<=q)*rnorm(n,2,0.1)+0*(runif(n)>q)
  }
  
  #Generate theta
  theta<- mu.mat+eta
  
  #Generate U_i~N(theta_i,sigma^2)
  var.u<- rep(1,n)
  U<- matrix(0,repts,n)
  for(reps in 1:repts){
    for(j in 1:n){
      set.seed(reps+j)
      U[reps,j]=rnorm(1,theta[reps,j],sqrt(var.u[j]))
    }
  }
  
  # Oracle threshold by minimizing true risk
  ttau<- seq(min(mu.mat),max(mu.mat),length.out = 50)
  temprisk.or<-matrix(0,repts,6)
  temp<- matrix(0,repts,length(ttau))
  for(reps in 1:repts){
    
    mu<- mu.mat[reps,]
    VV<-abs(mu)
    muhat.or<-matrix(0,n,1)
    for(k in 1:length(ttau)){
      tau<- ttau[k]
      i1<-(VV<=tau)
      i2<- (VV>tau)
      U1<- U[reps,i1]
      sigma21<-var.u[i1]
      U2<- U[reps,i2]
      sigma22<- var.u[i2]
      a1<- sureshrink(U1,sigma21)
      a2<- sureshrink(U2,sigma22)
      muhat.or[i1]<-a1
      muhat.or[i2]<-a2
      temp[reps,k]<- sum((theta[reps,]-muhat.or)^2)
    }
    index<- max(which(temp[reps,]==min(temp[reps,])))
    tau.min<- ttau[index]
    i1<-(VV<=tau.min)
    i2<- (VV>tau.min)
    U1<- U[reps,i1]
    sigma21<-var.u[i1]
    U2<- U[reps,i2]
    sigma22<- var.u[i2]
    a1<- sureshrink.mse(U1,sigma21,1,0)
    a2<- sureshrink.mse(U2,sigma22,1,0)
    temprisk.or[reps,]<- c(temp[reps,index],
                             a1$t,a2$t,tau.min,length(sigma21),length(sigma22))
  }
  or<- c(mean(temprisk.or[,1])/n,sd(temprisk.or[,1])/sqrt(repts))
  t0.or<- c(mean(temprisk.or[,2]),sd(temprisk.or[,2])/sqrt(repts))
  t1.or<-c(mean(temprisk.or[,3]),sd(temprisk.or[,3])/sqrt(repts))
  tau.or<- c(mean(temprisk.or[,4]),sd(temprisk.or[,4])/sqrt(repts))
  n1.or<- c(mean(temprisk.or[,5]),sd(temprisk.or[,5])/sqrt(repts))
  n2.or<-  c(mean(temprisk.or[,6]),sd(temprisk.or[,6])/sqrt(repts))
  
  # Calculate risk of the competing estimators
  temprisk.u<-matrix(0,repts,2)
  temprisk.ebt<- temprisk.u[,1]
  temprisk.ejs<-temprisk.u[,1]
  for(reps in 1:repts){
    
    a<- sureshrink.mse(U[reps,],var.u,1,0)
    temprisk.u[reps,]<- c(max(0,a$sure.est),a$t)
    muhat.ebt<- ebayesthresh(U[reps,], prior = "laplace", a = NA, bayesfac = FALSE,
                             sdev = NA, verbose = FALSE, threshrule = "median")
    
    temprisk.ebt[reps]<- sum((theta[reps,]-muhat.ebt)^2)
    temprisk.ejs[reps]<- sum((JS(U[reps,],var.u)-theta[reps,])^2)
  }
  t.u<- c(mean(temprisk.u[,2]),sd(temprisk.u[,2])/sqrt(repts))
  u<- c(mean(temprisk.u[,1])/n,sd(temprisk.u[,1])/sqrt(repts))
  ebt<- c(mean(temprisk.ebt)/n,sd(temprisk.ebt)/sqrt(repts))
  ejs<- c(mean(temprisk.ejs)/n,sd(temprisk.ejs)/sqrt(repts))
  
  
  #Generate perturbation for V.1
  for(reps in 1:repts){
    set.seed(reps)
    a<- matrix(rnorm(n*m,0,0.1),n,m)
    eta.1[reps,]<- rowMeans(a)
    rm('a')
  }
  
  #Generate V.1
  V.1<- matrix(0,repts,n)
  m1<- 0
  m2<- m1+sqrt(log(log(n)))
  set.seed(11)
  var.v<- runif(n,0.1,1)
  for(reps in 1:repts){
    mu<- mu.mat[reps,]
    n1<- length(mu[mu>temprisk.or[reps,4]])
    n2<- n-n1
    i1<-(mu>temprisk.or[reps,4])
    i2<-(mu<=temprisk.or[reps,4])
    set.seed(reps)
    V.1[reps,i1]<- rnorm(n1,m1,sqrt(var.v[i1]))
    set.seed(reps)
    V.1[reps,i2]<- rnorm(n2,m2,sqrt(var.v[i2]))
  }
  V.1<- V.1+eta.1
 
  #Side Information Risk estimation for V.1
  ttau<- seq(min(abs(V.1)),max(abs(V.1)),length.out = 20)
  temprisk.uv<-matrix(0,repts,6)
  temp<- matrix(0,repts,length(ttau))
  for(reps in 1:repts){
    VV<-abs(V.1[reps,])
    
    for(k in 1:length(ttau)){
      tau<- ttau[k]
      i1<-(VV<=tau)
      i2<- (VV>tau)
      U1<- U[reps,i1]
      sigma21<-var.u[i1]
      U2<- U[reps,i2]
      sigma22<- var.u[i2]
      a1<- sureshrink.mse(U1,sigma21,1,0)
      a2<- sureshrink.mse(U2,sigma22,1,0)
      temp[reps,k]<- max(0,a1$sure.est)+
        max(0,a2$sure.est)
    }
    index<- min(which(temp[reps,]==min(temp[reps,])))
    tau.min<- ttau[index]
    i1<-(VV<=tau.min)
    i2<- (VV>tau.min)
    U1<- U[reps,i1]
    sigma21<-var.u[i1]
    U2<- U[reps,i2]
    sigma22<- var.u[i2]
    a1<- sureshrink.mse(U1,sigma21,1,0)
    a2<- sureshrink.mse(U2,sigma22,1,0)
    temprisk.uv[reps,]<- c(max(0,a1$sure.est)+max(0,a2$sure.est),
                           a1$t,a2$t,tau.min,length(sigma21),length(sigma22))
    
  }
  uv.1<- c(mean(temprisk.uv[,1])/n,sd(temprisk.uv[,1])/sqrt(repts))
  t0.1<- c(mean(temprisk.uv[,2]),sd(temprisk.uv[,2])/sqrt(repts))
  t1.1<-c(mean(temprisk.uv[,3]),sd(temprisk.uv[,3])/sqrt(repts))
  tau.1<- c(mean(temprisk.uv[,4]),sd(temprisk.uv[,4])/sqrt(repts))
  n1.1<- c(mean(temprisk.uv[,5]),sd(temprisk.uv[,5])/sqrt(repts))
  n2.1<-  c(mean(temprisk.uv[,6]),sd(temprisk.uv[,6])/sqrt(repts))
  
  #Generate V.2
  V.2<- matrix(0,repts,n)
  m1<- 0
  m2<- m1+sqrt(log(n))
  set.seed(11)
  var.v<- runif(n,0.1,1)
  for(reps in 1:repts){
    mu<- mu.mat[reps,]
    n1<- length(mu[mu>temprisk.or[reps,4]])
    n2<- n-n1
    i1<-(mu>temprisk.or[reps,4])
    i2<-(mu<=temprisk.or[reps,4])
    set.seed(reps)
    V.2[reps,i1]<- rnorm(n1,m1,sqrt(var.v[i1]))
    set.seed(reps)
    V.2[reps,i2]<- rnorm(n2,m2,sqrt(var.v[i2]))
  }
  V.2<-V.2+eta.1
  
  #Side Information Risk estimation for V.2
  ttau<- seq(min(abs(V.2)),max(abs(V.2)),length.out = 20)
  temprisk.uv<-matrix(0,repts,6)
  temp<- matrix(0,repts,length(ttau))
  for(reps in 1:repts){
    VV<-abs(V.2[reps,])
    
    for(k in 1:length(ttau)){
      tau<- ttau[k]
      i1<-(VV<=tau)
      i2<- (VV>tau)
      U1<- U[reps,i1]
      sigma21<-var.u[i1]
      U2<- U[reps,i2]
      sigma22<- var.u[i2]
      a1<- sureshrink.mse(U1,sigma21,1,0)
      a2<- sureshrink.mse(U2,sigma22,1,0)
      temp[reps,k]<- max(0,a1$sure.est)+
        max(0,a2$sure.est)
    }
    index<- min(which(temp[reps,]==min(temp[reps,])))
    tau.min<- ttau[index]
    i1<-(VV<=tau.min)
    i2<- (VV>tau.min)
    U1<- U[reps,i1]
    sigma21<-var.u[i1]
    U2<- U[reps,i2]
    sigma22<- var.u[i2]
    a1<- sureshrink.mse(U1,sigma21,1,0)
    a2<- sureshrink.mse(U2,sigma22,1,0)
    temprisk.uv[reps,]<- c(max(0,a1$sure.est)+max(0,a2$sure.est),
                           a1$t,a2$t,tau.min,length(sigma21),length(sigma22))
    
  }
  uv.2<- c(mean(temprisk.uv[,1])/n,sd(temprisk.uv[,1])/sqrt(repts))
  t0.2<- c(mean(temprisk.uv[,2]),sd(temprisk.uv[,2])/sqrt(repts))
  t1.2<-c(mean(temprisk.uv[,3]),sd(temprisk.uv[,3])/sqrt(repts))
  tau.2<- c(mean(temprisk.uv[,4]),sd(temprisk.uv[,4])/sqrt(repts))
  n1.2<- c(mean(temprisk.uv[,5]),sd(temprisk.uv[,5])/sqrt(repts))
  n2.2<-  c(mean(temprisk.uv[,6]),sd(temprisk.uv[,6])/sqrt(repts))
  
  return(list("or"=or,"t0.or"=t0.or,"t1.or"=t1.or,"tau.or"=tau.or,
  "u"=u,"ebt"=ebt,"ejs"=ejs,
  "uv.1"=uv.1,"t0.1"=t0.1,"t1.1"=t1.1,"tau.1"=tau.1,"n1.1"=n1.1,"n2.1"=n2.1,
  "uv.2"=uv.2,"t0.2"=t0.2,"t1.2"=t1.2,"tau.2"=tau.2,"n1.2"=n1.2,"n2.2"=n2.2))

}
sfExport('repts','Ns','m')
out.sim1<- sfClusterApplyLB(1:length(Ns),wrapper)
sfStop()

for (i in 1:length(Ns)){
  
  risk.or[i]<- out.sim1[[i]]$or[1]
  risk.u[i]<- out.sim1[[i]]$u[1]
  risk.ejs[i]<- out.sim1[[i]]$ejs[1]
  risk.ebt[i]<- out.sim1[[i]]$ebt[1]
  risk.uv.1[i]<- out.sim1[[i]]$uv.1[1]
  risk.uv.2[i]<- out.sim1[[i]]$uv.2[1]
}

plotdata1<- as.data.frame(cbind(risk.u,risk.or))
names(plotdata1)<-c("SureShrink","OR")
plotdata1$Ns<- Ns

plotdata2<- as.data.frame(c(risk.uv.1,risk.uv.2,
                            risk.ebt,risk.ejs))
names(plotdata2)<-"risk"
plotdata2$legend<- as.factor(c(rep("ASUS.1",length(Ns)),
                               rep("ASUS.2",length(Ns)),rep("EBT",length(Ns)),
                               rep("EJS",length(Ns))))

plotdata2$Ns<- c(Ns,Ns,Ns,Ns)

g1<-ggplot(data=plotdata1, aes(x=Ns, y=SureShrink)) +
  geom_line(linetype="solid",size=1)+
  geom_line(data=plotdata1,aes(x=Ns,y=OR),linetype="dashed",size=1)+
  geom_line(data=plotdata2,aes(x=Ns,y=risk,color=legend))+
  geom_point(data=plotdata2,aes(x=Ns,y=risk,shape=legend),fill=NA)+
  #scale_x_continuous(breaks = round(seq(0, max(Ns), by = 20),1))+
  xlab(expression(n))+ylab("risk")+theme_bw()+
  theme(legend.position=c(0.9,0.85),legend.title=element_blank(),
        legend.background = element_rect(fill="white",
        size=1, linetype="solid", colour ="black"),
        legend.text=element_text(size=12),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))

save.image('exp2_sim1_new.RData')

#------ set 2 -------------------------------
risk.u<- matrix(0,length(Ns),1)
risk.uv.1<- risk.u
risk.uv.2<- risk.u
risk.or<- risk.u
risk.ebt<- risk.u
risk.ejs<- risk.u

sfInit(parallel=TRUE,cpu=6)
sfSource('sideinfo_lib.R')
sfLibrary(rmutil)

wrapper<- function(i){
  n<-Ns[i]
  size.1<-4*n/100
  size.2<- 16*n/100
  size<- size.1+size.2
  mu.mat<- matrix(0,repts,n)
  eta<- matrix(0,repts,n)
  eta.1<- matrix(0,repts,n)
  eta.2<- matrix(0,repts,n)
  theta<-matrix(0,repts,n)
  q<- n^{-0.5}
  
  #Generate sparse mu
  for(reps in 1:repts){
    set.seed(reps)
    mu.mat[reps,]<- c(runif(size.1,4,8),runif(size.2,1,3),rep(0,n-size))
  }
  #Generate sparse perturbation
  for(reps in 1:repts){
    set.seed(reps^2)
    eta[reps,]<- (runif(n)<=q)*rnorm(n,2,0.1)+0*(runif(n)>q)
  }
  
  #Generate theta
  theta<- mu.mat+eta
  
  #Generate U_i~N(theta_i,sigma^2)
  set.seed(reps)
  var.u<- runif(n,0.1,1)
  U<- matrix(0,repts,n)
  for(reps in 1:repts){
    for(j in 1:n){
      set.seed(reps+j)
      U[reps,j]=rnorm(1,theta[reps,j],sqrt(var.u[j]))
    }
  }
  
  # Oracle threshold by minimizing true risk
  ttau<- seq(min(mu.mat),max(mu.mat),length.out = 50)
  temprisk.or<-matrix(0,repts,6)
  temp<- matrix(0,repts,length(ttau))
  for(reps in 1:repts){
    
    mu<- mu.mat[reps,]
    VV<-abs(mu)
    muhat.or<-matrix(0,n,1)
    for(k in 1:length(ttau)){
      tau<- ttau[k]
      i1<-(VV<=tau)
      i2<- (VV>tau)
      U1<- U[reps,i1]
      sigma21<-var.u[i1]
      U2<- U[reps,i2]
      sigma22<- var.u[i2]
      a1<- sureshrink(U1,sigma21)
      a2<- sureshrink(U2,sigma22)
      muhat.or[i1]<-a1
      muhat.or[i2]<-a2
      temp[reps,k]<- sum((theta[reps,]-muhat.or)^2)
    }
    index<- max(which(temp[reps,]==min(temp[reps,])))
    tau.min<- ttau[index]
    i1<-(VV<=tau.min)
    i2<- (VV>tau.min)
    U1<- U[reps,i1]
    sigma21<-var.u[i1]
    U2<- U[reps,i2]
    sigma22<- var.u[i2]
    a1<- sureshrink.mse(U1,sigma21,1,0)
    a2<- sureshrink.mse(U2,sigma22,1,0)
    temprisk.or[reps,]<- c(temp[reps,index],
                           a1$t,a2$t,tau.min,length(sigma21),length(sigma22))
  }
  or<- c(mean(temprisk.or[,1])/n,sd(temprisk.or[,1])/sqrt(repts))
  t0.or<- c(mean(temprisk.or[,2]),sd(temprisk.or[,2])/sqrt(repts))
  t1.or<-c(mean(temprisk.or[,3]),sd(temprisk.or[,3])/sqrt(repts))
  tau.or<- c(mean(temprisk.or[,4]),sd(temprisk.or[,4])/sqrt(repts))
  n1.or<- c(mean(temprisk.or[,5]),sd(temprisk.or[,5])/sqrt(repts))
  n2.or<-  c(mean(temprisk.or[,6]),sd(temprisk.or[,6])/sqrt(repts))
  
  # Calculate risk of the competing estimators
  temprisk.u<-matrix(0,repts,2)
  temprisk.ebt<- temprisk.u[,1]
  temprisk.ejs<-temprisk.u[,1]
  for(reps in 1:repts){
    
    a<- sureshrink.mse(U[reps,],var.u,1,0)
    temprisk.u[reps,]<- c(max(0,a$sure.est),a$t)
    muhat.ebt<- ebayesthresh(U[reps,], prior = "laplace", a = NA, bayesfac = FALSE,
                             sdev = NA, verbose = FALSE, threshrule = "median")
    
    temprisk.ebt[reps]<- sum((theta[reps,]-muhat.ebt)^2)
    temprisk.ejs[reps]<- sum((JS(U[reps,],var.u)-theta[reps,])^2)
  }
  t.u<- c(mean(temprisk.u[,2]),sd(temprisk.u[,2])/sqrt(repts))
  u<- c(mean(temprisk.u[,1])/n,sd(temprisk.u[,1])/sqrt(repts))
  ebt<- c(mean(temprisk.ebt)/n,sd(temprisk.ebt)/sqrt(repts))
  ejs<- c(mean(temprisk.ejs)/n,sd(temprisk.ejs)/sqrt(repts))
  
  
  #Generate perturbation for V.1
  for(reps in 1:repts){
    set.seed(reps)
    a<- matrix(rnorm(n*m,0,0.1),n,m)
    eta.1[reps,]<- rowMeans(a)
    rm('a')
  }
  
  #Generate V.1
  V.1<- matrix(0,repts,n)
  m1<- 1
  m2<- m1+log(n)
  set.seed(11)
  var.v<- runif(n,0.1,1)
  for(reps in 1:repts){
    
    mu<- mu.mat[reps,]
    n1<- length(mu[mu>temprisk.or[reps,4]])
    n2<- n-n1
    i1<-(mu>temprisk.or[reps,4])
    i2<-(mu<=temprisk.or[reps,4])
    set.seed(reps)
    V.1[reps,i1]<- rchisq(n1,m1,ncp=0)
    set.seed(reps)
    V.1[reps,i2]<- rchisq(n2,m2,ncp=0)
  }
  V.1<- V.1+eta.1
  
  #Side Information Risk estimation for V.1
  ttau<- seq(min(abs(V.1)),max(abs(V.1)),length.out = 20)
  temprisk.uv<-matrix(0,repts,6)
  temp<- matrix(0,repts,length(ttau))
  for(reps in 1:repts){
    VV<-abs(V.1[reps,])
    
    for(k in 1:length(ttau)){
      tau<- ttau[k]
      i1<-(VV<=tau)
      i2<- (VV>tau)
      U1<- U[reps,i1]
      sigma21<-var.u[i1]
      U2<- U[reps,i2]
      sigma22<- var.u[i2]
      a1<- sureshrink.mse(U1,sigma21,1,0)
      a2<- sureshrink.mse(U2,sigma22,1,0)
      temp[reps,k]<- max(0,a1$sure.est)+
        max(0,a2$sure.est)
    }
    index<- min(which(temp[reps,]==min(temp[reps,])))
    tau.min<- ttau[index]
    i1<-(VV<=tau.min)
    i2<- (VV>tau.min)
    U1<- U[reps,i1]
    sigma21<-var.u[i1]
    U2<- U[reps,i2]
    sigma22<- var.u[i2]
    a1<- sureshrink.mse(U1,sigma21,1,0)
    a2<- sureshrink.mse(U2,sigma22,1,0)
    temprisk.uv[reps,]<- c(max(0,a1$sure.est)+max(0,a2$sure.est),
                           a1$t,a2$t,tau.min,length(sigma21),length(sigma22))
    
  }
  uv.1<- c(mean(temprisk.uv[,1])/n,sd(temprisk.uv[,1])/sqrt(repts))
  t0.1<- c(mean(temprisk.uv[,2]),sd(temprisk.uv[,2])/sqrt(repts))
  t1.1<-c(mean(temprisk.uv[,3]),sd(temprisk.uv[,3])/sqrt(repts))
  tau.1<- c(mean(temprisk.uv[,4]),sd(temprisk.uv[,4])/sqrt(repts))
  n1.1<- c(mean(temprisk.uv[,5]),sd(temprisk.uv[,5])/sqrt(repts))
  n2.1<-  c(mean(temprisk.uv[,6]),sd(temprisk.uv[,6])/sqrt(repts))
  
  #Generate V.2
  V.2<- matrix(0,repts,n)
  m1<- 1
  m2<- m1+sqrt(log(log(n)))
  set.seed(11)
  var.v<- runif(n,0.1,1)
  for(reps in 1:repts){
    
    mu<- mu.mat[reps,]
    n1<- length(mu[mu>temprisk.or[reps,4]])
    n2<- n-n1
    i1<-(mu>temprisk.or[reps,4])
    i2<-(mu<=temprisk.or[reps,4])
    set.seed(reps)
    V.2[reps,i1]<- rchisq(n1,m1,ncp=0)
    set.seed(reps)
    V.2[reps,i2]<- rchisq(n2,m2,ncp=0)
  }
  V.2<-V.2+eta.1
  
  #Side Information Risk estimation for V.2
  ttau<- seq(min(abs(V.2)),max(abs(V.2)),length.out = 20)
  temprisk.uv<-matrix(0,repts,6)
  temp<- matrix(0,repts,length(ttau))
  for(reps in 1:repts){
    VV<-abs(V.2[reps,])
    
    for(k in 1:length(ttau)){
      tau<- ttau[k]
      i1<-(VV<=tau)
      i2<- (VV>tau)
      U1<- U[reps,i1]
      sigma21<-var.u[i1]
      U2<- U[reps,i2]
      sigma22<- var.u[i2]
      a1<- sureshrink.mse(U1,sigma21,1,0)
      a2<- sureshrink.mse(U2,sigma22,1,0)
      temp[reps,k]<- max(0,a1$sure.est)+
        max(0,a2$sure.est)
    }
    index<- min(which(temp[reps,]==min(temp[reps,])))
    tau.min<- ttau[index]
    i1<-(VV<=tau.min)
    i2<- (VV>tau.min)
    U1<- U[reps,i1]
    sigma21<-var.u[i1]
    U2<- U[reps,i2]
    sigma22<- var.u[i2]
    a1<- sureshrink.mse(U1,sigma21,1,0)
    a2<- sureshrink.mse(U2,sigma22,1,0)
    temprisk.uv[reps,]<- c(max(0,a1$sure.est)+max(0,a2$sure.est),
                           a1$t,a2$t,tau.min,length(sigma21),length(sigma22))
    
  }
  uv.2<- c(mean(temprisk.uv[,1])/n,sd(temprisk.uv[,1])/sqrt(repts))
  t0.2<- c(mean(temprisk.uv[,2]),sd(temprisk.uv[,2])/sqrt(repts))
  t1.2<-c(mean(temprisk.uv[,3]),sd(temprisk.uv[,3])/sqrt(repts))
  tau.2<- c(mean(temprisk.uv[,4]),sd(temprisk.uv[,4])/sqrt(repts))
  n1.2<- c(mean(temprisk.uv[,5]),sd(temprisk.uv[,5])/sqrt(repts))
  n2.2<-  c(mean(temprisk.uv[,6]),sd(temprisk.uv[,6])/sqrt(repts))
  
  return(list("or"=or,"t0.or"=t0.or,"t1.or"=t1.or,"tau.or"=tau.or,
              "u"=u,"ebt"=ebt,"ejs"=ejs,
              "uv.1"=uv.1,"t0.1"=t0.1,"t1.1"=t1.1,"tau.1"=tau.1,"n1.1"=n1.1,"n2.1"=n2.1,
              "uv.2"=uv.2,"t0.2"=t0.2,"t1.2"=t1.2,"tau.2"=tau.2,"n1.2"=n1.2,"n2.2"=n2.2))
  
}
sfExport('repts','Ns','m')
out.sim2<- sfClusterApplyLB(1:length(Ns),wrapper)
sfStop()

for (i in 1:length(Ns)){
  
  risk.or[i]<- out.sim2[[i]]$or[1]
  risk.u[i]<- out.sim2[[i]]$u[1]
  risk.ejs[i]<- out.sim2[[i]]$ejs[1]
  risk.ebt[i]<- out.sim2[[i]]$ebt[1]
  risk.uv.1[i]<- out.sim2[[i]]$uv.1[1]
  risk.uv.2[i]<- out.sim2[[i]]$uv.2[1]
}

plotdata1<- as.data.frame(cbind(risk.u,risk.or))
names(plotdata1)<-c("SureShrink","OR")
plotdata1$Ns<- Ns

plotdata2<- as.data.frame(c(risk.uv.2,risk.uv.1,
                            risk.ebt,risk.ejs))
names(plotdata2)<-"risk"
plotdata2$legend<- as.factor(c(rep("ASUS.1",length(Ns)),
                               rep("ASUS.2",length(Ns)),rep("EBT",length(Ns)),
                               rep("EJS",length(Ns))))

plotdata2$Ns<- c(Ns,Ns,Ns,Ns)

g2<-ggplot(data=plotdata1, aes(x=Ns, y=SureShrink)) +
  geom_line(linetype="solid",size=1)+
  geom_line(data=plotdata1,aes(x=Ns,y=OR),linetype="dashed",size=1)+
  geom_line(data=plotdata2,aes(x=Ns,y=risk,color=legend))+
  geom_point(data=plotdata2,aes(x=Ns,y=risk,shape=legend),fill=NA)+
  #scale_x_continuous(breaks = round(seq(0, max(Ns), by = 20),1))+
  xlab(expression(n))+ylab("risk")+theme_bw()+
  theme(legend.position=c(0.9,0.9),legend.title=element_blank(),
        legend.background = element_rect(fill="white",
                                         size=1, linetype="solid", colour ="black"),
        legend.text=element_text(size=12),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))

save.image('exp2_sim2_new.RData')


grid.arrange(g1,g2,ncol=2)
