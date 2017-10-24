# Simulation results - Exp I - set1

require(ggplot2)
require(rmutil)
library(gridExtra)
library(snowfall)

#setwd('D:/Trambak/Side Information/simulations')
source('sideinfo_lib.R')

M<-seq(10,200,10)
n<-5000
repts<- 500

#------ set 1 -------------------------------
risk.u<- matrix(0,length(M),1)
risk.uv.1<- risk.u
risk.uv.2<- risk.u
risk.uv.3<- risk.u
risk.uv.4<- risk.u
risk.or<- risk.u
risk.ebt<- risk.u
risk.ejs<- risk.u

sfInit(parallel=TRUE,cpu=6)
sfSource('sideinfo_lib.R')
sfLibrary(rmutil)

wrapper<- function(i){
  m<-M[i]
  mu.mat<- matrix(0,repts,n)
  eta<- matrix(0,repts,n)
  eta.1<- matrix(0,repts,n)
  eta.2<- matrix(0,repts,n)
  theta<-matrix(0,repts,n)
  q<- n^{-0.5}
  
  #Generate sparse mu
  for(reps in 1:repts){
    set.seed(reps)
    mu.mat[reps,]<- c(runif(50,6,7),runif(200,2,3),rep(0,n-250))
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
  
  # oracle threshold calculation and true loss
  ttau<- seq(min(mu.mat),max(mu.mat),length.out = 50)
  temprisk.or<-matrix(0,repts,6)
  temp<- matrix(0,repts,length(ttau))
  for(reps in 1:repts){
    
    mu<- mu.mat[reps,]
    VV<-abs(mu)
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
      temp[reps,k]<- max(0,a1$sure.est)+max(0,a2$sure.est)
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
    temprisk.or[reps,]<- c(max(0,a1$sure.est)+max(0,a2$sure.est),
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
    a<- matrix(rlaplace(n*m,0,4),n,m)
    eta.1[reps,]<- rowMeans(a)
    rm('a')
  }
  
  #Generate V.1
  V.1<- matrix(0,repts,n)
  V.1<-mu.mat+eta.1
 
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
  
  #Generate perturbation for V.2
  for(reps in 1:repts){
    set.seed(reps)
    a<- matrix(rchisq(n*m,10,0),n,m)
    eta.2[reps,]<- rowMeans(a)
    rm('a')
  }
  
  #Generate V.2
  V.2<- matrix(0,repts,n)
  V.2<-mu.mat+eta.2
  
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
  
 #Generate V.3
  V.3<- matrix(0,repts,n)
  sig<-mu.mat
  nn<- length(which(mu.mat==0))
  set.seed(reps)
  sig[mu.mat==0]<-rlnorm(nn,0,5/sqrt(m))
  V.3<-sig+eta.1
  
  #Side Information Risk estimation for V.3
  ttau<- seq(min(abs(V.3)),max(abs(V.3)),length.out = 20)
  temprisk.uv<-matrix(0,repts,6)
  temp<- matrix(0,repts,length(ttau))
  for(reps in 1:repts){
    VV<-abs(V.3[reps,])
    
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
  uv.3<- c(mean(temprisk.uv[,1])/n,sd(temprisk.uv[,1])/sqrt(repts))
  t0.3<- c(mean(temprisk.uv[,2]),sd(temprisk.uv[,2])/sqrt(repts))
  t1.3<-c(mean(temprisk.uv[,3]),sd(temprisk.uv[,3])/sqrt(repts))
  tau.3<- c(mean(temprisk.uv[,4]),sd(temprisk.uv[,4])/sqrt(repts))
  n1.3<- c(mean(temprisk.uv[,5]),sd(temprisk.uv[,5])/sqrt(repts))
  n2.3<-  c(mean(temprisk.uv[,6]),sd(temprisk.uv[,6])/sqrt(repts))
  
  #Generate V.4
  V.4<- matrix(0,repts,n)
  sig<-mu.mat
  nn<- length(which(mu.mat==0))
  set.seed(reps)
  sig[mu.mat==0]<-rt(nn,round(2*m/10),0)
  V.4<-sig+eta.2
  
  #Side Information Risk estimation for V.4
  ttau<- seq(min(abs(V.4)),max(abs(V.4)),length.out = 20)
  temprisk.uv<-matrix(0,repts,6)
  temp<- matrix(0,repts,length(ttau))
  for(reps in 1:repts){
    VV<-abs(V.4[reps,])
    
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
  uv.4<- c(mean(temprisk.uv[,1])/n,sd(temprisk.uv[,1])/sqrt(repts))
  t0.4<- c(mean(temprisk.uv[,2]),sd(temprisk.uv[,2])/sqrt(repts))
  t1.4<-c(mean(temprisk.uv[,3]),sd(temprisk.uv[,3])/sqrt(repts))
  tau.4<- c(mean(temprisk.uv[,4]),sd(temprisk.uv[,4])/sqrt(repts))
  n1.4<- c(mean(temprisk.uv[,5]),sd(temprisk.uv[,5])/sqrt(repts))
  n2.4<-  c(mean(temprisk.uv[,6]),sd(temprisk.uv[,6])/sqrt(repts))
  
  
  return(list("or"=or,"t0.or"=t0.or,"t1.or"=t1.or,"tau.or"=tau.or,"n1.or"=n1.or,"n2.or"=n2.or,
              "u"=u,"ebt"=ebt,"ejs"=ejs,
              "uv.1"=uv.1,"t0.1"=t0.1,"t1.1"=t1.1,"tau.1"=tau.1,"n1.1"=n1.1,"n2.1"=n2.1,
              "uv.2"=uv.2,"t0.2"=t0.2,"t1.2"=t1.2,"tau.2"=tau.2,"n1.2"=n1.2,"n2.2"=n2.2,
              "uv.3"=uv.3,"t0.3"=t0.3,"t1.3"=t1.3,"tau.3"=tau.3,"n1.3"=n1.3,"n2.3"=n2.3,
              "uv.4"=uv.4,"t0.4"=t0.4,"t1.4"=t1.4,"tau.4"=tau.4,"n1.4"=n1.4,"n2.4"=n2.4))

}
sfExport('repts','M','n')
out.sim1<- sfClusterApplyLB(1:length(M),wrapper)
sfStop()

for (i in 1:length(M)){
  
  risk.or[i]<- out.sim1[[i]]$or[1]
  risk.u[i]<- out.sim1[[i]]$u[1]
  risk.ejs[i]<- out.sim1[[i]]$ejs[1]
  risk.ebt[i]<- out.sim1[[i]]$ebt[1]
  risk.uv.1[i]<- out.sim1[[i]]$uv.1[1]
  risk.uv.2[i]<- out.sim1[[i]]$uv.2[1]
  risk.uv.3[i]<- out.sim1[[i]]$uv.3[1]
  risk.uv.4[i]<- out.sim1[[i]]$uv.4[1]
}

plotdata1<- as.data.frame(cbind(risk.u,risk.or))
names(plotdata1)<-c("SureShrink","OR")
plotdata1$M<- M

plotdata2<- as.data.frame(c(risk.uv.1,risk.uv.2,risk.uv.3,risk.uv.4,
                            risk.ebt,risk.ejs))
names(plotdata2)<-"risk"
plotdata2$legend<- as.factor(c(rep("ASUS.1",length(M)),
                               rep("ASUS.2",length(M)),rep("ASUS.3",length(M)),
                               rep("ASUS.4",length(M)),rep("EBT",length(M)),
                               rep("EJS",length(M))))

plotdata2$M<- c(M,M,M,M,M,M)

g1<-ggplot(data=plotdata1, aes(x=M, y=SureShrink)) +
  geom_line(linetype="solid",size=1)+
  geom_line(data=plotdata1,aes(x=M,y=OR),linetype="dashed",size=1)+
  geom_line(data=plotdata2,aes(x=M,y=risk,color=legend))+
  geom_point(data=plotdata2,aes(x=M,y=risk,shape=legend),fill=NA)+
  scale_x_continuous(breaks = round(seq(0, max(M), by = 20),1))+
  xlab(expression(m))+ylab("risk")+theme_bw()+
  theme(legend.position=c(0.9,0.75),legend.title=element_blank(),
        legend.background = element_rect(fill="white",
        size=1, linetype="solid", colour ="black"),
        legend.text=element_text(size=12),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))

save.image('exp1_sim1.RData')

#------ set 2 -------------------------------
risk.u<- matrix(0,length(M),1)
risk.uv.1<- risk.u
risk.uv.2<- risk.u
risk.uv.3<- risk.u
risk.uv.4<- risk.u
risk.or<- risk.u
risk.ebt<- risk.u
risk.ejs<- risk.u

sfInit(parallel=TRUE,cpu=6)
sfSource('sideinfo_lib.R')
sfLibrary(rmutil)

wrapper<- function(i){
  m<-M[i]
  mu.mat<- matrix(0,repts,n)
  eta<- matrix(0,repts,n)
  eta.1<- matrix(0,repts,n)
  eta.2<- matrix(0,repts,n)
  theta<-matrix(0,repts,n)
  q<- n^{-0.5}
  
  #Generate sparse mu
  for(reps in 1:repts){
    set.seed(reps)
    mu.mat[reps,]<- c(runif(200,4,8),runif(800,1,3),rep(0,n-1000))
  }
  #Generate sparse perturbation
  for(reps in 1:repts){
    set.seed(reps^2)
    eta[reps,]<- (runif(n)<=q)*rnorm(n,2,0.1)+0*(runif(n)>q)
  }
  
  #Generate theta
  theta<- mu.mat+eta
  
  #Generate U_i~N(theta_i,sigma^2)
  set.seed(1)
  var.u<- rep(1,n)
  U<- matrix(0,repts,n)
  for(reps in 1:repts){
    for(j in 1:n){
      set.seed(reps+j)
      U[reps,j]=rnorm(1,theta[reps,j],sqrt(var.u[j]))
    }
  }
  
  # oracle threshold calculation and true loss
  ttau<- seq(min(mu.mat),max(mu.mat),length.out = 50)
  temprisk.or<-matrix(0,repts,6)
  temp<- matrix(0,repts,length(ttau))
  for(reps in 1:repts){
    
    mu<- mu.mat[reps,]
    VV<-abs(mu)
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
      temp[reps,k]<- max(0,a1$sure.est)+max(0,a2$sure.est)
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
    temprisk.or[reps,]<- c(max(0,a1$sure.est)+max(0,a2$sure.est),
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
    a<- matrix(rlaplace(n*m,0,4),n,m)
    eta.1[reps,]<- rowMeans(a)
    rm('a')
  }
  
  #Generate V.1
  V.1<- matrix(0,repts,n)
  V.1<-mu.mat+eta.1
  
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
  
  #Generate perturbation for V.2
  for(reps in 1:repts){
    set.seed(reps)
    a<- matrix(rchisq(n*m,5,0),n,m)
    eta.2[reps,]<- rowMeans(a)
    rm('a')
  }
  
  #Generate V.2
  V.2<- matrix(0,repts,n)
  V.2<-mu.mat+eta.2
  
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
  
  #Generate V.3
  V.3<- matrix(0,repts,n)
  sig<-mu.mat
  nn<- length(which(mu.mat==0))
  set.seed(reps)
  sig[mu.mat==0]<-rlnorm(nn,0,5/sqrt(m))
  V.3<-sig+eta.1
  
  #Side Information Risk estimation for V.3
  ttau<- seq(min(abs(V.3)),max(abs(V.3)),length.out = 20)
  temprisk.uv<-matrix(0,repts,6)
  temp<- matrix(0,repts,length(ttau))
  for(reps in 1:repts){
    VV<-abs(V.3[reps,])
    
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
  uv.3<- c(mean(temprisk.uv[,1])/n,sd(temprisk.uv[,1])/sqrt(repts))
  t0.3<- c(mean(temprisk.uv[,2]),sd(temprisk.uv[,2])/sqrt(repts))
  t1.3<-c(mean(temprisk.uv[,3]),sd(temprisk.uv[,3])/sqrt(repts))
  tau.3<- c(mean(temprisk.uv[,4]),sd(temprisk.uv[,4])/sqrt(repts))
  n1.3<- c(mean(temprisk.uv[,5]),sd(temprisk.uv[,5])/sqrt(repts))
  n2.3<-  c(mean(temprisk.uv[,6]),sd(temprisk.uv[,6])/sqrt(repts))
  
  #Generate V.4
  V.4<- matrix(0,repts,n)
  sig<-mu.mat
  nn<- length(which(mu.mat==0))
  set.seed(reps)
  sig[mu.mat==0]<-rt(nn,round(2*m/10),0)
  V.4<-sig+eta.2
  
  #Side Information Risk estimation for V.4
  ttau<- seq(min(abs(V.4)),max(abs(V.4)),length.out = 20)
  temprisk.uv<-matrix(0,repts,6)
  temp<- matrix(0,repts,length(ttau))
  for(reps in 1:repts){
    VV<-abs(V.4[reps,])
    
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
  uv.4<- c(mean(temprisk.uv[,1])/n,sd(temprisk.uv[,1])/sqrt(repts))
  t0.4<- c(mean(temprisk.uv[,2]),sd(temprisk.uv[,2])/sqrt(repts))
  t1.4<-c(mean(temprisk.uv[,3]),sd(temprisk.uv[,3])/sqrt(repts))
  tau.4<- c(mean(temprisk.uv[,4]),sd(temprisk.uv[,4])/sqrt(repts))
  n1.4<- c(mean(temprisk.uv[,5]),sd(temprisk.uv[,5])/sqrt(repts))
  n2.4<-  c(mean(temprisk.uv[,6]),sd(temprisk.uv[,6])/sqrt(repts))
  
  return(list("or"=or,"t0.or"=t0.or,"t1.or"=t1.or,"tau.or"=tau.or,"n1.or"=n1.or,"n2.or"=n2.or,
              "u"=u,"ebt"=ebt,"ejs"=ejs,
              "uv.1"=uv.1,"t0.1"=t0.1,"t1.1"=t1.1,"tau.1"=tau.1,"n1.1"=n1.1,"n2.1"=n2.1,
              "uv.2"=uv.2,"t0.2"=t0.2,"t1.2"=t1.2,"tau.2"=tau.2,"n1.2"=n1.2,"n2.2"=n2.2,
              "uv.3"=uv.3,"t0.3"=t0.3,"t1.3"=t1.3,"tau.3"=tau.3,"n1.3"=n1.3,"n2.3"=n2.3,
              "uv.4"=uv.4,"t0.4"=t0.4,"t1.4"=t1.4,"tau.4"=tau.4,"n1.4"=n1.4,"n2.4"=n2.4))
  
}
sfExport('repts','M','n')
out.sim2<- sfClusterApplyLB(1:length(M),wrapper)
sfStop()

for (i in 1:length(M)){
  
  risk.or[i]<- out.sim2[[i]]$or[1]
  risk.u[i]<- out.sim2[[i]]$u[1]
  risk.ejs[i]<- out.sim2[[i]]$ejs[1]
  risk.ebt[i]<- out.sim2[[i]]$ebt[1]
  risk.uv.1[i]<- out.sim2[[i]]$uv.1[1]
  risk.uv.2[i]<- out.sim2[[i]]$uv.2[1]
  risk.uv.3[i]<- out.sim2[[i]]$uv.3[1]
  risk.uv.4[i]<- out.sim2[[i]]$uv.4[1]
}

plotdata1<- as.data.frame(cbind(risk.u,risk.or))
names(plotdata1)<-c("SureShrink","OR")
plotdata1$M<- M

plotdata2<- as.data.frame(c(risk.uv.1,risk.uv.2,risk.uv.3,risk.uv.4,
                            risk.ebt,risk.ejs))
names(plotdata2)<-"risk"
plotdata2$legend<- as.factor(c(rep("ASUS.1",length(M)),
                               rep("ASUS.2",length(M)),rep("ASUS.3",length(M)),
                               rep("ASUS.4",length(M)),rep("EBT",length(M)),
                               rep("EJS",length(M))))

plotdata2$M<- c(M,M,M,M,M,M)

g2<-ggplot(data=plotdata1, aes(x=M, y=SureShrink)) +
  geom_line(linetype="solid",size=1)+
  geom_line(data=plotdata1,aes(x=M,y=OR),linetype="dashed",size=1)+
  geom_line(data=plotdata2,aes(x=M,y=risk,color=legend))+
  geom_point(data=plotdata2,aes(x=M,y=risk,shape=legend),fill=NA)+
  scale_x_continuous(breaks = round(seq(0, max(M), by = 20),1))+
  xlab(expression(m))+ylab("risk")+theme_bw()+
  theme(legend.position=c(0.9,0.65),legend.title=element_blank(),
        legend.background = element_rect(fill="white",
        size=1, linetype="solid", colour ="black"),
        legend.text=element_text(size=12),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))

save.image('exp1_sim2.RData')


grid.arrange(g1,g2,ncol=2)
