# Simulation results - Exp II - set1

require(ggplot2)
require(rmutil)
library(gridExtra)
library(snowfall)

#setwd('D:/Trambak/Side Information/simulations')
source('sideinfo_lib.R')

Ns<-seq(500,5000,100)
repts<- 1000

#------ set 1 -------------------------------
risk.u<- matrix(0,length(Ns),1)
risk.uv.1<- risk.u
risk.or<- risk.u
risk.ebt<- risk.u
risk.ejs<- risk.u

sfInit(parallel=TRUE,cpu=6)
sfSource('sideinfo_lib.R')
sfLibrary(rmutil)

wrapper<- function(i){
  n<-Ns[i]
  
  ksi.x<- matrix(0,repts,n)
  eta.x<- matrix(0,repts,n)
  theta.x<-matrix(0,repts,n)
  p.x<- n^{-0.6}
  
  #Generate sparse ksi.x
  for(reps in 1:repts){
    set.seed(reps)
    r<- runif(n)
    ksi.x[reps,]<- (r<=p.x)*runif(n,3,7)+0*(r>p.x)
  }
  #Generate perturbation eta.x
  for(reps in 1:repts){
    set.seed(reps^2)
    eta.x[reps,]<- rnorm(n,0,0.1)
  }
  #Generate theta.x
  theta.x<- ksi.x+eta.x
  
  ksi.y<- matrix(0,repts,n)
  eta.y<- matrix(0,repts,n)
  theta.y<-matrix(0,repts,n)
  p.y<- n^{-0.3}
  
  for(reps in 1:repts){
    set.seed(reps^2)
    r<- runif(n)
    ksi.y[reps,]<- (r<=p.y)*runif(n,4,4)+0*(r>p.y)
  }
  #Generate perturbation eta.y
  for(reps in 1:repts){
    set.seed(reps^3)
    eta.y[reps,]<- rnorm(n,0,0.1)
  }
  #Generate theta.y
  theta.y<- ksi.y+eta.y
  
  theta<- theta.x-theta.y
  ksi<- ksi.x-ksi.y
  
  #Generate X_i,Y_i
  var.x<- rep(1,n)
  var.y<-rep(1,n)
  X<- matrix(0,repts,n)
  Y<- X
  U<-X
  for(reps in 1:repts){
    for(j in 1:n){
      set.seed(reps+j)
      X[reps,j]=rnorm(1,theta.x[reps,j],sqrt(var.x[j]))
      Y[reps,j]=rnorm(1,theta.y[reps,j],sqrt(var.y[j]))
    }
  }
  U<-X-Y
  var.u<- sqrt(var.x+var.y)
  V<- abs(X+Y)#here k = 1
  
  # Oracle threshold by minimizing true risk
  ttau<- seq(min(ksi),max(ksi),length.out = 50)
  temprisk.or<-matrix(0,repts,6)
  temp<- matrix(0,repts,length(ttau))
  for(reps in 1:repts){
    
    mu<- ksi[reps,]
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
  
  
  #Side Information Risk estimation for V
  ttau<- seq(0,max(V),length.out = 20)
  temprisk.uv<-matrix(0,repts,6)
  temp<- matrix(0,repts,length(ttau))
  for(reps in 1:repts){
    VV<-abs(V[reps,])
    
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
  
  return(list("or"=or,"t0.or"=t0.or,"t1.or"=t1.or,"tau.or"=tau.or,"n1.or"=n1.or,"n2.or"=n2.or,
              "u"=u,"ebt"=ebt,"ejs"=ejs,
              "uv.1"=uv.1,"t0.1"=t0.1,"t1.1"=t1.1,"tau.1"=tau.1,"n1.1"=n1.1,"n2.1"=n2.1))

}
sfExport('repts','Ns')
out.sim1<- sfClusterApplyLB(1:length(Ns),wrapper)
sfStop()

for (i in 1:length(Ns)){
  
  risk.or[i]<- out.sim1[[i]]$or[1]
  risk.u[i]<- out.sim1[[i]]$u[1]
  risk.ejs[i]<- out.sim1[[i]]$ejs[1]
  risk.ebt[i]<- out.sim1[[i]]$ebt[1]
  risk.uv.1[i]<- out.sim1[[i]]$uv.1[1]
}

plotdata1<- as.data.frame(cbind(risk.u,risk.or))
names(plotdata1)<-c("SureShrink","OR")
plotdata1$Ns<- Ns

plotdata2<- as.data.frame(c(risk.uv.1,
                            risk.ebt,risk.ejs))
names(plotdata2)<-"risk"
plotdata2$legend<- as.factor(c(rep("ASUS",length(Ns)),
                               rep("EBT",length(Ns)),
                               rep("EJS",length(Ns))))

plotdata2$Ns<- c(Ns,Ns,Ns)

g1<-ggplot(data=plotdata1, aes(x=Ns, y=SureShrink)) +
  geom_line(linetype="solid",size=1)+
  geom_line(data=plotdata1,aes(x=Ns,y=OR),linetype="dashed",size=1)+
  geom_line(data=plotdata2,aes(x=Ns,y=risk,color=legend))+
  geom_point(data=plotdata2,aes(x=Ns,y=risk,shape=legend),fill=NA)+
  #scale_x_continuous(breaks = round(seq(0, max(Ns), by = 20),1))+
  xlab(expression(n))+ylab("risk")+theme_bw()+
  theme(legend.position=c(0.9,0.75),legend.title=element_blank(),
        legend.background = element_rect(fill="white",
                                         size=1, linetype="solid", colour ="black"),
        legend.text=element_text(size=12),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))

save.image('exp3_sim1.RData')

#------ set 2 -------------------------------
risk.u<- matrix(0,length(Ns),1)
risk.or<- risk.u
risk.uv.1<- risk.u
risk.ebt<- risk.u
risk.ejs<- risk.u

sfInit(parallel=TRUE,cpu=6)
sfSource('sideinfo_lib.R')
sfLibrary(rmutil)

wrapper<- function(i){
  n<-Ns[i]
  
  ksi.x<- matrix(0,repts,n)
  eta.x<- matrix(0,repts,n)
  theta.x<-matrix(0,repts,n)
  p.x<- n^{-0.6}
  
  #Generate sparse ksi.x
  for(reps in 1:repts){
    set.seed(reps)
    r<- runif(n)
    ksi.x[reps,]<- (r<=p.x)*runif(n,3,7)+0*(r>p.x)
  }
  #Generate perturbation eta.x
  for(reps in 1:repts){
    set.seed(reps^2)
    eta.x[reps,]<- rnorm(n,0,0.1)
  }
  #Generate theta.x
  theta.x<- ksi.x+eta.x
  
  ksi.y<- matrix(0,repts,n)
  eta.y<- matrix(0,repts,n)
  theta.y<-matrix(0,repts,n)
  p.y<- n^{-0.3}
  
  for(reps in 1:repts){
    set.seed(reps^2)
    r<- runif(n)
    ksi.y[reps,]<- (r<=p.y)*runif(n,4,4)+0*(r>p.y)
  }
  #Generate perturbation eta.y
  for(reps in 1:repts){
    set.seed(reps^3)
    eta.y[reps,]<- rnorm(n,0,0.1)
  }
  #Generate theta.y
  theta.y<- ksi.y+eta.y
  
  theta<- theta.x-theta.y
  ksi<- ksi.x-ksi.y
  
  #Generate X_i,Y_i
  set.seed(n)
  var.x<- runif(n,0.1,1)
  set.seed(n^2)
  var.y<-runif(n,0.1,1)
  X<- matrix(0,repts,n)
  Y<- X
  U<-X
  for(reps in 1:repts){
    for(j in 1:n){
      set.seed(reps+j)
      X[reps,j]=rnorm(1,theta.x[reps,j],sqrt(var.x[j]))
      Y[reps,j]=rnorm(1,theta.y[reps,j],sqrt(var.y[j]))
    }
  }
  U<-X-Y
  var.u<- sqrt(var.x+var.y)
  k = sqrt(var.x/var.y)
  V<- abs(X+k*Y)
  
  # Oracle threshold by minimizing true risk
  ttau<- seq(min(ksi),max(ksi),length.out = 50)
  temprisk.or<-matrix(0,repts,6)
  temp<- matrix(0,repts,length(ttau))
  for(reps in 1:repts){
    
    mu<- ksi[reps,]
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
  
  
  #Side Information Risk estimation for V
  ttau<- seq(0,max(V),length.out = 20)
  temprisk.uv<-matrix(0,repts,6)
  temp<- matrix(0,repts,length(ttau))
  for(reps in 1:repts){
    VV<-abs(V[reps,])
    
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
  
  return(list("or"=or,"t0.or"=t0.or,"t1.or"=t1.or,"tau.or"=tau.or,"n1.or"=n1.or,"n2.or"=n2.or,
              "u"=u,"ebt"=ebt,"ejs"=ejs,
              "uv.1"=uv.1,"t0.1"=t0.1,"t1.1"=t1.1,"tau.1"=tau.1,"n1.1"=n1.1,"n2.1"=n2.1))
  
}
sfExport('repts','Ns')
out.sim2<- sfClusterApplyLB(1:length(Ns),wrapper)
sfStop()

for (i in 1:length(Ns)){
  
  risk.or[i]<- out.sim2[[i]]$or[1]
  risk.u[i]<- out.sim2[[i]]$u[1]
  risk.ejs[i]<- out.sim2[[i]]$ejs[1]
  risk.ebt[i]<- out.sim2[[i]]$ebt[1]
  risk.uv.1[i]<- out.sim2[[i]]$uv.1[1]
}

plotdata1<- as.data.frame(cbind(risk.u,risk.or))
names(plotdata1)<-c("SureShrink","OR")
plotdata1$Ns<- Ns

plotdata2<- as.data.frame(c(risk.uv.1,
                            risk.ebt,risk.ejs))
names(plotdata2)<-"risk"
plotdata2$legend<- as.factor(c(rep("ASUS",length(Ns)),
                               rep("EBT",length(Ns)),
                               rep("EJS",length(Ns))))

plotdata2$Ns<- c(Ns,Ns,Ns)

g2<-ggplot(data=plotdata1, aes(x=Ns, y=SureShrink)) +
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

save.image('exp3_sim2.RData')


grid.arrange(g1,g2,ncol=2)
