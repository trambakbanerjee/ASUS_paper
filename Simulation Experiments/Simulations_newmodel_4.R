# Simulation experiment 4 - Choice of K

#setwd('D:/Trambak/Side Information/simulations')
source('sideinfo_lib.R')
require(ggplot2)
library(gridExtra)
library(snowfall)

#-------------------------------
# Set 1 -- How many groups?

# Generate data
repts<- 10
n<-10^4
X<-matrix(0,n,repts)
Y<-X
theta<-X

for(reps in 1:repts){
  set.seed(reps)
  theta.x<-c(runif(2000,4,6),runif(2000,2,3),rep(0,n-4000))
  theta.y<-c(runif(2000,1,2),runif(2000,1,6),rep(0,n-4000))
  theta[,reps]<-theta.x-theta.y
  for(i in 1:n){
    set.seed(i^2)
    X[i,reps]<- rnorm(1,theta.x[i],0.5)
    Y[i,reps]<-rnorm(1,theta.y[i],0.5)
  }
}
U<- X-Y
var.u<-rep(0.5,n)
V<- abs(X+Y)

group.count<- 1:10
mse<- matrix(0,length(group.count),1)
n.min<- matrix(0,length(group.count),1)
n.min[1]<- n
for (gc in 1:length(group.count)){
  
  out<- ASUS_K(U,var.u,V,group.count[gc])
  mse[gc]<- out$sure.est
  if(gc>1){
    n.min[gc]<-min(out$nn) 
  }
  print(gc)
}

plotdata1<- as.data.frame(mse)
names(plotdata1)<-"risk"
plotdata1$gc<- group.count

g1<-ggplot(data=plotdata1, aes(x=gc, y=risk)) +
  geom_line(linetype="solid",size=1)+
  geom_point(data=plotdata1,aes(x=gc,y=risk),fill=NA)+
  scale_x_continuous(breaks = round(seq(1, max(group.count), by = 1),1))+
  xlab('K')+ylab("risk of ASUS")+theme_bw()+
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))

save.image('exp4_sim1.RData')

#---------------------------------------
# Set 2 --- groups with log n separation

Ns<-seq(500,10000,100)
repts<- 1000

risk.u<- matrix(0,length(Ns),1)
risk.block.u<- risk.u
risk.uv<- risk.u

sfInit(parallel=TRUE,cpu=6)
sfSource('sideinfo_lib.R')
sfLibrary(rmutil)

wrapper<- function(i){
  n<-Ns[i]
  # Generate data
  repts<- 10
  X<-matrix(0,n,repts)
  Y<-X
  theta<-X
  n1<- 0.2*n
  n2<- n-2*n1
  for(reps in 1:repts){
    set.seed(reps)
    theta.x<-c(runif(n1,4,6),runif(n1,2,3),rep(0,n2))
    theta.y<-c(runif(n1,1,2),runif(n1,1,6),rep(0,n2))
    theta[,reps]<-theta.x-theta.y
    for(i in 1:n){
      set.seed(i^2)
      X[i,reps]<- rnorm(1,theta.x[i],0.5)
      Y[i,reps]<-rnorm(1,theta.y[i],0.5)
    }
  }
  U<- X-Y
  var.u<-rep(0.5,n)
  V<- abs(X+Y)
  
  # Calculate risk of the competing estimators
  temprisk.u<-matrix(0,repts,2)
  for(reps in 1:repts){
    
    a<- sureshrink.mse(U[,reps],var.u,1,0)
    temprisk.u[reps,]<- c(max(0,a$sure.est),a$t)
  }
  t.u<- c(mean(temprisk.u[,2]),sd(temprisk.u[,2])/sqrt(repts))
  u<- c(mean(temprisk.u[,1])/n,sd(temprisk.u[,1])/sqrt(repts))
  
  # Calculate risk of block estimator
  temprisk.u<-matrix(0,repts,1)
  for(reps in 1:repts){
    U.sliced <- slice(U[,reps],log(n))
    var.u.sliced<- slice(var.u,log(n))
    a<- matrix(0,length(U.sliced),1)
    for (k in 1:length(U.sliced)){
      
      a[k]<- max(0,sureshrink.mse(U.sliced[[k]],var.u.sliced[[k]],1,0)$sure.est)
    }
    temprisk.u[reps]<- sum(a)
  }
  block.u<- mean(temprisk.u)/n
  
  #Side Information Risk estimation
  temprisk.uv<-matrix(0,repts,6)
  ttau<- seq(min(V),max(V),length.out = 20)
  temp<- matrix(0,repts,length(ttau))
  for(reps in 1:repts){
    VV<-abs(V[,reps])
    
    for(k in 1:length(ttau)){
      tau<- ttau[k]
      i1<-(VV<=tau)
      i2<- (VV>tau)
      U1<- U[i1,reps]
      sigma21<-var.u[i1]
      U2<- U[i2,reps]
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
    U1<- U[i1,reps]
    sigma21<-var.u[i1]
    U2<- U[i2,reps]
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
  
  return(list("u"=u,"block.u" = block.u,
              "uv.1"=uv.1,"t0.1"=t0.1,"t1.1"=t1.1,"tau.1"=tau.1,"n1.1"=n1.1,"n2.1"=n2.1))
  
}
sfExport('repts','Ns')
out.sim2<- sfClusterApplyLB(1:length(Ns),wrapper)
sfStop()

for (i in 1:length(Ns)){
  
  risk.u[i]<- out.sim2[[i]]$u[1]
  risk.block.u[i]<- out.sim2[[i]]$block.u[1]
  risk.uv[i]<- out.sim2[[i]]$uv.1[1]
}

plotdata1<- as.data.frame(risk.u)
names(plotdata1)<-"risk"
plotdata1$legend<- as.factor(rep("SureShrink",length(Ns)))
plotdata1$Ns<- Ns

plotdata2<- as.data.frame(c(risk.uv,risk.block.u))
names(plotdata2)<-"risk"
plotdata2$legend<- as.factor(c(rep("ASUS K = 2",length(Ns)),
                               rep("K = n/log(n)",length(Ns))))

plotdata2$Ns<- c(Ns,Ns)

g2<-ggplot(data=plotdata1, aes(x=Ns, y=risk,color=legend)) +
  geom_line(linetype="solid",size=1)+
  geom_point(data=plotdata1,aes(x=Ns,y=risk,shape=legend),fill=NA)+
  geom_line(data=plotdata2,aes(x=Ns,y=risk,color=legend))+
  geom_point(data=plotdata2,aes(x=Ns,y=risk,shape=legend),fill=NA)+
  scale_x_continuous(breaks = round(seq(0, max(Ns), by = 1000),1))+
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2))+
  xlab(expression(n))+ylab("risk")+theme_bw()+
  theme(legend.position=c(0.8,0.9),legend.title=element_blank(),
        legend.background = element_rect(fill="white",
                                         size=1, linetype="solid", colour ="black"),
        legend.text=element_text(size=12),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))

save.image('exp4_sim2.RData')

grid.arrange(g1,g2,ncol=2)
