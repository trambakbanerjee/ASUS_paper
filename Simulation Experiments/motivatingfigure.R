#setwd('D:/Trambak/Side Information/simulations')
source('sideinfo_lib.R')
require(ggplot2)

n<-10^4
X<-matrix(0,n,1)
Y<-X
set.seed(1)
theta.x<-c(runif(2000,4,6),runif(2000,2,3),rep(0,n-4000))
theta.y<-c(runif(2000,1,2),runif(2000,1,6),rep(0,n-4000))
theta<-theta.x-theta.y
for(i in 1:n){
  set.seed(i^2)
X[i]<- rnorm(1,theta.x[i],0.5)
Y[i]<-rnorm(1,theta.y[i],0.5)
}
U<- X-Y
var.u<-rep(0.5,n)
a<-sureshrink.mse(U,var.u,1,0)
V<- abs(X+Y)
ttau<- seq(min(abs(V)),max(abs(V)),length.out = 20)
temp<-matrix(0,length(ttau),1)
for(k in 1:length(ttau)){
    tau<- ttau[k]
    i1<-(V<=tau)
    i2<- (V>tau)
    U1<- U[i1]
    sigma21<-var.u[i1]
    U2<- U[i2]
    sigma22<- var.u[i2]
    temp[k]<- max(0,sureshrink.mse(U1,sigma21,1,0)$sure.est)+
      max(0,sureshrink.mse(U2,sigma22,1,0)$sure.est)
  }
tau<- ttau[which(temp==min(temp))]
b<-min(temp)
n1<-length(U[V<=tau])
n2<-n-n1

i1<-(V<=tau)
i2<- (V>tau)
U1<- U[i1]
sigma21<-var.u[i1]
U2<- U[i2]
sigma22<- var.u[i2]
c<- sureshrink.mse(U1,sigma21,1,0)
d<-sureshrink.mse(U2,sigma22,1,0)

U<- as.data.frame(U/sqrt(var.u))
names(U)<-"V1"
U$group<- rep(0,n)
U$truth<-rep(0,n)
i1<-which(V<=tau)
i2<-which(theta==0)
U$group[i1]<-1
U$truth[i2]<-1

dat.1<-as.data.frame(U$V1[U$truth==1])
names(dat.1)<-"V1"
dat.2<-as.data.frame(U$V1[U$truth==1 & U$group==0])
names(dat.2)<-"V1"


g1<-ggplot(U,aes(x=V1)) + 
  geom_histogram(data=subset(U,truth == 1),color="white",fill = "red", alpha = 0.1,bins=75)+
  geom_histogram(data=subset(U,truth == 0),color="white",fill = "red", alpha = 0.3,bins=75)+
  #geom_point(data=dat.1,aes(x=V1,y=rep(0,length(dat.1))),shape=4,color='blue')+
  scale_x_continuous(breaks = round(seq(-10.5, 10.5, by = 1.5),1))+
  scale_y_continuous(breaks = round(seq(0, 700, by = 100),1))+
  theme_bw()+xlab(expression(Y))+ylab("counts")


g2<-ggplot(U,aes(x=V1)) + 
  geom_histogram(data=subset(U,truth == 1),color="white",fill = "red", alpha = 0.1,bins=75)+
  geom_histogram(data=subset(U,truth == 0),color="white",fill = "red", alpha = 0.3,bins=75)+
  geom_histogram(data=subset(U,group == 1),color="white",fill="blue",alpha=0.5,bins=75) +
  #geom_point(data=dat.2,aes(x=V1,y=rep(0,length(dat.2))),shape=4,color = 'green')+
  scale_x_continuous(breaks = round(seq(-10.5, 10.5, by = 1.5),1))+
  scale_y_continuous(breaks = round(seq(0, 700, by = 100),1))+
  theme_bw()+xlab(expression(Y))+ylab("counts")

g3<-ggplot(U,aes(x=V1)) + 
  geom_histogram(data=subset(U,truth == 1),color="white",fill = "red", alpha = 0.1,bins=75)+
  geom_histogram(data=subset(U,truth == 0),color="white",fill = "red", alpha = 0.3,bins=75)+
  geom_histogram(data=subset(U,group == 0),color="white",fill = "green", alpha = 0.5,bins=75) +
  #geom_point(data=dat.2,aes(x=V1,y=rep(0,length(dat.2))),shape=4,color = 'green')+
  scale_x_continuous(breaks = round(seq(-10.5, 10.5, by = 1.5),1))+
  scale_y_continuous(breaks = round(seq(0, 700, by = 100),1))+
  theme_bw()+xlab(expression(Y))+ylab("counts")

g4<-ggplot(U,aes(x=V1)) + 
  geom_histogram(data=subset(U,truth == 1),color="white",fill = "red", alpha = 0.1,bins=75)+
  geom_histogram(data=subset(U,truth == 0),color="white",fill = "red", alpha = 0.3,bins=75)+
  geom_histogram(data=subset(U,group == 0),color="white",fill = "green", alpha = 0.5,bins=75)+
  geom_histogram(data=subset(U,group == 1),color="white",fill="blue",alpha=0.5,bins=75)+
  #geom_point(data=dat.2,aes(x=V1,y=rep(0,length(dat.2))),shape=4,color = 'green')+
  scale_x_continuous(breaks = round(seq(-10.5, 10.5, by = 1.5),1))+
  scale_y_continuous(breaks = round(seq(0, 700, by = 100),1))+
  theme_bw()+xlab(expression(Y))+ylab("counts")

theta<- as.data.frame(theta)
names(theta)<-"V1"
theta$group<- rep(0,n)
i1<-which(V<=tau)
theta$group[i1]<-1

ggplot(theta,aes(x=V1)) + 
  geom_histogram(data=subset(theta,group == 0),fill = "green", alpha = 0.3,bins=75) +
  geom_histogram(data=subset(theta,group == 1),fill="blue",alpha=0.3,bins=75) +
  geom_histogram(data=theta,fill = "red", alpha = 0.1,bins=75)
