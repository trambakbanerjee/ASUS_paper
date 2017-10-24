require(ggplot2)
library(gridExtra)

setwd('D:/Trambak/Side Information/Real Data/timecourse')
source('D:/Trambak/Side Information/simulations/sideinfo_lib.R')

# read data
dat<- read.csv('ImmuTimeCourse.txt', header = TRUE, sep = " ")
dat$ProbeSetId<- as.character(dat$ProbeSetId)
dat.0<- dat[,c(1,3,4,5,6)]
dat.2<- dat[,c(1,11,12,13,14)]
dat.4<- dat[,c(1,18,19,20,21)]
dat.6<- dat[,c(1,25,26,27,28)]
dat.9<- dat[,c(1,33,34,35,36)]
dat.24<-dat[,c(1,41,42,43,44)]

id<-c(1,2,3,1,2,4,1,3,4)


v.1<-as.matrix(apply(asinh(dat.4[,2:5]-dat.0[,2:5]),1,var))
v.1[which(v.1==0)]<-0.1

v.2<-as.matrix(apply(asinh(dat.24[,2:5]-dat.0[,2:5]),1,var))
v.2[which(v.2==0)]<-0.1

v.y<-v.2+v.1 
kappa<- sqrt(v.2/v.1)

ttau<- seq(0,12,length.out = 500)
ntau<- length(ttau)

risk.uv<- matrix(0,ntau,3)
risk.u<-matrix(0,3,1)
t.u<- matrix(0,3,1)
t.uv<-matrix(0,3,2)
tau.uv<-matrix(0,3,1)

d4<- dat.4[,2:5]
d0<- dat.0[,2:5]
d24<- dat.24[,2:5]

for (k in 1:3){
  s<-3*k-2
  e<-3*k
  c<- id[s:e]
  Y.1<-as.matrix(rowMeans(asinh(d4[c]-d0[c]))) 
  Y.2<-as.matrix(rowMeans(asinh(d24[c]-d0[c]))) 
  Y<- Y.2-Y.1
  
  X<- abs(Y.2+kappa*Y.1)
  n<-length(X)
  # Calculate risk of SureShrink estimator
  ss<- sureshrink.mse(Y,v.y,1,0)
  risk.u[k]<- 100*(max(0,ss$sure.est))/n
  t.u[k]<- ss$t
  
  #Side Information Risk estimation
  temp<- matrix(0,length(ttau),2)
  for(kk in 1:length(ttau)){
    tau<- ttau[kk]
    i1<-(X<=tau)
    i2<- (X>tau)
    U1<- Y[i1]
    sigma21<-v.y[i1]
    U2<- Y[i2]
    sigma22<- v.y[i2]
    temp.1<-sureshrink.mse(U1,sigma21,1,0)
    temp.2<- sureshrink.mse(U2,sigma22,1,0)
    temp[kk,]<- c(temp.1$t,temp.2$t)
    risk.uv[kk,k]<- 100*(max(0,temp.1$sure.est)+max(0,temp.2$sure.est))/n
  }
  index<- min(which(risk.uv[,k]==min(risk.uv[,k])))
  tau.uv[k]<- ttau[index]
  t.uv[k,]<- temp[index,]
}

risk.u.test<-matrix(0,3,1)
risk.uv.test<- risk.u.test
p<-matrix(0,3,2)
gain<- risk.u.test

for(k in 1:3){
  s<-3*k-2
  e<-3*k
  c<- id[s:e]
  Y.1<-as.matrix(rowMeans(asinh(d4[-c]-d0[-c]))) 
  Y.2<-as.matrix(rowMeans(asinh(d24[-c]-d0[-c]))) 
  Y<- Y.2-Y.1
  
  X<- abs(Y.2+kappa*Y.1)
  n<-length(X)
  out<- sureshrink.mse(Y,v.y,2,t.u[k])
  risk.u.test[k]<-100*out$sure.est/n
  U1<- Y[X<=tau.uv[k]]
  var.U1<- v.y[X<=tau.uv[k]]
  U2<- Y[X>tau.uv[k]]
  var.U2<- v.y[X>tau.uv[k]]
  p[k,]<- c(length(U1),length(U2))
  out1<- sureshrink.mse(U1,var.U1,2,t.uv[k,1])
  out2<- sureshrink.mse(U2,var.U2,2,t.uv[k,2])
  risk.uv.test[k]<- 100*(out1$sure.est+out2$sure.est)/n
  
  gain[k]<- 100*(risk.u.test[k]-risk.uv.test[k])/risk.u.test[k]
}

