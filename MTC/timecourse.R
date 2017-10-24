require(ggplot2)
library(gridExtra)

#setwd('D:/Trambak/Side Information/Real Data/timecourse')
source('sideinfo_lib.R')

# read data
dat<- read.csv('ImmuTimeCourse.txt', header = TRUE, sep = " ")
dat$ProbeSetId<- as.character(dat$ProbeSetId)
dat.0<- dat[,c(1,3,4,5,6)]
dat.2<- dat[,c(1,11,12,13,14)]
dat.4<- dat[,c(1,18,19,20,21)]
dat.6<- dat[,c(1,25,26,27,28)]
dat.9<- dat[,c(1,33,34,35,36)]
dat.24<-dat[,c(1,41,42,43,44)]

Y.1<-as.matrix(rowMeans(asinh(dat.4[,2:5]-dat.0[,2:5]))) 
v.1<-as.matrix(apply(asinh(dat.4[,2:5]-dat.0[,2:5]),1,var))
v.1[which(v.1==0)]<-0.1

Y.2<-as.matrix(rowMeans(asinh(dat.24[,2:5]-dat.0[,2:5]))) 
v.2<-as.matrix(apply(asinh(dat.24[,2:5]-dat.0[,2:5]),1,var))
v.2[which(v.2==0)]<-0.1

Y<- Y.2-Y.1
v.y<-v.2+v.1 
kappa<- sqrt(v.2/v.1)
X<- abs(Y.2+kappa*Y.1)
n<-length(X)

# Calculate risk of SureShrink estimator
ss<- sureshrink.mse(Y,v.y,1,0)
risk.u<- 100*(max(0,ss$sure.est))/n
t.u<- ss$t

#Side Information Risk estimation
ttau<- seq(0,12,length.out = 500)
ntau<- length(ttau)
risk.uv<-matrix(0,length(ttau),3)
t.uv<- matrix(0,length(ttau),2)
for(k in 1:length(ttau)){
  tau<- ttau[k]
  i1<-(X<=tau)
  i2<- (X>tau)
  U1<- Y[i1]
  sigma21<-v.y[i1]
  U2<- Y[i2]
  sigma22<- v.y[i2]
  temp.1<-sureshrink.mse(U1,sigma21,1,0)
  temp.2<- sureshrink.mse(U2,sigma22,1,0)
  t.uv[k,]<- c(temp.1$t,temp.2$t)
  risk.uv[k,]<- c(100*(max(0,temp.1$sure.est)+max(0,temp.2$sure.est))/n,temp.1$sure.est,
                  temp.2$sure.est)
}
index<- min(which(risk.uv[,1]==min(risk.uv[,1])))
risk.uv.min<- risk.uv[index,1]
tau.min<- ttau[index]
t.uv.min<- t.uv[index,]

plotdata1<- as.data.frame(c((risk.u[1]*matrix(1,ntau,1)),risk.uv[,1]))
names(plotdata1)<-"risk"
plotdata1$legend<- as.factor(c(rep("SS",ntau),rep("SI",ntau)))
plotdata1$N<- c(ttau,ttau)
taudata1<-as.data.frame(cbind(tau*matrix(1,5,1),seq(0,min(risk.uv[,1]),length.out = 5)))

g1<- ggplot(data=plotdata1, aes(x=N, y=risk)) +
  geom_line(aes(linetype=legend),size=1)+
  geom_point(data=taudata1,aes(x=tau.min,y=risk.uv.min),shape=16,
             size=4,fill="red",color="red")+
  xlab(expression(tau))+ylab("risk(%)")+theme_bw()+
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15))

plot(g1)

i1<-(X<=tau.min)
i2<- (X>tau.min)
U1<- Y[i1]
sigma21<-v.y[i1]
U2<- Y[i2]
sigma22<- v.y[i2]
c<- sureshrink.mse(U1,sigma21,1,0)
d<-sureshrink.mse(U2,sigma22,1,0)

est.1<- sureshrink(U1,sigma21)
est.2<- sureshrink(U2,sigma22)
est<-matrix(0,n,1)
est[i1]<- est.1
est[i2]<-est.2

U<- as.data.frame(Y)
names(U)<-"V1"
U$group<- rep(0,n)
i1<-which(X<=tau.min)
U$group[i1]<-1

plotdata2<-as.data.frame(U$V1[U$group==1])
names(plotdata2)<-"V1"
plotdata3<-as.data.frame(U$V1[U$group==0])
names(plotdata3)<-"V1"

g2<-ggplot(U,aes(x=V1)) + 
  geom_histogram(data=subset(U,group == 0),fill = "green", alpha = 0.5,bins=50) +
  geom_histogram(data=subset(U,group == 1),fill="blue",alpha=0.5,bins=50) +
  geom_histogram(data=U,color="white",fill = "red", alpha = 0.1,bins=50)+
  #coord_cartesian(ylim=c(0, 800))+
  #scale_x_continuous(breaks = round(seq(min(U$V1), max(U$V1), by = 0.5),1))+
  #scale_y_continuous(breaks = round(seq(0, 700, by = 100),1))+
  theme_bw()+xlab("Y")+ylab("counts")+
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14))

plot(g2)

g3<-ggplot(U,aes(x=V1)) + 
  geom_histogram(data=subset(U,group == 1),fill="blue",alpha=0.5,bins=50) +
  geom_histogram(data=subset(U,group == 0),fill = "green", alpha = 0.5,bins=50) +
  geom_histogram(data=U,color="white",fill = "red", alpha = 0.1,bins=50)+
  coord_cartesian(ylim=c(0, 100))+
  #scale_x_continuous(breaks = round(seq(min(U$V1), max(U$V1), by = 0.5),1))+
  scale_y_continuous(breaks = round(seq(0, 100, by = 25),1))+
  theme_bw()+xlab("Y")+ylab("counts")+
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14))

plot(g3)

grid.arrange(g1,g2,ncol=2)

#------- Heatmap ----------------------
# par(mfrow=c(1,1))
ii<- order(abs(Y),decreasing = FALSE)
YY<-abs(Y)
U<-YY[ii]
V<-X[ii]
heat.data<-cbind(U,V)

colors.1 = c(seq(0,0.5,length=3),seq(0.5,1,length=3),seq(1,1.5,length=3),
             seq(1.5,2,length=3),seq(2,4,length=5))
mycolor<- colorRampPalette(c("black","yellow","orange","red"))(16);

image(1,1:n,t(heat.data[,1]),
      xlab='U',ylab='',col=mycolor,xaxt='n',yaxt='n',breaks=colors.1);

colors.2 = c(seq(0,1,length=3),seq(1,2,length=3),seq(2,3,length=3),
             seq(3,4,length=3),seq(4,60,length=5))
mycolor<- colorRampPalette(c("black","yellow","orange","red"))(16);
image(1,1:n,t(heat.data[,2]),
      xlab='V',ylab='',col=mycolor,xaxt='n',yaxt='n',breaks=colors.2);


