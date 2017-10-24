
require(utils)
require(ggplot2)
require(plotly)
library(gridExtra)
library(grid)
#setwd('D:/Trambak/Side Information/Real Data/rnaseq')
source('sideinfo_lib.R')

options(stringsAsFactors = FALSE)
dat<- read.csv('all samples-rna-seq.csv', header = TRUE, sep = ",")
dat$id<- as.character(seq(1,(nrow(dat))))
N<- nrow(dat)

X<- asinh(cbind(dat$VZV.1,dat$VZV.2,dat$VZV.3))
VZV<- rowMeans(X)
v.y<-read.csv('variances.csv', header = TRUE, sep = ",")

V<- asinh(cbind(dat[,2:4],dat[,8:15]))
V<- rowMeans(V)
tt<- seq(-10^{-6},4,length.out = 100)
ntau<- length(tt)

risk.uv<- matrix(0,ntau,3)

out<- sureshrink.mse(VZV,v.y,1,0)
risk.u<- 100*out$sure.est/N
t.u<-out$t
t.uv<- matrix(0,ntau,2)

for (i in 1:ntau){
  
  tau<- tt[i]
  U1<- VZV[V<=tau]
  var.U1<- v.y[V<=tau]
  U2<- VZV[V>tau]
  var.U2<- v.y[V>tau]
  out1<- sureshrink.mse(U1,var.U1,1,0)
  out2<- sureshrink.mse(U2,var.U2,1,0)
  t.uv[i,]<-c(out1$t,out2$t)
  risk.uv[i,]<- c(100*(max(0,out1$sure.est)+max(0,out2$sure.est))/N,
                   out1$sure.est,out2$sure.est)
  print(i)
}
idx<- min(which(risk.uv[,1]==min(risk.uv[,1])))
tau.min<- tt[idx]
risk.uv.min<- risk.uv[idx,1]
t.uv.min<- t.uv[idx,]

plotdata1<- as.data.frame(c((risk.u[1]*matrix(1,ntau,1)),risk.uv[,1]))
names(plotdata1)<-"risk"
plotdata1$legend<- as.factor(c(rep("SS",ntau),rep("SI",ntau)))
plotdata1$N<- c(tt,tt)
taudata1<-as.data.frame(cbind(tau*matrix(1,5,1),seq(0,min(risk.uv[,1]),length.out = 5)))

g1<- ggplot(data=plotdata1, aes(x=N, y=risk)) +
  geom_line(aes(linetype=legend),size=1)+
  geom_point(data=taudata1,aes(x=tau.min,y=risk.uv.min),shape=16,
             size=4,fill="red",color="red")+
  xlab(expression(tau))+ylab("risk(%)")+theme_bw()+
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))

i1<-(V<=tau.min)
i2<- (V>tau.min)
U1<- VZV[i1]
sigma21<-v.y[i1]
U2<- VZV[i2]
sigma22<- v.y[i2]
c<- sureshrink.mse(U1,sigma21,1,0)
d<-sureshrink.mse(U2,sigma22,1,0)
est.1<- sureshrink(U1,sigma21)
est.2<- sureshrink(U2,sigma22)
est<-matrix(0,N,1)
est[i1]<- est.1
est[i2]<-est.2

top50<- order(est,decreasing = TRUE)[1:50]
top50<-as.data.frame(dat$gene_id[top50])

bottom50<- order(est,decreasing = FALSE)[1:50]
bottom50<-as.data.frame(dat$gene_id[bottom50])


U<- as.data.frame(VZV)
names(U)<-"V1"
U$group<- rep(0,N)
i1<-which(V<=tau.min)
U$group[i1]<-1

plotdata2<-as.data.frame(U$V1[U$group==1])
names(plotdata2)<-"V1"
plotdata3<-as.data.frame(U$V1[U$group==0])
names(plotdata3)<-"V1"

g2<-ggplot(U,aes(x=V1)) + 
  geom_histogram(data=subset(U,group == 0),fill = "green", alpha = 0.5,bins=75) +
  geom_histogram(data=subset(U,group == 1),fill="blue",alpha=0.5,bins=75) +
  geom_histogram(data=U,color="white",fill = "red", alpha = 0.1,bins=75)+
  coord_cartesian(xlim=c(0, 8),ylim=c(0, 5500))+
  scale_x_continuous(breaks = round(seq(min(U$V1), 8, by = 1),1))+
  scale_y_continuous(breaks = round(seq(0, 25000, by = 1000),1))+
  theme_bw()+xlab("VZV")+ylab("counts")+
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14))

grid.arrange(g1,g2,ncol=2)

#------- Heatmap ----------------------
# par(mfrow=c(1,1))

U<-VZV

V.1<- rowMeans(dat[,2:4])
V.1<- as.matrix(asinh(V.1))
V.2<- rowMeans(dat[,8:10])
V.2<- as.matrix(asinh(V.2))
V.3<- rowMeans(dat[,11:12])
V.3<- as.matrix(asinh(V.3))
V.4<- rowMeans(dat[,13:15])
V.4<- as.matrix(asinh(V.4))
V<-cbind(V.1,V.2,V.3,V.4)

heat.data<-cbind(U,V)

colors = c(seq(0,1,length=1),seq(1,3,length=3),seq(3,6,length=4),seq(6,8,length=3),
           seq(8,14,length=7))
mycolor<- colorRampPalette(c("black","yellow","orange","red"))(17);

image(1:dim(heat.data)[2],1:N,t(heat.data),
      xlab='U.train',ylab='',col=mycolor,xaxt='n',yaxt='n',breaks=colors);
axis(1, at=c(1,2,3,4,5),labels=c("VZV (3)","HELF (3)",
                                   "HT1080 (3)","IFNG (2)","IFNA (3)"))
abline(v=c(1.5,2.5,3.5,4.5,5.5),col='white',lwd=10)

test<-as.matrix(seq(0,14,length=18))
image(1:dim(test)[2],1:18,t(test),
      xlab='',col=mycolor,xaxt='n',breaks=colors)

#---------- comparison with SureShrink ---------------
V<- asinh(cbind(dat[,2:4],dat[,8:15]))
V<- rowMeans(V)
theta.u<- sureshrink(VZV,v.y)/sqrt(v.y)
theta.uv<- matrix(0,N,1)

tau<- tau.min
i1<- (V<=tau)
i2<-(V>tau)
U1<- VZV[i1]
var.U1<- v.y[i1]
U2<- VZV[i2]
var.U2<- v.y[i2]
theta1.uv<- sureshrink(U1,var.U1)/sqrt(var.U1)
theta2.uv<- sureshrink(U2,var.U2)/sqrt(var.U2)
theta.uv[i1]<- theta1.uv
theta.uv[i2]<- theta2.uv

p <- pt(theta.u, 2, lower.tail = TRUE, log.p = FALSE)
pp<- cbind(p,1-p)
p.u<- 2*apply(pp,1,min)
padj.u<-as.matrix(cbind(p.u,p.adjust(p.u,method="fdr")))

p <- pt(theta.uv, 2, lower.tail = TRUE, log.p = FALSE)
pp<- cbind(p,1-p)
p.uv<- 2*apply(pp,1,min)
padj.uv<-as.matrix(cbind(p.uv,p.adjust(p.uv,method="fdr")))

pvals<- cbind(padj.u,padj.uv)   
#write.csv(pvals,file="C:/Users/trambakb/Documents/Study/Research/Side Information/Real Data/rnaseq/pavls.csv")


