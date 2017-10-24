require(MASS)
require(rwt)
require(wavethresh)
require(EbayesThresh)
require(snowfall)


# SURE Shrink estimator of mean from Donoho (1994)
sureshrink<- function(d,var){
  p<- length(d)
  if (p>0){
    dd<- d/sqrt(var)
    td<- sqrt(2*log(p))
    gam<- ((log2(p))^(1.5))/sqrt(p)
    ts<- sure(dd)
    sd2<- (1/p)*(sum(dd^2)-p)
    tt<- ts
    if(sd2<=gam){
      tt<-td
    }
    d.ss<-sqrt(var)*softTh(dd,tt)
    return(d.ss)
  }
  if (p==0){
    return(0)
  }
}

# Stein's Unbiased Risk Estimate using SURE threshold
sureshrink.mse<- function(d,var,type,t){
  
  pp<- length(d)
  if (pp==0){
    return(list("sure.est"=0,"t"=0))
  }
  if(type==1){
    if (pp>0){
      dd<-d/sqrt(var)
      td<- sqrt(2*log(pp))
      #gamma<- (var*((log2(pp))^(1.5)))/sqrt(pp)
      gam<- ((log2(pp))^(1.5))/sqrt(pp)
      ts<- sure(dd)
      sd2<- (1/pp)*(sum(dd^2)-pp)
      tt<- ts
      if (sd2<=gam){
        tt<-td
      }
      sure.est<- sum(var)-2*sum(var*(abs(dd)<=tt))+
        sum(var*(pmin(tt,abs(dd)))^2)
      return(list("sure.est"=sure.est,"t"=tt))
    }
  }
  if(type==2){
    if (pp>0){
      dd<-d/sqrt(var)
      tt<-t
      sure.est<- sum(var)-2*sum(var*(abs(dd)<=tt))+
        sum(var*(pmin(tt,abs(dd)))^2)
      return(list("sure.est"=sure.est,"t"=tt))
    } 
  }
}

# Covariate Assisted Sure Shrink
sureshrink.cass<- function(d,covar,vard,tau,mu){
  
  ntau<- length(tau)
  if (ntau>1){
    # do cass estimation by finding the minimizing tau
    mse.cass<- matrix(0,ntau,1)
    temp<- matrix(0,ntau,2)
    mse.or<- matrix(0,ntau,1)
    for (i in 1:ntau){
      
      tt<- tau[i]
      z1<- d[covar<=tt]
      z2<- d[covar>tt]
      p1<- length(z1)
      p2<- length(z2)
      temp[i,]<-c(sureshrink.mse(z1,p1,vard),sureshrink.mse(z2,p2,vard))
      mse.cass[i]<- sureshrink.mse(z1,p1,vard)+
        sureshrink.mse(z2,p2,vard)
      muhat.cass<- matrix(NA,length(d),1)
      muhat.cass[covar<=tt]<- sureshrink(z1,p1,vard)
      muhat.cass[covar>tt]<- sureshrink(z2,p2,vard)
      mse.or[i]<- sum((muhat.cass-mu)^2)
    }
    taumin<- tau[min(which(mse.cass==min(mse.cass)))]
    z1<- d[covar<=taumin]
    z2<- d[covar>taumin]
    p1<- length(z1)
    p2<- length(z2)
    mse.min<- sureshrink.mse(z1,p1,vard)+sureshrink.mse(z2,p2,vard)
    muhat.cass<- matrix(NA,length(d),1)
    muhat.cass[covar<=taumin]<- sureshrink(z1,p1,2)
    muhat.cass[covar>taumin]<- sureshrink(z2,p2,2)
    mseor.min<- sum((muhat.cass-mu)^2)
    return(list("est"=muhat.cass,"mse.min"=mse.min,"mseor.min"=mseor.min,
                "mse"=mse.cass,"mseor"=mse.or,"taumin"=taumin,
                "temp"=temp))
  }
  
}

# True loss
truerisk<- function(U,var,mu,theta,params){
  
  sfInit(parallel=TRUE,cpu=8)
  sfLibrary(rwt)
  wrapper<-function(idx){
    repts<- dim(U)[1]
    p<- dim(U)[2]
    tt<- params[idx,]
    p1<-length(mu[mu<=tt[1]])
    p2<- p-p1
    i1<-(mu<=tt[1])
    i2<- (mu>tt[1])
    U1<- as.matrix(U[,i1])
    var1<-var[i1]
    U2<- as.matrix(U[,i2])
    var2<- var[i2]
    loss1<-matrix(0,repts,1)
    loss2<-loss1
    if(p1>0){
      theta1<-theta[i1]
      t1<-tt[2]
      for(reps in 1:repts){
        UU<- U1[reps,]/sqrt(var1)
        
        g1<-sqrt(var1)*softTh(UU,t1)
        
        loss1[reps]<- sum((theta1-g1)^2)
      }
    }
    if(p1==0){
      loss1<- rep(0,repts)
    }
    if(p2>0){
      theta2<-theta[i2]
      t2<-tt[3]
      for(reps in 1:repts){
        UU<- U2[reps,]/sqrt(var2)
        
        g2<-sqrt(var2)*softTh(UU,t2)
        
        loss2[reps]<- sum((theta2-g2)^2)
      }
    }
    if(p2==0){
      loss2<- rep(0,repts)
    }
    return(mean(loss1+loss2))
    
  }
  sfExport('U','var','mu','theta','params')
  out<-sfClusterApplyLB(1:2000,wrapper)
  out<-matrix(unlist(out),ncol=1,byrow=TRUE)
  sfStop()
  index<- min(which(out==min(out)))
  tt<- params[index,]
  risk<- min(out)
  
  return(list("risk"=risk,"a"=tt[1],"t1"=tt[2],"t2"=tt[3]))
}

# Group-Linear (Weinstein et al 2015:
# https://github.com/MaZhuang/grouplinear/blob/master/functions.R)

grouplinear <- function(x,v,nbreak=floor(length(x)^(1/3)) ){
  # default: bin log(v) into same NUMBER (=n^(1/3) of intervals
  n <- length(x)
  splitby=cut(log(v),breaks=nbreak, labels=F)
  xsub <- split(x,splitby)
  vsub <- split(v,splitby)
  indexsub <- split(1:n,splitby)
  thetahatsub <- mapply(spher,xsub,vsub)
  indexsub.unlist <- as.vector( unlist(indexsub) )
  thetahatsub.unlist <- as.vector( unlist(thetahatsub) )
  thetahat <- thetahatsub.unlist[order(indexsub.unlist)]	
  return(thetahat)
}

## spherically symmetric estimator with c_n = c^*_n
spher <- function(x.,v.){
  n. <- length(x.)
  if ( (n.==1) | (var(x.)==0) ) x. else {
    cstar <- max( 1-2*( max(v.)/mean(v.) )/(n.-1), 0)
    bhat <- min( cstar*mean(v.)/var(x.), 1 )
    x. - bhat*(x. - mean(x.))
  }
}

# extended James-Stein (equation (7.3) in Xie et al, 2012)
JS <- function(X,A){
  p <- length(X)
  muhat <- sum( X*(1/A) )/sum(1/A)
  b <- 1 - (p-3) / sum( (X-muhat)^2/A )
  return( muhat + max(0,b) * (X-muhat) )
}


# ASUS calculation with K groups
ASUS_K <- function(U,var.U,V,K){
  
  repts<- dim(U)[2]
  n<- dim(U)[1]
  sure.est<- matrix(0,1,repts)
  
  if (K==1){
    # simple SureShrink (no side info)
    for (reps in 1:repts){
      sure.est[reps]<-(sureshrink.mse(U[,reps],var.u,1,0)$sure.est)
    }
    return(list("sure.est"=mean(sure.est)/n,"n1"=n,"n2"=0))
  }
  if (K==2){
    n1 = matrix(0,1,repts)
    n2 = n1
    # 2 group ASUS
    for (reps in 1:repts){
      ttau<- seq(min(abs(V)),max(abs(V)),length.out = 20)
      temp<-matrix(0,length(ttau),1)
      for(k in 1:length(ttau)){
        tau<- ttau[k]
        i1<-(V[,reps]<=tau)
        i2<- (V[,reps]>tau)
        U1<- U[i1,reps]
        sigma21<-var.U[i1]
        U2<- U[i2,reps]
        sigma22<- var.U[i2]
        temp[k]<- max(0,sureshrink.mse(U1,sigma21,1,0)$sure.est)+
          max(0,sureshrink.mse(U2,sigma22,1,0)$sure.est)
      }
      tau<- ttau[which(temp==min(temp))]
      sure.est[reps]<-min(temp)
      n1[reps]<-length(U[V[,reps]<=tau,reps])
      n2[reps]<-n-n1[reps]
    }
    return(list("sure.est"=mean(sure.est)/n,"nn"=c(mean(n1),mean(n2))))
    
  }
  if (K>2){
    # a bit of snowfall will help!
    
    vec<- 1:20
    len<- K-1
    a<-t(combn(length(vec), len, function(x) vec[x]))
    ttau.temp<- seq(1,(max(abs(V))-1),length.out = 20)
    ttau<- matrix(0,nrow(a),ncol(a))
    for (k in 1:ncol(a)){
      ttau[,k]<-ttau.temp[a[,k]] 
    }
    rm("ttau.temp","a","vec")
    
    sfInit(parallel=TRUE,cpu=6)
    sfSource('sideinfo_lib.R')
    sfLibrary(rmutil)
    wrapper<- function(reps){
      
      nn = matrix(0,K,1)
      temp<-matrix(0,nrow(ttau),1)
      UU<- U[,reps]
      VV<- V[,reps]
      
      for(i in 1:nrow(ttau)){
        tau<- ttau[i,]
        temp[i]<-check_index(UU,var.U,VV,tau)$mse
      }
      tau<- ttau[which(temp==min(temp)),]
      sure.est<-min(temp)
      nn[1]<-length(UU[VV<=tau[1]])
      nn[K]<-length(UU[VV>tau[(K-1)]])
      for (k in 2:(K-1)){
        
        nn[k]<- length(UU[(VV>tau[k-1] & VV<=tau[k])])
      } 
      return(list("sure.est"=sure.est,"nn"=nn))
      
    }
    # for (reps in 1:repts){
    #   
    #   temp<-matrix(0,nrow(ttau),1)
    #   UU<- U[,reps]
    #   VV<- V[,reps]
    #   
    #   for(i in 1:nrow(ttau)){
    #     tau<- ttau[i,]
    #     temp[i]<-check_index(UU,var.U,VV,tau)$mse
    #   }
    #   tau<- ttau[which(temp==min(temp)),]
    #   sure.est[reps]<-min(temp)
    #   nn[1,reps]<-length(UU[VV<=tau[1]])
    #   nn[K,reps]<-length(UU[VV>tau[(K-1)]])
    #   for (k in 2:(K-1)){
    #     
    #     nn[k,reps]<- length(UU[(VV>tau[k-1] & VV<=tau[k])])
    #   }
    # }
    # return(list("sure.est"=mean(sure.est)/n,"nn"=rowMeans(nn)))
  sfExport('ttau','n','U','var.U','V','K')
  out.sim1<- sfClusterApplyLB(1:repts,wrapper)
  sfStop()
  sure.est<- matrix(unlist(out.sim1),ncol=(K+1),byrow=TRUE)
  return(list("sure.est"=mean(sure.est[,1])/n,"nn"=colMeans(sure.est[,-1])))
  }
}

# check index if K>2
check_index<- function(UU,var.UU,VV,cutpoints){
  
  K = length(cutpoints)+1
  temp<- matrix(0,K,1)
  i.1 = matrix(FALSE,length(VV),K)
  i.1[,1] <- (VV<=cutpoints[1])
  temp[1]<- max(0,sureshrink.mse(UU[i.1[,1]],var.UU[i.1[,1]],1,0)$sure.est)
  i.1[,K]<- (VV>cutpoints[K-1])
  temp[K]<- max(0,sureshrink.mse(UU[i.1[,K]],var.UU[i.1[,K]],1,0)$sure.est)
  for (k in 2:(K-1)){
    
    i.1[,k]<- (VV>cutpoints[k-1] & VV<=cutpoints[k])
    temp[k]<- max(0,sureshrink.mse(UU[i.1[,k]],var.UU[i.1[,k]],1,0)$sure.est)
    
  }
  
  return(list("index" =i.1,"mse"=sum(temp)))
}

# slice
slice<-function(x,n) {
  N<-length(x);
  lapply(seq(1,N,n),function(i) x[i:min(i+n-1,N)])
}

