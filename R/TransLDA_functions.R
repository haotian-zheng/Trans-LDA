library(lassoshooting)
library(pROC)
library(numDeriv) 

# Q-aggregation
agg.fun<- function(B, X.test,y.test, lasso=F, total.step=30, lda=F, lda.exp=F, tuning = 1){
  if(sum(B==0)==ncol(B)*nrow(B)){
    return(rep(0,nrow(B)))
  }
  if(nrow(X.test==1)){
    Sig.test = t(X.test) %*% X.test
  }else{
    mu1.test<- apply(X.test[y.test<0,],2,mean)
    mu2.test<- apply(X.test[y.test>0,],2,mean)
    S1 = (sum(y.test<0)-1)*stats::cov(t(t(X.test[y.test<0,])-mu1.test))
    S2 = (sum(y.test>0)-1)*stats::cov(t(t(X.test[y.test>0,])-mu2.test))
    Sig.test = (S1 + S2)/(length(y.test)-2)
  }
  colnames(B)<-NULL
  if(lasso){
    zero.col<-apply(B,2, function(x) sum(abs(x))==0)
    B<-B[,which(!zero.col)]
    B<-apply(B, 2, function(x) x/as.numeric((t(x)%*%Sig.test%*%x)^(1/2)))
    cor.mat<-cor(X.test%*%B)
    cat(min(abs(cor.mat)),max(abs(cor.mat[row(cor.mat)!=col(cor.mat)])),'\n')
    lam.theta=cv.glmnet(X.test%*%B, y.test, lambda=seq(0.2,2,0.2)*sqrt(2*log(ncol(B))/length(y.test)))$lambda.min
    theta.hat<-as.numeric(glmnet(X.test%*%B, y.test, lambda=lam.theta)$beta)
    beta<-as.numeric(B%*%theta.hat)
  }else if(lda){
    if(lda.exp){
      theta.pow = apply(B,2, function(x) 1/2*t(x)%*%Sig.test%*%x*length(y.test)-t(x)%*%t(X.test)%*%y.test)
      theta.pow = theta.pow-min(theta.pow)
      theta.hat<- exp(
        
        -theta.pow
        
        /2)
      if(sum(theta.hat)!=0) theta.hat=theta.hat/sum(theta.hat)
      beta<-as.numeric(B%*%theta.hat)
      theta.old=theta.hat
      for(ss in 1:total.step){
        #exponential aggregating
        theta.pow = apply(B,2, function(x) 1/2*t(x)%*%Sig.test%*%x*length(y.test)-t(x)%*%t(X.test)%*%y.test)/4 -
          apply(X.test%*%B,2, function(x) sum((x-X.test%*%beta)^2))/8
        theta.pow = theta.pow-min(theta.pow)
        theta.hat<- exp(
          
          -theta.pow)
        if(sum(theta.hat)!=0) theta.hat=theta.hat/sum(theta.hat)
        beta=3*beta/4+ as.numeric(B%*%theta.hat)/4
        if(sum(abs(theta.hat-theta.old))<10^(-6)){break}
        theta.old=theta.hat
      }
      if(ss==total.step){
        theta.hat<- exp(
          
          -apply(B,2, function(x) 1/2*t(x)%*%Sig.test%*%x*length(y.test)-t(x)%*%t(X.test)%*%y.test)
          
          /2)
        theta.hat=theta.hat/sum(theta.hat)
        beta<-as.numeric(B%*%theta.hat)
      }
      cat('ss.lda=',ss,'\n')
    }else{
      XtX=t(B)%*%Sig.test%*%B
      Xty=t(B)%*%t(X.test)%*%y.test/length(y.test)
      lambda=sqrt(2*2*log(ncol(B))/length(y.test))
      
      if(tuning == 1){
        simpleFun = function(x) {
          return(1/2*t(x)%*%XtX%*%x-t(x)%*%Xty+lambda*sum(abs(x)))
        }
        x0 = rep(0, ncol(B))
        gd = coor.descent(simpleFun,x0,step.size=0.05, max.iter = 600)
        theta.hat = gd$x
      }
      if(tuning == 0){
        simpleFun = function(x) {
          return(1/2*t(x)%*%XtX%*%x-t(x)%*%Xty)
        }
        x0 = rep(0, ncol(B))
        gd = grad.descent(simpleFun,x0,step.size=0.05, max.iter = 200)
        theta.hat = gd$x
        print("yes")
      }
      beta<-as.numeric(B%*%theta.hat)
      
    }
  }else{
    theta.hat<- exp(-colSums((y.test-X.test%*%B)^2)/2)
    if(sum(theta.hat)!=0) theta.hat=theta.hat/sum(theta.hat)
    beta<-as.numeric(B%*%theta.hat)
    theta.old=theta.hat
    for(ss in 1:total.step){
      #exponential aggregating
      theta.hat<- exp(-colSums((y.test-X.test%*%B)^2)/4 + apply(X.test%*%B,2, function(x) sum((x-X.test%*%beta)^2))/8)
      if(sum(theta.hat)!=0) theta.hat=theta.hat/sum(theta.hat)
      beta=3*beta/4+ as.numeric(B%*%theta.hat)/4
      if(sum(abs(theta.hat-theta.old))<10^(-6)){break}
      theta.old=theta.hat
    }
    if(ss==total.step){
      theta.hat<- exp(-colSums((y.test-X.test%*%B)^2)/2)
      theta.hat=theta.hat/sum(theta.hat)
      beta<-as.numeric(B%*%theta.hat)
    }
    cat('ss=',ss,'\n')
  }
  
  list(theta=theta.hat, beta=beta)
}

# Trans-LDA with adaptive steps for unknown informative studies
Trans.lasso <- function(X, y, n.vec, I.til, l1=T, lasso=F, lda = F, lda.exp = F, tuning = 1, lambdas = NA){
  M= length(n.vec)-1
  #step 1
  if(is.na(I.til)){
    X0.til<-X[ind.set(n.vec,1),]
    y0.til<-y[ind.set(n.vec,1)]
  }else{
    X0.til<-X[I.til,] #used for aggregation
    y0.til<-y[I.til]
    X<- X[-I.til,]
    y<-y[-I.til]
  }
  #step 2
  Rhat <- rep(0, M+1)
  p<- ncol(X)
  n.vec.I = n.vec
  n.vec.I[1]<- n.vec[1]-length(I.til)
  ind.1<-ind.set(n.vec.I,1)
  y.hat<-y[ind.1]
  sub.1<- sample(ind.1, floor(length(ind.1)/2), replace=F)
  sub.1.c<-setdiff(ind.1,sub.1)
  for(k in 2: (M+1)){
    ind.k<-ind.set(n.vec.I,k)
    Xty.k <- t(X[ind.k,])%*%y[ind.k]/n.vec.I[k] - t(X[ind.1,])%*%y[ind.1]/ n.vec.I[1]
    margin.T<-sort(abs(Xty.k),decreasing=T)[1:(43)]
    Rhat[k] <-  sum(margin.T^2)
  }
  Tset<- list()
  k0=0
  kk.list<-unique(rank(Rhat[-1]))
  for(kk in 1:length(kk.list)){#use pi.hat as selection rule
    Tset[[k0+kk]]<- which(rank(Rhat[-1]) <= kk.list[kk])
  }
  k0=length(Tset)
  
  Tset<- unique(Tset)
  cat("Tset=",length(Tset),'\n')
  
  beta.T<-list()
  result.T = las.kA.tune(X=X, y=y, A0=NULL, n.vec=n.vec.I, l1=l1, lambda=lambdas) #211119 change: atrans-lda each candidate same tuning parameter
  beta.T[[1]] <- result.T$beta
  for(kk in 1:length(Tset)){#use pi.hat as selection rule
    T.k <- Tset[[kk]]
    beta.T[[kk+1]] <-las.kA.tune(X=X, y=y, A0=T.k, n.vec=n.vec.I, l1=l1, lambda=result.T$lambda)$beta
  }
  beta.T<-beta.T[!duplicated((beta.T))]
  beta.T<- as.matrix(as.data.frame(beta.T))
  
  agg.re <- agg.fun(B=beta.T, X.test=X0.til, y.test=y0.til, lasso=lasso, lda=lda, lda.exp=lda.exp, tuning = tuning)
  
  return(list(beta.hat=agg.re$beta, theta.hat=agg.re$theta, rank.pi=rank(Rhat[-1]), beta.T=beta.T))
}

# Naive Trans-LDA with BIC for parameter tuning
las.kA.tune<-function(X, y, A0, n.vec, l1=T, lambda=NA){
  p<-ncol(X)
  ind.kA<- ind.set(n.vec, c(1, A0+1))
  size.A0<- length(A0)
  ind.1<-1:n.vec[1]
  #compute Sig.X.h (sigma hat)
  Sig.X.h = matrix(0,p,p)
  xty = 0
  if(is.na(lambda)){
    lambda.min = bic.las.kA(max.lambda = 3, X.test = X, y.test = y, A0=A0, n.vec=n.vec, l1=T)
    cat(lambda.min, '\n')
  }
  for (k in 1:(size.A0+1)) {
    ind.k<-ind.set(n.vec, k)
    y.k<- y[ind.k]
    x.k<- X[ind.k,]
    mu1.h<- apply(x.k[y.k<0,],2,mean)
    mu2.h<- apply(x.k[y.k>0,],2,mean)
    S1 = (sum(y.k<0)-1)*stats::cov(t(t(x.k[y.k<0,])-mu1.h))
    S2 = (sum(y.k>0)-1)*stats::cov(t(t(x.k[y.k>0,])-mu2.h))
    Sig.X.k.h = (S1 + S2)/(n.vec[k]-2)
    
    if(k==1){
      if(!is.na(lambda)){
        lambda.min = lambda
      }
      Sig.X.0.h = Sig.X.k.h
      mu1.0.h = mu1.h
      mu2.0.h = mu2.h
      result.ini = lassoshooting(XtX=Sig.X.0.h, Xty=(mu2.0.h-mu1.0.h), maxit=1000, lambda=lambda.min*sqrt(2*2*log(p)/length(ind.1)))
      cat(result.ini$iterations, '\n')
      if(result.ini$iterations>=1000) return("error")
      beta.ini<-as.numeric(result.ini$coef)
    }
    result.k = lassoshooting(XtX=Sig.X.k.h, Xty=(mu2.h-mu1.h-Sig.X.k.h%*%beta.ini),lambda=lambda.min*sqrt(2*2*log(p)/min(length(ind.1),length(ind.k))))
    cat(result.k$iterations, '\n')
    # if(result.k$iterations>=1000) return("error")
    delta.k<-as.numeric(result.k$coef)
    
    alpha.k = n.vec[k]/length(ind.kA)
    Sig.X.h = Sig.X.h + alpha.k*Sig.X.k.h
    xty = xty + alpha.k*(mu2.h-mu1.h-Sig.X.k.h%*%delta.k)
  }
  
  if(size.A0 >= 0){
    ind.1<-1:n.vec[1]
    y.A<- y[ind.1]
    if(l1){
      y.A<-y[ind.kA]
    }else{
    }
    result.kA = lassoshooting(XtX=Sig.X.h, Xty=xty,lambda=lambda.min*sqrt(2*2*log(p)/length(ind.kA)))
    cat(result.kA$iterations, '\n')
    # if(result.kA$iterations>=1000) return("error")
    beta.kA=as.numeric(result.kA$coef)
  }
  return(list(beta = as.numeric(beta.kA), lambda = lambda.min))
}

# BIC parameter tuning 
bic.las.kA <- function(max.lambda = 1.9, X.test, y.test, A0, n.vec, l1=T) {
  seq.lambda = seq(max.lambda,0,-0.02)
  err = rep(Inf, length(seq.lambda))
  for (i in 1:length(seq.lambda)) {
    ind.1<-1:n.vec[1]
    y.0<- y.test[ind.1]
    x.0<- X.test[ind.1,]
    mu1.0.h<- apply(x.0[y.0<0,],2,mean)
    mu2.0.h<- apply(x.0[y.0>0,],2,mean)
    S1 = (sum(y.0<0)-1)*stats::cov(t(t(x.0[y.0<0,])-mu1.0.h))
    S2 = (sum(y.0>0)-1)*stats::cov(t(t(x.0[y.0>0,])-mu2.0.h))
    Sig.X.0.h = (S1 + S2)/(n.vec[1]-2)
    result.bic = lassoshooting(XtX=Sig.X.0.h, Xty=(mu2.0.h-mu1.0.h), maxit = 1000, lambda=seq.lambda[i]*sqrt(2*2*log(p)/length(ind.1)))
    if(result.bic$iterations>=1000){
      err[i]=Inf
      break
    }else{
      beta.bic <- as.numeric(result.bic$coef)
      tll.bic = min(-t(beta.bic)%*%Sig.X.0.h%*%beta.bic + 2 * t(beta.bic)%*%(mu2.0.h-mu1.0.h) - sum(y.0^2/n.vec[1]),0) + sum(y.0^2/n.vec[1])
      k <- sum(as.numeric(beta.bic)!=0)
      n <- n.vec[1]
      BIC<-1*log(n)*k/n - tll.bic
      
      err[i] = BIC
    }
    cat("iter",i,"err",err[i],'\n')
  }
  if(sum(err!=Inf) == 0) return("max lambda too small")
  return(seq.lambda[which.min(err)])
}

# project vector to a simplex
projsimplx = function (y) {
  n = length(y)
  s = sort(y, decreasing = TRUE)
  tmpsum = cumsum(s[1:(n - 1)])
  tmp = (tmpsum - 1)/(1:(n - 1))
  ind = which(tmp >= s[2:n])[1]
  if (!is.na(ind)) {
    t = tmp[ind]
  }
  else {
    t = (tmpsum[n - 1] + s[n] - 1)/n
  }
  x = apply(cbind(y - t, 0), 1, max)
  return(x)
}

# Gradient descent method
grad.descent = function(f, x0, max.iter=200, step.size=0.05,
                        stopping.deriv=0.01, ...) {
  n = length(x0)
  xmat = matrix(0,nrow=n,ncol=max.iter)
  xmat[,1] = projsimplx(x0)
  for (k in 2:max.iter) {
    # Calculate the gradient
    grad.cur = grad(f,xmat[,k-1],...)
    if (all(abs(grad.cur) < stopping.deriv)) {
      k = k-1; break
    }
    # Move in the opposite direction of the grad
    xmat[,k] = projsimplx(xmat[,k-1] - step.size * grad.cur)
  }
  xmat = xmat[,1:k] # Trim
  return(list(x=xmat[,k], xmat=xmat, k=k))
}

# Coordinate descent method
coor.descent = function(f, x0, max.iter=200, step.size=0.05,
                        stopping.deriv=0.01, ...) {
  n = length(x0)
  xmat = matrix(0,nrow=n,ncol=max.iter)
  xmat[,1] = projsimplx(x0)
  i = 1
  for (k in 2:max.iter) {
    # Calculate the gradient
    grad.cur = grad(f,xmat[,k-1],...)
    if (all(abs(grad.cur) < stopping.deriv)) {
      k = k-1; break
    }
    # Move in the opposite direction of the grad
    xmat[,k] = xmat[,k-1]
    xmat[i,k] = xmat[i,k-1] - step.size * grad.cur[i]
    xmat[,k] = projsimplx(xmat[,k])
    if(i==n){
      i=1
    }else{
      i=i+1
    }
  }
  xmat = xmat[,1:k] # Trim
  return(list(x=xmat[,k], xmat=xmat, k=k))
}

# For evaluation: Estimation of mean
mu.est <- function(X, y, n.vec){
  ind.1<-1:n.vec[1]
  y.k<- y[ind.1]
  x.k<- X[ind.1,]
  pi1.0.h<- sum(y.k<0)/n.vec[1]
  pi2.0.h<- sum(y.k>0)/n.vec[1]
  mu1.0.h<- apply(x.k[y.k<0,],2,mean)
  mu2.0.h<- apply(x.k[y.k>0,],2,mean)
  S1 = (sum(y.k<0)-1)*stats::cov(t(t(x.k[y.k<0,])-mu1.0.h))
  S2 = (sum(y.k>0)-1)*stats::cov(t(t(x.k[y.k>0,])-mu2.0.h))
  Sig.X.0.h = (S1 + S2)/(n.vec[1]-2)
  return(list(mu1.0.h=mu1.0.h, mu2.0.h=mu2.0.h, pi1.0.h=pi1.0.h, pi2.0.h=pi2.0.h, Sig.X.0.h =  Sig.X.0.h))
}

# Get the set of indexes of different studies
ind.set<- function(n.vec, k.vec){
  ind.re <- NULL
  for(k in k.vec){
    if(k==1){
      ind.re<-c(ind.re,1: n.vec[1])
    }else{
      ind.re<- c(ind.re, (sum(n.vec[1:(k-1)])+1): sum(n.vec[1:k]))
    }
  }
  ind.re
}

# Repeat column
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}


   
  

