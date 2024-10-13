#EM algorithm for GMR model - equal variance
EMGuasMix<- function(Y,alp_int,X,beta_int,sig_int) { 
  #alp,beta,sig start with initial values
  alp<-alp_int
  beta<-beta_int
  sig<-sig_int
  n<-length(Y)
  g<-length(alp)
  nbeta<-nrow(beta) #number of beta
  mu=X%*%beta
  pifx<- matrix(0, nrow=n, ncol=g) 
  for(j in 1:g ){
    pifx[,j] <- alp[j]*dnorm(Y, mean = mu[,j], sd = sig) 
  }
  
  loglik<-sum(apply(pifx,1,function(x) log(sum(x)))) #l(1)
  llogliks<-c(loglik)
  delta=loglik
  while(abs(delta) > 1e-5) {# stopping criteria
    Z<-t(apply(pifx, 1, function(x) x / sum(x)))
    alp<-apply(Z,2,function(x) sum(x)/n)
    beta<- matrix(0, nrow=nbeta, ncol=g) 
    for(j in 1:g ){
      beta[,j] <- (lm(Y~X[,-1],weight=Z[,j]))$coefficients
    }
    mu=X%*%beta
    si<- matrix(0, nrow=n, ncol=g) 
    for(j in 1:g ){
      si[,j] <- (Y-mu[,j])^2*Z[,j]
    }
    sig02<-sum(si)/n
    sig<-sqrt(sig02)
    pifx<- matrix(0, nrow=n, ncol=g) 
    for(j in 1:g ){
      pifx[,j] <- alp[j]*dnorm(Y, mean = mu[,j], sd = sig) 
    }
    loglik<-sum(apply(pifx,1,function(x) log(sum(x)))) #l(k)
    llogliks<-c(llogliks,loglik)
    delta<- llogliks[length(llogliks)]  - llogliks[length(llogliks)-1]}
  #return(alp, beta, sigma, group, incomplete-loglikilihood)
  return(em=list(alp,beta,sig,apply(Z,1,function(x) (which(x==max(x)))),llogliks))#sig02
}
#EM algorithm for GMR model - unequal variance
EMGuasMix_v<- function(Y,alp_int,X,beta_int,sig_int) {
  alp<-alp_int
  beta<-beta_int
  sig<-sig_int
  n<-length(Y)
  g<-length(alp)
  nbeta<-nrow(beta) #number of beta
  Z<-t(rmultinom(n,1,alp))
  mu=X%*%beta
  pifx<- matrix(0, nrow=n, ncol=g) 
  
  for(j in 1:g ){
    pifx[,j] <- alp[j]*dnorm(Y, mean = mu[,j], sd = sig[j]) 
  }
  
  loglik<-sum(apply(pifx,1,function(x) log(sum(x)))) #l(1)
  llogliks<-c(loglik)
  delta=loglik
  while(abs(delta) > 1e-5) {
    Z<-t(apply(pifx, 1, function(x) x / sum(x)))
    alp<-apply(Z,2,function(x) sum(x)/n)
    beta<- matrix(0, nrow=nbeta, ncol=g) 
    for(j in 1:g ){
      beta[,j] <- (lm(Y~X[,-1],weight=Z[,j]))$coefficients
    }
    mu=X%*%beta
    si<- matrix(0, nrow=n, ncol=g) 
    for(j in 1:g ){
      si[,j] <- (Y-mu[,j])^2*Z[,j]
    }
    sig02<-apply(si,2,sum)/apply(Z,2,sum)
    sig<-sqrt(sig02)
    pifx<- matrix(0, nrow=n, ncol=g) 
    for(j in 1:g ){
      pifx[,j] <- alp[j]*dnorm(Y, mean = mu[,j], sd = sig[j]) 
    }
    loglik<-sum(apply(pifx,1,function(x) log(sum(x)))) #l(k)
    llogliks<-c(llogliks,loglik)
    delta<- llogliks[length(llogliks)]  - llogliks[length(llogliks)-1]}
  return(em=list(alp,beta,sig,apply(Z,1,function(x) (which(x==max(x)))),llogliks))#sig02
}
#EM algorithm for Homo-GMR model - equal variance
EMGuasMix_fix1<- function(Y,alp_int,X,beta_int,sig_int) {
  alp<-alp_int
  beta<-beta_int
  beta0<-beta[1,]
  beta1<-beta[2,1]
  sig<-sig_int
  n<-length(Y)
  g<-length(alp)
  nbeta<-nrow(beta) #number of beta
  Z<-t(rmultinom(n,1,alp))
  mu=X%*%beta
  pifx<- matrix(0, nrow=n, ncol=g) 
  
  for(j in 1:g ){
    pifx[,j] <- alp[j]*dnorm(Y, mean = mu[,j], sd = sig) 
  }
  
  loglik<-sum(apply(pifx,1,function(x) log(sum(x)))) #l(1)
  llogliks<-c(loglik)
  delta=loglik
  while(abs(delta) > 1e-5) {
    Z<-t(apply(pifx, 1, function(x) x / sum(x)))
    alp<-apply(Z,2,function(x) sum(x)/n)
    ybar<-apply(Z,2,function(t) sum(t*Y)/sum(t))
    xbar<-apply(Z,2,function(t) sum(t*X[,-1])/sum(t))
    ydiffZ<-outer(Y,ybar,'-')*Z
    xdiffZ<-outer(X[,-1],xbar,'-')*Z
    beta1<-sum(apply(ydiffZ,2,function(t) t*X[,-1]))/sum(apply(xdiffZ,2,function(t) t*X[,-1]))
    beta0<-ybar-beta1*xbar
    beta<- rbind(beta0,beta1)
    mu=X%*%beta
    si<- matrix(0, nrow=n, ncol=g) 
    for(j in 1:g ){
      si[,j] <- (Y-mu[,j])^2*Z[,j]
    }
    sig02<-sum(si)/n
    sig<-sqrt(sig02)
    pifx<- matrix(0, nrow=n, ncol=g) 
    for(j in 1:g ){
      pifx[,j] <- alp[j]*dnorm(Y, mean = mu[,j], sd = sig) 
    }
    loglik<-sum(apply(pifx,1,function(x) log(sum(x)))) #l(k)
    llogliks<-c(llogliks,loglik)
    delta<- llogliks[length(llogliks)]  - llogliks[length(llogliks)-1]}
  return(em=list(alp,beta,sig,apply(Z,1,function(x) (which(x==max(x)))),llogliks))#sig02
}

#EM algorithm for Homo-GMR model - unequal variance
EMGuasMix_fix1_v<- function(Y,alp_int,X,beta_int,sig_int) {
  alp<-alp_int
  beta<-beta_int
  beta0<-beta[1,]
  beta1<-beta[2,1]
  sig<-sig_int
  n<-length(Y)
  g<-length(alp)
  nbeta<-nrow(beta) #number of beta
  #logliks<-c()
  #llc<- matrix(0, nrow=n, ncol=g) #complete data loglikelihood
  Z<-t(rmultinom(n,1,alp))
  mu=X%*%beta
  pifx<- matrix(0, nrow=n, ncol=g) 
  
  for(j in 1:g ){
    pifx[,j] <- alp[j]*dnorm(Y, mean = mu[,j], sd = sig[j]) 
  }
  
  loglik<-sum(apply(pifx,1,function(x) log(sum(x)))) #l(1)
  llogliks<-c(loglik)
  delta=loglik
  while(abs(delta) > 1e-5) {
    Z<-t(apply(pifx, 1, function(x) x / sum(x)))
    alp<-apply(Z,2,function(x) sum(x)/n)
    ybar<-apply(Z,2,function(t) sum(t*Y)/sum(t))
    xbar<-apply(Z,2,function(t) sum(t*X[,-1])/sum(t))
    ydiffZ<-outer(Y,ybar,'-')*Z
    xdiffZ<-outer(X[,-1],xbar,'-')*Z
    beta1<-sum(apply(ydiffZ,2,function(t) t*X[,-1]))/sum(apply(xdiffZ,2,function(t) t*X[,-1]))
    beta0<-ybar-beta1*xbar
    beta<- rbind(beta0,beta1)
    mu=X%*%beta
    si<- matrix(0, nrow=n, ncol=g) 
    for(j in 1:g ){
      si[,j] <- (Y-mu[,j])^2*Z[,j]
    }
    sig02<-apply(si,2,sum)/apply(Z,2,sum)
    sig<-sqrt(sig02)
    pifx<- matrix(0, nrow=n, ncol=g) 
    for(j in 1:g ){
      pifx[,j] <- alp[j]*dnorm(Y, mean = mu[,j], sd = sig[j]) 
    }
    loglik<-sum(apply(pifx,1,function(x) log(sum(x)))) #l(k)
    llogliks<-c(llogliks,loglik)
    delta<- llogliks[length(llogliks)]  - llogliks[length(llogliks)-1]}
  return(em=list(alp,beta,sig,apply(Z,1,function(x) (which(x==max(x)))),llogliks))#sig02
}

#EM algorithm for GM model - equal variance
EMGuasMix_NoX<- function(Y,alp_int,mu_int,sig_int) {
  alp<-alp_int
  mu<-mu_int
  sig<-sig_int
  n<-length(Y)
  g<-length(alp)
  Z<-t(rmultinom(n,1,alp))
  pifx<- matrix(0, nrow=n, ncol=g) 
  for(j in 1:g ){
    pifx[,j] <- alp[j]*dnorm(Y, mean = mu[j], sd = sig) 
  }
  loglik<-sum(apply(pifx,1,function(x) log(sum(x)))) #l(1)
  llogliks<-c(loglik)
  delta=loglik
  while(abs(delta) > 1e-5)   {
    Z<-t(apply(pifx, 1, function(x) x / sum(x)))
    alp<-apply(Z,2,function(x) sum(x)/n)
    mu<-apply(Z,2,function(x) sum(x*Y)/sum(x))
    si<- matrix(0, nrow=n, ncol=g) 
    for(j in 1:g){
      si[,j] <- (Y-mu[j])^2*Z[,j]
    }
    sig02<-sum(si)/n
    sig<-sqrt(sig02)
    pifx<- matrix(0, nrow=n, ncol=g) 
    for(j in 1:g ){
      pifx[,j] <- alp[j]*dnorm(Y, mean = mu[j], sd = sig) 
    }
    loglik<-sum(apply(pifx,1,function(x) log(sum(x))))#l(k)
    llogliks<-c(llogliks,loglik)
    delta<- llogliks[length(llogliks)]  - llogliks[length(llogliks)-1]}
  return(em=list(alp,mu,sig,apply(Z,1,function(x) (which(x==max(x)))),llogliks))#sig2
}

#EM algorithm for GM model - unequal variance
EMGuasMix_NoX_v<- function(Y,alp_int,mu_int,sig_int) {
  alp<-alp_int
  mu<-mu_int
  sig<-sig_int
  n<-length(Y)
  g<-length(alp)
  Z<-t(rmultinom(n,1,alp))
  pifx<- matrix(0, nrow=n, ncol=g) 
  for(j in 1:g ){
    pifx[,j] <- alp[j]*dnorm(Y, mean = mu[j], sd = sig[j]) 
  }
  loglik<-sum(apply(pifx,1,function(x) log(sum(x)))) #l(1)
  llogliks<-c(loglik)
  delta=loglik
  while(abs(delta) > 1e-5)   {
    Z<-t(apply(pifx, 1, function(x) x / sum(x)))
    alp<-apply(Z,2,function(x) sum(x)/n)
    mu<-apply(Z,2,function(x) sum(x*Y)/sum(x))
    si<- matrix(0, nrow=n, ncol=g) 
    for(j in 1:g){
      si[,j] <- (Y-mu[j])^2*Z[,j]
    }
    sig02<-apply(si,2,sum)/apply(Z,2,sum)
    sig<-sqrt(sig02)
    pifx<- matrix(0, nrow=n, ncol=g) 
    for(j in 1:g ){
      pifx[,j] <- alp[j]*dnorm(Y, mean = mu[j], sd = sig[j]) 
    }
    loglik<-sum(apply(pifx,1,function(x) log(sum(x))))#l(k)
    llogliks<-c(llogliks,loglik)
    delta<- llogliks[length(llogliks)]  - llogliks[length(llogliks)-1]}
  return(em=list(alp,mu,sig,apply(Z,1,function(x) (which(x==max(x)))),llogliks))#sig2
}



#Likelihood ratio test for overall test when G is known
LRT<-function(MG,MGN){#MG-GMR model, MGN-GM model
  n<-length(unlist(MG[4]))
  g<-length(unlist(MG[1]))
  nbeta<-length(unlist(MG[2]))/g
  niter<-length(unlist(MG[5]))
  niterN<-length(unlist(MGN[5]))
  loglikic<-unlist(MG[5])
  loglikicN<-unlist(MGN[5])
  Lrt<-2*(loglikic[niter]-loglikicN[niterN])
  pvalue<-pchisq(Lrt,(nbeta-1)*g,lower=F)
  return(list(loglikic[niter],loglikicN[niterN],Lrt,pvalue)) 
}

#Likelihood ratio test for heterogenous test when G is known
LRT_fix<-function(MG,MGN){ #MG-GMR model, MGN-HomoGMR model
  n<-length(unlist(MG[4]))
  g<-length(unlist(MG[1]))
  nbeta<-length(unlist(MG[2]))/g
  niter<-length(unlist(MG[5]))
  niterN<-length(unlist(MGN[5]))
  loglikic<-unlist(MG[5])
  loglikicN<-unlist(MGN[5])
  Lrt<-2*(loglikic[niter]-loglikicN[niterN])
  pvalue<-pchisq(Lrt,(g-1),lower=F)
  return(list(loglikic[niter],loglikicN[niterN],Lrt,pvalue))
}

#Calculate BIC for GMR/GM model - equal variance
BICem<- function(MG) {
  n<-length(unlist(MG[4]))
  g<-length(unlist(MG[1]))
  nbeta<-length(unlist(MG[2]))/g
  niter<-length(unlist(MG[5]))
  loglikic<-unlist(MG[5])
  if (nbeta==1)
  {BICEM<--2*loglikic[niter]+g*2*log(n)}
  else
  {BICEM<--2*loglikic[niter]+g*(nbeta+1)*log(n)}
  return(BICEM)}

#Calculate BIC for GMR/GM model - unequal variance
BICem_v<- function(MG) {
  n<-length(unlist(MG[4]))
  g<-length(unlist(MG[1]))
  nbeta<-length(unlist(MG[2]))/g
  niter<-length(unlist(MG[5]))
  loglikic<-unlist(MG[5])
  if (nbeta==1)
  {BICEM<--2*loglikic[niter]+(g*3-1)*log(n)}
  else
  {BICEM<--2*loglikic[niter]+(g*(nbeta+2)-1)*log(n)}
  return(BICEM)}

#Calculate BIC for Homo-GMR model - equal variance
BICem_fix<- function(MG) {
  n<-length(unlist(MG[4]))
  g<-length(unlist(MG[1]))
  nbeta<-length(unlist(MG[2]))/g
  niter<-length(unlist(MG[5]))
  loglikic<-unlist(MG[5])
  BICEM<--2*loglikic[niter]+(2*g+1)*log(n)
  return(BICEM)}

#Calculate BIC for Homo-GMR model - unequal variance
BICem_fix_v<- function(MG) {
  n<-length(unlist(MG[4]))
  g<-length(unlist(MG[1]))
  nbeta<-length(unlist(MG[2]))/g
  niter<-length(unlist(MG[5]))
  loglikic<-unlist(MG[5])
  BICEM<--2*loglikic[niter]+(3*g)*log(n)
  return(BICEM)}


