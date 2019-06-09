#ASSIGNMENT 4 
vectorizedBootstrap1<-
  function(vec0,nboot=1000,alpha=0.1)
  {
    #vec0 = c(1:10)
    #pretend this is our theoretical sample data 
    n<-length(vec0)
    mu0<-mean(vec0)
    sd0<-sd(vec0)
    
    #generate a vector of length 'nboot' that contains iterative means
    #vectorized implementation
    #Im changing 2 things, creating function that can be applied directly with apply
    # and changing the input matrix to be full length
    my.sample<-function(x){sample(x,replace=T)}
    bootdist0<-matrix(rep(vec0,nboot),byrow=F)
    #  bootdist<-apply(bootdist0,2,sample,x=vec0,size=nboot*n,replace=T)
    bootdist<-apply(bootdist0,2,my.sample)
    
    bootdist<-matrix(bootdist,nrow = n,ncol = nboot)
    print(bootdist[,1:10])
    #convert the sampled data into a studentized vector 
    helper <- function(x) ((mean(x)-mu0)/(sd(x)/sqrt(n)))
    bootvec<-apply(bootdist,2,helper)
    print(bootvec[1:100])
  #  print(mean(bootvec))
    
    #create bootstrapped CI
    lowerq<-quantile(bootvec,alpha/2,na.rm=T)
    upperq<-quantile(bootvec,1-alpha/2,na.rm=T)
    print(lowerq)
    print(upperq)
    
    LB<-mu0-(sd0/sqrt(n))*upperq #switch just like last project 
    UB<-mu0-(sd0/sqrt(n))*lowerq
    print(LB)
    print(UB)
    
    #create normal CI
    NLB<-mu0-(sd0/sqrt(n))*qnorm(1-alpha/2)
    NUB<-mu0+(sd0/sqrt(n))*qnorm(1-alpha/2)
    print(NLB)
    print(NUB)
    
    list(bootstrap.CI=c(LB,UB),normal.CI=c(NLB,NUB))
  }

vectorizedBootstrap <- function(vec0,nboot=10000,alpha=0.1)
{
  #vec0 = c(1:10)
  #pretend this is our theoretical sample data 
  n<-length(vec0)
  mu0<-mean(vec0)
  sd0<-sd(vec0)
  
  #generate a vector of length 'nboot' that contains iterative means
  #vectorized implementation
  bootdist<-matrix(vec0)
  bootdist<-apply(bootdist,2,sample,x=vec0,size=nboot*n,replace=T)
  bootdist<-matrix(bootdist,nrow = n,ncol = nboot)
  print(bootdist[,1:10])
  #convert the sampled data into a studentized vector 
  helper <- function(x) ((mean(x)-mu0)/(sd(x)/sqrt(n)))
  bootvec<-apply(bootdist,2,helper)
  print(bootvec[1:100])
  
  #create bootstrapped CI
  lowerq<-quantile(bootvec,alpha/2,na.rm=T)
  upperq<-quantile(bootvec,1-alpha/2,na.rm=T)
  print(lowerq)
  print(upperq)
  
  LB<-mu0-(sd0/sqrt(n))*upperq #switch just like last project 
  UB<-mu0-(sd0/sqrt(n))*lowerq
  print(LB)
  print(UB)
  
  #create normal CI
  NLB<-mu0-(sd0/sqrt(n))*qnorm(1-alpha/2)
  NUB<-mu0+(sd0/sqrt(n))*qnorm(1-alpha/2)
  print(NLB)
  print(NUB)
  
  list(bootstrap.CI=c(LB,UB),normal.CI=c(NLB,NUB))
}

simulate <- function(mu.theoretical=3,n=30,nsim=1000,alpha=0.1)
{
  #create coverage indicator vectors for both bootstrapped and normal dists.
  cvec.boot<-NULL
  cvec.norm<-NULL
  
  #real population mean 
  mu<-(exp(mu.theoretical+0.5))
  
  #simulate multiple bootstraps
  for (i in 1:nsim)
  {
    mu0vec<-rlnorm(n,mu.theoretical) #vector of sample means 
    tempobj<-vectorizedBootstrap1(mu0vec,alpha=alpha)
    boot.CI<-tempobj$bootstrap.CI
    norm.CI<-tempobj$normal.CI
    cvec.boot<-c(cvec.boot,(boot.CI[1]<mu)*(boot.CI[2]>mu))
    cvec.norm<-c(cvec.norm,(norm.CI[1]<mu)*(norm.CI[2]>mu))

  }
  
  #percentage results
  list(boot.coverage=(sum(cvec.boot)/nsim),norm.coverage=(sum(cvec.norm)/nsim))
}

print(simulate(n=15))
#something's wrong here 