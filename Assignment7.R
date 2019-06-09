my.boot.xy.conf<-function(mat.train=auto.q.train,mat.test=auto.q.pred,yind=1,xstring="lars",brep=10000,pred.int=T,alpha=.05){
  #specialized version of bootstrap
  if(xstring=="lars"){
    library(lars)
    func0<-lars
    reduce.function<-my.cp.extract.lars1
  }
  if(xstring=="leaps"){
    library(leaps)
    func0<-leaps
    reduce.function<-my.cp.extract.leaps1
  }
  ypredmat<-NULL
  residmat<-NULL
  betamat<-NULL
  n1<-length(mat.train[,1])
  #
  out0<-reduce.function(func0(mat.train[,-yind],mat.train[,yind]),mat.train,mat.test,yind)
  ypred0<-c(out0$ypredict0)
  resid0<-c(out0$resid)
  beta0<-c(out0$beta)
  ypredmat<-rbind(ypredmat,ypred0)
  residmat<-rbind(residmat,resid0)
  betamat<-rbind(betamat,beta0)
  #
  for(i in 1:brep){
    if((i/200)==floor(i/200)){
      print(c(i,brep))
    }
    #
    v1<-sort(sample(n1,n1,replace=T))
    #
    m1<-(mat.train[,-yind])[v1,]
    y1<-(mat.train[,yind])[v1]
    #
    out1<-reduce.function(func0(m1,y1),mat.train,mat.test,yind)
    #
    ypred1<-c(out1$ypredict0)
    resid1<-c(out1$resid)
    beta1<-c(out1$beta)
    ypredmat<-rbind(ypredmat,ypred1)
    residmat<-rbind(residmat,resid1)
    betamat<-rbind(betamat,beta1)

  }
  #
  bagged.pred<-apply(ypredmat,2,mean)
  bagged.beta<-apply(betamat,2,mean)
  quant.boot<-function(x){quantile(x,c(alpha/2,1-alpha/2))}
  #
  if(pred.int){
    main1<-paste("Prediction interval",xstring,"alpha=",alpha)
    qboot<-apply(ypredmat+residmat,2,quant.boot)
  }
  else{
    main1=paste("Confidence interval,",xstring,"alpha=",alpha)
    qboot<-apply(ypredmat,2,quant.boot)
  }
  #
  y0<-mat.test[,yind]
  plot(rep(bagged.pred,5),c(y0,ypred0,bagged.pred,qboot[1,],qboot[2,]),xlab="Bagged prediction",ylab="Data and intervals",type="n",main=main1)
  points(bagged.pred,y0)
  lines(bagged.pred,bagged.pred)
  o1<-order(bagged.pred)
  lines(bagged.pred[o1],ypred0[o1],col=2)
  lines(bagged.pred[o1],smooth(qboot[1,o1]),col=3)
  lines(bagged.pred[o1],smooth(qboot[2,o1]),col=3)
  #
  list(bpred=bagged.pred,ypred0=ypred0,type=xstring,bagged.beta=bagged.beta,orig.beta=beta0,pred.int=pred.int)
}

my.cp.extract.lars1<-function(str,matrix.train=auto.q.train,matrix.test=auto.q.pred,yindex=1)
{

  I1<-(str$Cp==min(str$Cp))
  s1<-c(1:length(I1))[I1]
  xmat.train<-matrix.train[,-yindex]
  ymat.train<-matrix.train[,yindex]
  yp0<-predict(str,xmat.train,s1)$fit
  resid<-ymat.train-yp0
  xmat.pred<-matrix.test[,-yindex]
  yp1<-predict(str,xmat.pred,s1)$fit
  npred<-length(yp1)
  resp<-sample(resid,npred,replace=T)
  ypout<-yp1+resp
  list(ypredict0=yp1,resid=resp,beta=str$beta[I1,])

}

my.cp.extract.leaps1<- function(str,matrix.train=auto.q.train,matrix.test=auto.q.pred,yindex=1)
{
  I1<-(str$Cp==min(str$Cp))
  which1<-str$which[I1,]
  xmat.train<-(matrix.train[,-yindex])[,which1]
  ymat.train<-matrix.train[,yindex]
  ls.train<-lsfit(xmat.train,ymat.train)
  coef0<-ls.train$coef
  resid<-ls.train$resid
  xmat.pred<-(matrix.test[,-yindex])[,which1]
  yp1<-xmat.pred%*%coef0[-1]+coef0[1]
  npred<-length(yp1)
  resp<-sample(resid,npred,replace=T)
  list(ypredict0=yp1,resid=resp,beta=str$which)
  
}

auto.mat.q<-matrix.2ndorder.make(Auto.mat,T)
sampleindicies<-sample(392,292)
auto.q.train<-auto.mat.q[sampleindicies,]
auto.q.pred<-auto.mat.q[-sampleindicies,]
my.boot.xy.conf(xstring="leaps",brep=5000)