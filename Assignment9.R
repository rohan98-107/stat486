my.cp.extract.lars0<-function(str,matrix.train,matrix.test,yindex=1)
{
  I1<-(str$Cp==min(str$Cp))
  s1<-c(1:length(I1))[I1]
  s1.ls<-length(I1)
  xmat.train<-matrix.train[,-yindex]
  ymat.train<-matrix.train[,yindex]
  yp0<-predict(str,xmat.train,s1)$fit
  resid<-ymat.train-yp0
  xmat.pred<-matrix.test[,-yindex]
  yp1<-predict(str,xmat.pred,s1)$fit
  yp2<-predict(str,xmat.pred,s1.ls)$fit
  list(ypredict0=yp1,ypredict.full=yp2)
  
}
my.chaos.reg.test<-function(nsel,ncol=20,trim0=.1,reps=10){
  library(lars)
  library(caret)
  Lorenz2<-matlag1(lorenz1,c(50,100,150,200,250,300))
  D0<-disjoint.delaymap.make1(reps,ncol,nsel)
  y<-Lorenz2[,21]
  D1a<-Lorenz2[c(1:3000),]
  D1b<-Lorenz2[-c(1:3000),]
  Y1a<-y[1:3000]
  Y2a<-y[-c(1:3000)]
  PM1<-NULL
  PM2<-NULL
  PM3<-NULL
  for(i in 1:length(D0)){
    M1a<-D1a[,D0[[i]]]
    M2a<-D1b[,D0[[i]]]
    Lars.str<-lars(M1a,Y1a)
    Pred.out<-my.cp.extract.lars0(Lars.str,cbind(Y1a,M1a),cbind(Y2a,M2a))
    PM1<-rbind(PM1,Pred.out$ypredict0)
    PM2<-rbind(PM2,Pred.out$ypredict.full)
  }
  cv.lars(M1a,Y1a)
  tmean<-function(x){mean(x,trim=trim0)}
  cor1<-cor(apply(PM1,2,tmean),Y2a)
  cor2<-cor(apply(PM2,2,tmean),Y2a)
  t1<-cor1/sqrt((1-cor1^2)/(length(Y2a)-2))
  p1<-1-pt(t1,length(Y2a)-2)
  t2<-cor1/sqrt((1-cor2^2)/(length(Y2a)-2))
  p2<-1-pt(t2,length(Y2a)-2)
  list(lars=list(cor=cor1,p=p1),ls=list(cor=cor2,p=p2))
}

matlag0<-function(mat,n)
{
  n2<-length(mat[,1])
  nvec<-c((n2-n+1):n2)
  matA<-mat[-c(1:n),]
  matB<-mat[-nvec,]
  matC<-cbind(matB,matA)
  matC
}

matlag1<-function(mat,nvec)
{
  n3<-length(nvec)
  for(i in 1:n3){
    j<-nvec[i]
    m2<-matlag0(mat,j)
    if(i==1){
      m3<-m2
    }
    else{
      n1<-length(m2[1,])
      m2a<-m2[,-c(1:(n1/2))]
      n2<-length(m2a[,1])
      m3<-cbind(m3[1:n2,],m2a)	
    }
  }
  m3
}

disjoint.delaymap.make1<-function(ntrial, ncol, nsel)
{
  nmod<-ntrial*ceiling(ncol/nsel)
  nmid<-ceiling(ncol/nsel)
  outlist <- rep(list(), nmod)
  k<-1
  for(i in 1:ntrial) {
    dum <- sample(c(1:ncol),  replace = F)
    dum0<-dum
    for(j in 1:nmid){
      if(j<nmid){
        selvec<-((j-1)*nsel)+c(1:nsel)
        outlist[[k]] <- dum[selvec]
        dum0<-dum0[-c(1:nsel)]
        k<-k+1
      }
      else{
        if(length(dum0)<nsel){
          nplus<-nsel-length(dum0)
          dum0<-c(dum0,sample(c(1:ncol)[-dum0],nplus,replace=F))
        }
        outlist[[k]]<-dum0
        k<-k+1
      }
    }
  }
  outlist
}

temp<-read.csv("//Users//rohanrele//Downloads//lorenz1.csv")
lorenz1<-as.matrix(temp[-1])
print(my.chaos.reg.test(2))