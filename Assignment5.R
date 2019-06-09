#ASSIGNMENT 5 

leaps.pck<-c("leaps.pck", "leaps.then.press.plot", "regpluspress", "matrix.2ndorder.make")

leaps.then.press.plot<-function(xmat0,yvec,ncheck=4,print.ls=F, resid.plot=F)
  {
    n1<-ceiling(sqrt(ncheck)) # make n1 into an integer
    par(mfrow=c(n1,n1)) # create array of plots 
    if(resid.plot) # if residual plots are included, we want sqrt(2n)
    {
      n1<-ceiling(sqrt(2*ncheck))
    }
    par(mfrow=c(n1,n1))
    xmat<-matrix.2ndorder.make(xmat0) # create a matrix of only quadratic terms 
    #xmat<-xmat0
    leaps.str<-leaps(xmat,yvec) #apply leaps function - exhaustive search for best regression - variable selection
    z1<-leaps.str$Cp # assign cp values to a vector z1
    o1<-order(z1) # will use later to sort cp values 
    matwhich<-(leaps.str$which[o1,])[1:ncheck,] # create matrix of relevant data terms 
    z2<-z1[o1][1:ncheck] # sorted cp values in a vector 
    #zed<-list()
    for(i in 1:ncheck)
      {
      ls.str0<-regpluspress(xmat[,matwhich[i,]],yvec) # apply projection function 
      if(print.ls)
      {
        ls.print(ls.str0)
      }
      print(paste("Press=",ls.str0$press))
      parvec<-matwhich[i,] #relevant parameter vector 
      npar<-sum(parvec) 
      print(paste("MPSE=",ls.str0$press/(length(yvec)-(npar+1)))) #mean squared prediction error 
      MPSE<-floor(1000*ls.str0$press/(length(yvec)-(npar+1)))/1000 #exactifies the MPSE (truncates)
      print(paste("Cp-p=",z2[i]))
      ypred<-cbind(1,xmat[,matwhich[i,]])%*%ls.str0$coef #predictive vector
      #zed<-c(zed,list(coef=ls.str0$coef))
      
      plot(ypred,yvec,main=paste("I=",i,"MPSE=",MPSE,"Cp-p=",z2[i])) #plotting the regression 
      if(resid.plot)
      {
        plot(ypred,ls.str0$resid) #plotting residuals if tag = T
      }
    }
    #print(zed)
  }
regpluspress<-function(x,y) #function to project the regressions onto the 2d space 
  {
    ls.str<-lsfit(x,y) #create a least squares fit object 
    reg.influence<-1/(1-hat(x)) #influence factor 
    press<-sum((ls.str$resid/(1-hat(x)))^2) #press regression line onto the 2d space and store it in variable 
    ls.str$leverage<-hat(x) #diff between old and new press value 
    ls.str$press<-press #assign new press value 
    ls.str # return object 
  }
matrix.2ndorder.make<-function(x, only.quad=T)
  {
    x0<-x
    dimn<-dimnames(x)[[2]] #extract the names of the variables
    num.col<-length(x[1,]) # how many columns
    for(i in 1:num.col){
      # if we are doing all 2nd order - set to default 
      if(!only.quad){
        for(j in i:num.col){
          x0<-cbind(x0,x[,i]*x[,j]) # recursively add only the terms of degree 2 
          dimn<-c(dimn,paste(dimn[i],dimn[j],sep="")) # name each of these features 
          #create interaction dimnames
          
        }
      }
      else
      {
        #in here only if doing only squared terms
        x0<-cbind(x0,x[,i]*x[,i])
        dimn<-c(dimn,paste(dimn[i],"2",sep="")) # squared dimmension names
      }
    }
    dimnames(x0)[[2]]<-dimn
    x0
  }

library(leaps)
Auto.mat<-data.matrix(Auto)
#plot(Auto.mat)
x = 5 # Choose which column to model 
leaps.then.press.plot(Auto.mat[,-x],Auto.mat[,x],ncheck=5) #top 5 regressions 
