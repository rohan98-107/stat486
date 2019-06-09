myFunc <- function(dframe,var1_name,var2_name) {
  
  #extract column vectors 
  v1 = dframe[[which(names(dframe)==var1_name)]]
  v2 = dframe[[which(names(dframe)==var2_name)]]
  
  #print(v1)
  #print(v2)

  #plot smoothing splines 
  plot(v1,v2)
  l1<-smooth.spline(v1,v2)
  l2<-smooth.spline(v1,v2,df=2)
  #print(summary(l1))
  
  
  #difference of fit 
  l1.resid<-resid(l1)
  l2.resid<-resid(l2)
  SSF<-sum(l1.resid^2)
  SSN<-sum(l2.resid^2)
  F = ((SSN-SSF)/(l1$df-2))/(SSF/(length(v1)-2))
  p_val<-1-pf(F,floor((l1$df-2)),length(v1)-2)
  #print(p_val)
  #print(F)

  #plotting Q-Q graphs for resids
  par(mfrow=c(2,2))
  qqnorm(l1.resid,main="Smooth Residuals")
  qqnorm(l2.resid,main="Linear Residuals")
  plot(v1,v2)
  lines(l1)
  plot(v1,v2)
  lines(l2)
}

NOAA<-read.csv("//Users//rohanrele//Downloads//NOAA+GISS.csv")
#myFunc(NOAA,"delta.temp","X.disaster")
library(ISLR)
Auto <- read.table("http://www-bcf.usc.edu/~gareth/ISL/Auto.data", header=TRUE)
print(Auto)
myFunc(Auto,"horsepower","weight")
#myFunc(NOAA,"delta.temp","X.disaster")

