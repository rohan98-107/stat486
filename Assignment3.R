helper <- function(data,xindex,yindex,flag=F) #this is my own code 
{
  if(flag) {l1<-smooth.spline(data[[xindex]],sqrt(data[[yindex]])) }
  else {l1<-smooth.spline(data[[xindex]],data[[yindex]]) }
  x.resid<-resid(l1) #raw resids
  std.resid<-sqrt(sum(x.resid^2))/(length(data[[xindex]]) - l1$df) #standard resids 
  stud.resid<-x.resid/std.resid #studentized resids 
  D<-ks.test(stud.resid,pnorm)$statistic
  if(flag) {my.smooth<-sqrt(data[[yindex]])-x.resid }
  else {my.smooth<-data[[yindex]]-x.resid }
  list(D=D,raw.resid=x.resid,std.resid=std.resid,smooth=my.smooth)
}

myFunc <- function(data=NOAA,xindex=3,yindex=2,flag=F,nboot=1000,confidence=0.95)
{
  par(mfrow=c(1,1)) #plot only one graph 
  #conditional statement assigning sqrt to one of the variables if flag = T
  var1=data[[xindex]]
  if(flag)
    var2=sqrt(data[[yindex]])
  else
    var2=data[[yindex]]
  
  #initial values 
  init<-helper(data,xindex,yindex,flag)
 
  init.smooth<-init$smooth
  init.sdresid<-init$sd.resid
  init.resid<-init$raw.resid
  bootdata<-data
  n1<-length(init.smooth)
  smooth.dist<-NULL
  
  #bootstrapping
  for(i in 1:nboot)
  {
    temp<-sample(init.resid,length(init.resid),replace=T) #take a random sample of residuals with replacement
    boot.newdata<-((init.smooth+temp)) #create a temporary variable for bootstrapped data
    bootdata[[yindex]]<-boot.newdata #assign it 
    initp<-helper(bootdata,xindex,yindex) #iterate (new initial value)
    boot.smooth<-initp$smooth #assign new smoothing spline to boot.smooth 
    smooth.dist<-rbind(smooth.dist,boot.smooth-init.smooth) #appends the smooth data in data frame format and assigns it to smooth.data
  
    #this loop continuously updates the data that is being used for analysis by iterating random samples nboot times 
    #because the sample() function is pseudo random, this loop will generate new (or almost new) data sets every time the function is called 
    
  }
  
  #confidence interval
  alpha<-1-confidence
  LB<-NULL
  UB<-NULL
  for( i in 1:length(smooth.dist[1,]) )
  {
    #this loop approximates an interval at every x for y-alpha/2 < y < y+alpha/2
    s1<-sort(smooth.dist[,i])
    n2<-length(s1)
    v1<-c(1:n2)/n2
    bvec<-approx(v1,s1,c(alpha/2,1-alpha/2))$y
    LB<-c(LB,init.smooth[i]-bvec[2]) #bvec[2] > 0 thus every element of LB should be less than the initial spline
    UB<-c(UB,init.smooth[i]-bvec[1]) #bvec[1] < 0 thus every element of UB should be greater than the initial spline
    
    #provided the spline is a good fit - the lower bounds and upper bounds should never be intersected by the spline 
    
  }
  
  #plotting
  plot(rep(var1,4),c(LB,init.smooth,UB,var2),xlab="X",ylab="Y",type="n")
  points(var1,var2)
  o<-order(var1)
  lines(var1[o],LB[o],col=2)
  lines(var1[o],UB[o],col=2)
  lines(var1[o],init.smooth[o],col=3)	
  lines(smooth.spline(var1,var2,df=2),col=4)
}

NOAA<-read.csv("//Users//rohanrele//Downloads//NOAA+GISS.csv")
myFunc(flag=T)

