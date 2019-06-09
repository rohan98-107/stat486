#Assignment 8 - Interesting Hypotheses 

myfunc <- function(v,Q,dependent=F)
{
  pvals<-sort(v)
  m<-length(pvals)
  if(dependent)
    for(i in 1:m)
      qline<-Q*c(1:m)/(m*(sum(1/i)))
  else
    qline<-Q*c(1:m)/m
  plot(c(c(1:m),c(1:m)),c(qline,pvals),type="n",xlab="index",ylab="pvalue")
  lines(c(1:m),qline)
  points(c(1:m),pvals)
  temp_ind<-max(which(pvals<=qline))
  pstar<-pvals[temp_ind]
  indices<-(pvals<=pstar) #interesting indices
  res<-v[indices] #interesting values 
  points(pvals[indices],col="red")
}

p<-c(1e-6*runif(20),runif(400))
myfunc(p,0.95)