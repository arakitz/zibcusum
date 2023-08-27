### This code evaluates the OoC performace of the upper-sided ZIB-CUSUM control chart
require(VGAM)
###############
ZIBCUSUMout<-function(ds,H,h,n,p0,theta){
# ds is the shift in IC proportion p0
# H is the decision interval of the CUSUM chart
# h is the reference value of the CUSUM chart
# n is the sample size
# p0 is the in-control proportion
# theta is equal to 1-phi
###################
## OoC proportion
  p<-ds*p0;
## define phi
  pstr0<-1-theta;
  m<-100;
  tn<-round(m*H+1);
# Creating the transition probabilities matrix
  Q<-matrix(0,ncol=tn,nrow=tn);
  for(i in 1:tn){
    x<-(i-1)/100;
    k<-0;
    while(TRUE){
      j0<-max(0,x+k-h)
      j<-round(100*j0+1)
      if(j>tn){break}
      else{
        Q[i,j]<-Q[i,j]+dzibinom(k,size=n,prob=p,pstr0,log=FALSE);
        k<-k+1;
      }
    }
  }
# Identity Matrix
  I<-diag(1,tn);
# initial probabilities vector
  a1<-array(0,tn);
  a1[1]<-1
# a vector with 1s
  l1<-array(1,tn)
# Inverse of the (I-Q) matrix
  W<-solve(I-Q)
  nu1<-a1%*%W%*%l1
# zero-state ARL
  zsarl<-nu1
  nu2<-2*a1%*%W%*%W%*%Q%*%l1
# zero-state SDRL
  SDRL<-sqrt(nu2-(nu1^2)+nu1)
#  ssarl<-qss%*%solve(ID-Q)%*%l1  
  cat("ds:",ds,"  ARL:",zsarl,"  SDRL:",SDRL,"\n")
}
##########################################################
## An example for various upward shifts in p0
for(ds in c(1.0,1.1,1.2,1.3,1.5,1.7,2.0,2.5,3.0)){
  ZIBCUSUMout(ds,H=12.87,h=0.68,n=250,p0=0.01,theta=0.2)
}
########################## END ################################