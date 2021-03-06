options(encoding="UTF-8")
                                        #
#
# source("nlm_testing_1.R")
#
#library(ggplot)
library(optimx)
library(flexsurv)


flincorr <- function(x,para) {
    a <- para[1] ; b <- para[2] ; c <- para[3] ; d <- para[4] ; h <- para[5] ;
    l <- para[6] ;

    log(b) + a*(x + c*((x-h)^3+3*(x-h))/(l+exp(d*abs(x-h))))


    }


f <- function(x) {
  out <- -exp(-0.5 * x^2)
  attr(out, 'gradient') <- -x * out 
  attr(out, 'hessian') <-  (x^2 - 1) * out
  return(out)
}

nlm(f, 1.3, hessian = TRUE, check.analyticals = TRUE)


fpara <- function(para) {
  x <- para[1] ;  
  out <- -exp(-0.5 * x^2)
  attr(out, 'gradient') <- -x * out 
  attr(out, 'hessian') <-  (x^2 - 1) * out
  return(out)
}

#nlm(f, 1.3, hessian = TRUE, check.analyticals = TRUE)

fpara2 <- function(para) {
  x <- para[1] ; y <- para[2] ; 
  out <- -exp(-0.5 * (x^2 + y^2))
  attr(out, 'gradient') <- c(-x * out, -y * out) 
 # attr(out, 'hessian') <-  (x^2 - 1) * out
  return(out)
}

# nlm(fpara2, c(1.3,2), hessian = TRUE, check.analyticals = TRUE)

fpara3 <- function(para) {
  x <- para[1] ; y <- para[2] ; z <- para[3] ; 
  out <- exp(-0.5 * (x^2 + y^2 + z^2))
  attr(out, 'gradient') <- c(-x * out, -y * out, -z * out) 
 # attr(out, 'hessian') <-  (x^2 - 1) * out
  return(out)
}

# nlm(fpara3, c(1.3,2,10), hessian = TRUE, check.analyticals = TRUE)

fGompAugm5p0  <- function(x,para) {
    a <- para[1] ; b <- para[2] ;   c <- para[3] ; d <- para[4] ; h <- para[5] ; xepsi <- 1e-7 ;

    hxfunc <-(c*(x-h)^3 + xepsi*(x-h)) * exp(-d*(x-h)^2) ;
    ldxfunc <- 1 + (-2*d*(x-h)*(c*(x-h)^3 + xepsi*(x-h)) + 3*c*(x-h)^2 + xepsi)*exp(-d*(x-h)^2) ;
    b*ldxfunc*exp(a*(x + hxfunc))*exp(-b/a*exp(a*(x+hxfunc)) + b/a)
   
}

fGompMak  <- function(x,para,pPrt=0) {
    a <- para[1] ; b <- para[2] ;  lmbd <- para[3] ; xepsi <- 2 ; l <- 1 ;
    (b*exp(a*x) + lmbd)*exp(-lmbd*x  - b/a*exp(a*x) + b/a)
}

fGomp3p  <- function(x,para,pPrt=0) {
    a <- para[1] ; b <- para[2] ;  k <- para[3] ; xepsi <- 2 ; l <- 1 ;
    b*k*exp(a*x^k)*x^(k-1)*exp( - b/a*exp(a*x^k) + b/a)
}

fGompMak4p  <- function(x,para,pPrt=0) {
    a <- para[1] ; b <- para[2] ;  lmbd <- para[3] ; k <- para[4] ; xepsi <- 2 ; l <- 1 ;
    (b*k*exp(a*x)*x^(k-1) + lmbd)*exp(-lmbd*x  - b/a*exp(a*x^k) + b/a)
}


fGompAugm5p  <- function(x,para,pPrt=0) {
    a <- para[1] ; b <- para[2] ;   c <- para[3] ; d <- para[4] ; h <- para[5] ; xepsi <- 2 ; l <- 1 ;
    diffact <- ifelse(x<h,-1,1) ;
    
    hxfunc <-c*((x-h)^3 + 3*xepsi*(x-h))/(l + exp(d*abs(x-h))) ; 
    ldxfunc <- 1 + 3*c*((x-h)^2 + 1*xepsi)/(l + exp(d*abs(x-h))) + -d*c*diffact*((x-h)^3 + 3*xepsi*(x-h))*exp(d*abs(x-h))/(l + exp(d*abs(x-h)))^2 ;
    if (pPrt==1)  print(paste('--hxfunc: ',hxfunc,' --ldxfunc: ',ldxfunc)) ;
    b*ldxfunc*exp(a*(x + hxfunc))*exp(-b/a*exp(a*(x+hxfunc)) + b/a)
   
}

fGompMakAugm6p  <- function(x,para,pPrt=0) {
    a <- para[1] ; b <- para[2] ;   c <- para[3] ; d <- para[4] ; h <- para[5] ; lmbd <- para[6] ; xepsi <- 2 ; l <- 1 ;
    diffact <- ifelse(x<h,-1,1) ;
    
    hxfunc <-c*((x-h)^3 + 3*xepsi*(x-h))/(l + exp(d*abs(x-h))) ;  #print(paste('--hxfunc: ',hxfunc)) ;
    ldxfunc <- 1 + 3*c*((x-h)^2 + 1*xepsi)/(l + exp(d*abs(x-h))) + -d*c*diffact*((x-h)^3 + 3*xepsi*(x-h))*exp(d*abs(x-h))/(l + exp(d*abs(x-h)))^2 ;
    if (pPrt==1)  print(paste('--hxfunc: ',hxfunc,' --ldxfunc: ',ldxfunc)) ;
   # ldxfunc <- 1 ;
    (b*ldxfunc*exp(a*(x + hxfunc)) + lmbd)*exp(-lmbd*x  - b/a*exp(a*(x+hxfunc)) + b/a)
   
}



lfGompAugm5p  <- function(x,para) {
    a <- para[1] ; b <- para[2] ;   c <- para[3] ; d <- para[4] ; h <- para[5] ;

    log(b) + log(1 +( -2*d*(x-h)*(c*(x-h)^3 + (x-h))+(3*c*(x-h)^2 + 1))*exp(-d*(x-h)^2)) 
    -(b/a*exp(a*( x + (c*(x-h)^3 + (x-h)) * exp(-d*(x-h)^2))))  
    + a*(x + (c*(x-h)^3 + (x-h))*exp(-d*(x-h)^2)) + b/a ;
    ;
}

minusLogLikGompAugm5p0  <- function(para) {

    a <- para[1] ; b <- para[2] ;   c <- para[3] ; d <- para[4] ; h <- para[5] ;
    n <- sum(dx) ; xepsi <- 1e-2 ;

    hxfunc <-(c*(x-h)^3 + xepsi*(x-h)) * exp(-d*(x-h)^2)
    ldxfunc <- 1 + (-2*d*(x-h)*(c*(x-h)^3 + xepsi*(x-h)) + 3*c*(x-h)^2 + xepsi)*exp(-d*(x-h)^2) ;
    
    logdxsum <- sum(dx*log(ldxfunc)) ;      
    xsum <- sum(dx*(x + hxfunc)) ; expsum <- sum(dx*exp(a*(x + hxfunc))) ;

 - ( n*log(b) + n*b/a + logdxsum  + a*xsum  - b/a*expsum  )
    
}

minusLogLikGompAugm5p  <- function(para) {

    a <- para[1] ; b <- para[2] ;   c <- para[3] ; d <- para[4] ; h <- para[5] ;
    n <- sum(dx) ; xepsi <- 2 ; l <- 1 ;

    hxfunc <-c*((x-h)^3 + 3*xepsi*(x-h)) /(l+ exp(d*abs(x-h)))
    diffact <- ifelse(x<h,-1,1)
    ldxfunc <- 1 +  3*l*c*((x-h)^2 + 1*xepsi)/ (l + exp(d*abs(x-h))) -d*c*diffact*((x-h)^3 + 3*xepsi*(x-h))*exp(d*abs(x-h))/(l + exp(d*abs(x-h)))^2 ;
    
    logdxsum <- sum(dx*log(ldxfunc)) ;      
    xsum <- sum(dx*(x + hxfunc)) ; expsum <- sum(dx*exp(a*(x + hxfunc))) ;

 - ( n*log(b) + n*b/a + logdxsum  + a*xsum  - b/a*expsum  )
    
}

minusLogLikGompMakAugm6p  <- function(para) {

    a <- para[1] ; b <- para[2] ;   c <- para[3] ; d <- para[4] ; h <- para[5] ; lmbd <- para[6] ;
    n <- sum(dx) ; xepsi <- 2 ; l <- 1 ;

    hxfunc <-c*((x-h)^3 + 3*xepsi*(x-h)) /(l+ exp(d*abs(x-h)))
    diffact <- ifelse(x<h,-1,1)
    ldxfunc <- 1 +  3*l*c*((x-h)^2 + 1*xepsi)/ (l + exp(d*abs(x-h))) -d*c*diffact*((x-h)^3 + 3*xepsi*(x-h))*exp(d*abs(x-h))/(l + exp(d*abs(x-h)))^2 ;
    
    logdxsum <- sum(dx*log(ldxfunc)) ;      
    xsum <- sum(dx*x) ; expsum <- sum(dx*exp(a*(x + hxfunc))) ;

    logexpplusl <- log(b*exp(a*(x + hxfunc))*ldxfunc + lmbd) ; logexppluslsum <- sum(dx*logexpplusl) ;
    -(logexppluslsum - lmbd*xsum + ((-b)/a)*expsum + n*b/a);
  

    
}



minusLogLikGompAugm6p  <- function(para) {

    a <- para[1] ; b <- para[2] ;   c <- para[3] ; d <- para[4] ;  e <- para[5] ; h <- para[6] ;
    n <- sum(dx) ;

    hxfunc <-(c*(x-h)^3 + d*(x-h)) * exp(-e*(x-h)^2)
    ldxfunc <- 1 + (-2*e*(x-h)*(c*(x-h)^3 + d*(x-h)) + 3*c*(x-h)^2 + d)*exp(-e*(x-h)^2) ;
    
    logdxsum <- sum(dx*log(ldxfunc)) ;      
    xsum <- sum(dx*(x + hxfunc)) ; expsum <- sum(dx*exp(a*(x + hxfunc))) ;

 - ( n*log(b) + n*b/a + logdxsum  + a*xsum  - b/a*expsum  )
    
}



minusLogLikWeib <- function(para) {
    b <- para[1] ; k <- para[2] ;  n <- sum(dx) ;
    logsum <- sum(dx*log(x)) ; xksum <- sum(dx*x^k) ; xklogxsum <-  sum(dx*x^k*log(x)) ;
    out <- -( n*log(k) + n*k*log(b) + (k-1)*logsum  - b^k*xksum )
#    attr(out, 'gradient') <- c(-(n * k / b - k * b^(k-1) * xksum) ,
#       -( n/k + n * log(b) + logsum - b^k * log(b) * xksum - b^k * xklogxsum) ) 
    return(out)
}


minusLogLikGomp <- function(para) {
    b <- para[2] ; a <- para[1] ; n <- sum(dx) ;
    xsum <- sum(dx*x) ; expsum <- sum(dx*exp(a*x)) ;
    out <-  - (  n*log(b) + n*b/a + a*xsum  - b/a*expsum )
    return(out)
}



minusLogLikGompMak <- function(para) {

    b <- para[2] ; a <- para[1] ; l <- para[3] ;  n <- sum(dx) ;
    xsum <- sum(dx*x) ; expsum <- sum(dx*exp(a*x)) ;
    logexpplusl <- log(b*exp(a*x)+l) ; logexppluslsum <- sum(dx*logexpplusl) ;
    out <- -(logexppluslsum - l*xsum+((-b)/a)*expsum + n*b/a);
  
    return(out) 
}

minusLogLikGomp3p <- function(para) {
    b <- para[2] ; a <- para[1] ; k <- para[3] ; n <- sum(dx) ;
    x <- ifelse(x==0,1,x) ; 
    xksum <- sum(dx*x^k) ; expksum <- sum(dx*exp(a*x^k)) ;
    logxsum <-sum(dx*log(x)) ;
    out <- -(n*log(b)+ n*b/a + a*xksum+ n*log(k)+(k-1)*logxsum +(-b)/a*expksum  )
    return(out)
}


minusLogLikGompMak4p <- function(para) {

    b <- para[2] ; a <- para[1] ; l <- para[3] ; k <- para[4]; n <- sum(dx) ;
    x <- ifelse(x==0,1,x) ;
    xsum <- sum(dx*x) ; expksum <- sum(dx*exp(a*x^k)) ;
    logexpplusl <- log(b*k*exp(a*x^k)*x^(k-1)+l) ; logexppluslsum <- sum(dx*logexpplusl) ;
    xk1expksum <- sum(dx*exp(a*x^k)*x^(k-1)) ;
    xksum <- sum(dx*x^k) ; expsum <- sum(dx*exp(a*x)) ;
    logxsum <-sum(dx*log(x)) ;
  
    out <- -(logexppluslsum -l*xsum+((-b)/a)*expksum + n*b/a )

    
    return(out) 
}

minusLogLikWeib.gr <- function(para) {
    b <- para[1] ; k <- para[2] ; n <- sum(dx) ;
    x <- ifelse(x==0,1,x) ; 
    xsum <- sum(dx*x) ; expsum <- sum(dx*exp(a*x)) ;
    xexpsum <- sum(dx*exp(b*x)*x) ;
    xksum <- sum(dx*x^k) ; logsum <- sum(dx*log(x)) ; 

    c(-(k/b - b^(k-1)*k*xksum),
      -(-b^k*x^k*logsum +logsum - b^k*log(b)*xksum + 1/k + log(b))
      )
   
 
}

minusLogLikWeib.hess <- function(para) {
    b <- para[1] ; k <- para[2] ; n <- sum(dx) ;
    
    x <- ifelse(x==0,1,x) ; 
    xsum <- sum(dx*x) ;  xksum <- sum(dx*x^k) ;
    logsum <- sum(dx*log(x)) ; logsum2 <- sum(dx*log(x)^2) ;  

    c( -((-b^(k-2)*k^2*xksum) + b^(k-2)*k*xksum - k/b^2),
       -((-b^(k-1)*k*xksum*logsum) - b^(k-1)*log(b)*k*xksum - b^(k-1)*xksum + 1/b),
       -((-b^(k-1)*k*xksum*logsum) - b^(k-1)*log(b)*k*xksum - b^(k-1)*xksum + 1/b),
       -((-b^k*xksum*logsum2) - 2*b^k*log(b)*xksum*logsum - b^k*log(b)^2*xksum - 1/k^2)
      )
  
}

minusLogLikGomp.gr <- function(para) {
    b <- para[2] ; a <- para[1] ; n <- sum(dx) ;
    x <- ifelse(x==0,1,x) ; 
    xsum <- sum(dx*x) ; expsum <- sum(dx*exp(a*x)) ;
    xexpsum <- sum(dx*exp(b*x)*x) ;
    c( -( xsum - n*b/a^2 + b/a^2*expsum -b/a*xexpsum),
       -( n*(a+b)/(a*b) - 1/a*expsum )   
      ) 
 
}

minusLogLikGomp.hess <- function(para) {
    b <- para[2] ; a <- para[1] ; n <- sum(dx) ;
    x <- ifelse(x==0,1,x) ; 
    xsum <- sum(dx*x) ; x2sum <- sum(dx*x*x) ;
    expsum <- sum(dx*exp(a*x)) ;    xexpsum <- sum(dx*exp(b*x)*x) ;
     
     a11 <- -((-(b*x2sum*expsum)/a) - (2*b*xexpsum)/a^2 - (2*b*expsum)/a^3 + (2*b)/a^3) ;
     a12 <- -((xexpsum)/a + expsum/a^2 - 1/a^2) ;
     a21 <- a12 ;      a22 <-  1/b^2 ;
     c(a11,a12,a21,a22,a11*a22-a12*a12) 
    
 
}




minusLogLikGomp0<- function(para) {
    b <- para[2] ; ksi <- para[1] ; n <- sum(dx) ; xmax <-max(x) ;
  #  if (b<0) b <- 1e-5 ;  if (eta<0) eta <- 1e-5 ;
    xsum <- sum(dx*x) ; expsum <- sum(dx*exp(b*(x-xmax))) ;
  #  print(paste('log(eta): ',log(eta),'log(b): ',log(b),'b*xsum:  ',b*xsum)) ;
    h1l <- log(-((-n*log(ksi) + n*log(b) + n/ksi + b*xsum)*exp(-b*xmax) - 1/ksi*expsum)) ;
  #  print(paste('h1l: ',h1l)) ;
   
   # xexpsum <- sum(dx*exp(b*x)*x) ;
    out <- exp(h1l+b*xmax)
  #  out <-  - ( n*log(eta) + n*log(b) + n*eta + b*xsum  - eta*expsum )
 # attr(out, 'gradient') <- c(-(n / eta + n - expsum) ,
 #    -( n/b + xsum - eta *  xexpsum) ) 
    return(out)
}




minusLogLikGomp1 <- function(para) {
    eta <- para[1] ; n <- sum(dx) ;  #  b <- para[2] ;
                                        #  if (b<0) b <- 1e-5 ;  if (eta<0) eta <- 1e-5 ;
    b <- 0.11 ;
    xsum <- sum(dx*x) ; expsum <- sum(dx*exp(b*x)) ;
    xexpsum <- sum(dx*exp(b*x)*x) ;
    out <-  - ( n*log(eta) + n*log(b) + n*eta + b*xsum  - eta*expsum )
    attr(out, 'gradient') <- c(-(n / eta + n - expsum)  ) 
    return(out)
}


minusLogLikGompMak.gr <- function(para) {

    b <- para[2] ; a <- para[1] ; l <- para[3] ;  n <- sum(dx) ;
   
    xsum <- sum(dx*x) ; expsum <- sum(dx*exp(a*x)) ; xexpsum <- sum(dx*exp(a*x)*x) ;

    
    etaDsum <- sum(dx*b*exp(b*x)/(b*eta*exp(b*x)+lambda)) ;
    bDsum <- sum(dx*eta*(1+b*x)*exp(b*x)/(b*eta*exp(b*x)+lambda)) ;
    lambdaDsum <- sum(dx*1/(b*eta*exp(b*x)+lambda)) ;
    
    logexpsum <- sum(dx*log(1+lambda/(b*eta*exp(b*x)))) ;
    c(-(etaDsum + n / eta + n - expsum) ,
                              -(bDsum - eta *  xexpsum),
      -(lambdaDsum -xsum ))
    
    c(b*xexpsum/(b*expsum+l)-b*xexpsum/a+b*exp(a*x)/a^2-b/a^2,
      exp(a*x)/(b*exp(a*x)+l)-exp(a*x)/a+1/a,
      1/(b*exp(a*x)+l)-x);



    
}


loghxWeib <- function(x,para) {
    a <- para[1] ; b <- para[2] ;
   log(k) + k*log(b) + (k-1)*log(x)  
   
}

loghxGomp <- function(x,para) {
    a <- para[1] ; b <- para[2] ;
    a*x + log(b) ;
}

loghxGomp3p <- function(x,para) {
    a <- para[1] ; b <- para[2] ;  k <- para[3] ;
    log(b) + a*x^k + log(k)+(k-1)*log(x)
    
}

loghxGomp3p <- function(x,para) {
    a <- para[1] ; b <- para[2] ;  k <- para[3] ;
    log(b) + a*x^k + log(k)+(k-1)*log(x)
    
}


loghxGompM <- function(x,para) {
    a <- para[1] ; b <- para[2] ;  l <- para[3] ;
    log(b*exp(a*x) + l) 
    
}

loghxGompM4p <- function(x,para) {
    a <- para[1] ; b <- para[2] ;  l <- para[3] ;  k <- para[4] ;
    log(b*k*exp(a*x^k)*x^(k-1)+l) 
}


checkAllModels <- function(dX) {
    x <- dX$x ; dx <- dX$dx ;
    startW <- c( 0.02, 8) ;  
    startG <- c(0.11, 3e-06) ;   startGM <- c( 0.1, 1e-06,0.0003) ;
    startG3 <- c(0.11,5e-6,1) ;   startGM4 <- c(0.05, 4e-06, 0.0001,1) ;
     startG5 <- c(0.11,5e-6,0.5,1,85) ;

    w2x<-optimx(startW,fn=minusLogLikWeib,method=c("nlminb"),lower=c(0.01,1))
    g2x<-optimx(startG,fn=minusLogLikGomp,method=c("nlminb"),lower=c(0-01,1e-6))
    gMx<-optimx(startGM,fn=minusLogLikGompMak,method=c("nlminb"),lower=c(0.06,1e-12,1e-12))
    g3x<-optimx(startG3,fn=minusLogLikGomp3p,method=c("nlminb"),lower=c(0.01,1e-12,0.9))
    gM4x<-optimx(startGM4,fn=minusLogLikGompMak4p,method=c("nlminb"),lower=c(0.01,1e-12,1e-7,0.8))
    mloglik <- c(w2x$value[1],g2x$value[1],gMx$value[1],g3x$value[1],gM4x$value[1]) ; nparm <- c(2,2,3,3,4) ;
    mres <- list(w2x=w2x[1,],g2x=g2x[1,],gMx=gMx[1,],g3x=g3x[1,],gM4x=gM4x[1,])
    list(loglik=loglik,nparm=nparm,mres=mres) 
}

#   startW <- c( 0.02, 8) ;    w2x<-optimx(startW,fn=minusLogLikWeib,method=c("nlminb"),lower=c(0.01,1))

compPlot3p <- function(x,dx,p2,p3,pm3,p2w=c(0.01127964,9.132166),savePng=0,ylim=c(0,5000)) {
    
    if (savePng==0) X11() ;
    if (savePng==1) png(filename='overhaul_gomp_f1.png',width=800,height=800) ;
    
    curve(100000*fGomp3p(x,p3),from = 0, to = 106,ylim=ylim,col=3,lwd=2,xlab="Age",ylab="Deaths")
    points(x,dx,type="l",col=1,lty=2,lwd=2)
    curve(100000*fGompMak(x,pm3),col=2,lwd=2,add=T)

    curve(100000*dgompertz(x,p2[1],p2[2]),col=4,lwd=2,add=T)
    curve(100000*dweibull(x,p2w[2],1/p2w[1]),col=6,lty=3,lwd=2,add=T)
    legend(2,4900,lty=c(2,1,1,1,3),col=c(1,4,2,3,6), lwd=c(2,2,2,2,2),legend=c("Deaths","Gompertz","Gompertz-Makeham","Gompertz exp(ax^k)","Weibull"))

    if (savePng==1) dev.off() ;
}



compPlot5p <- function(x,dx,p2,p5,p2w=c(0.01127964,9.132166),savePng=0,ylim=c(0,5000),clr=2,drawLeg=0,yr=2017) {
    
 
   # if (savePng==0) X11() ;
    if (savePng==1) png(filename='overhaul_gomp_f2.png',width=800,height=800) ;
    if (clr==2) mainT <- ""  else mainT <- paste("Year:",yr  ) ;
    curve(100000*fGompAugm5p(x,p5),from = 0, to = 106,ylim=ylim,col=clr,lwd=2,xlab="Age",ylab="Deaths", main=mainT)
    points(x,dx,type="l",col=1,lty=2,lwd=2)
    curve(100000*dgompertz(x,p2[1],p2[2]),lty=3,col=clr,lwd=2,add=T)
    if (drawLeg==1) legend(2,4900,lty=c(2,1,3),col=c(1,clr,clr), lwd=c(2,2,2),legend=c("Deaths","Gompertz overhauled","Gompertz basic")) ;

    if (savePng==1) dev.off() ;

}

compPlot6p <- function(x,dx,p2,p6) {
    
    X11() ;
    curve(100000*fGompMakAugm6p(x,p6),from = 0, to = 106,col=3)
    points(x,dx,type="l",col=1,lty=2)
    curve(100000*dgompertz(x,p2[1],p2[2]),col=4,add=T,lty=3)


}


overhaulPanel <- function(savePng=0) {
    graphics.off() ;
   
    if (savePng==0) X11(width=7,height=12) ;
    if (savePng==1) png(filename='overhaul_gomp_f3.png',width=800,height=1300) ;
    par(mfrow=c(4,2)) ;
    
    dx1967 <- prepareYearMFDx(sT,"1967") ;
 # M 1967
    x <- dx1967$dxM$x ;  dx <- dx1967$dxM$dx ;
    p2 <- c(0.0867, 0.000092); p5 <- c(0.088, 0.00008, 0.001, 0.1330523, 85.5) ; compPlot5p(x,dx,p2,p5,clr=4,drawLeg=1,yr=1967) ;
 # F 1967
    x <- dx1967$dxF$x ;  dx <- dx1967$dxF$dx ;
    p2 <- c(0.108, 1.245738e-05) ; p5 <- c(0.108, 1.266084e-05, 0.015, 0.5, 82) ; compPlot5p(x,dx,p2,p5,yr=1967) ;

    dx1997 <- prepareYearMFDx(sT,"1997") ;
 # M 1997
    x <- dx1997$dxM$x ;  dx <- dx1997$dxM$dx ;
    p2 <- c(0.09635933, 3.68891e-05) ; p5 <- c( 0.0985, 2.95e-05, 0.032, 0.5, 83) ; compPlot5p(x,dx,p2,p5,clr=4,drawLeg=0,yr=1997) ;
 # F 1997
    x <- dx1997$dxF$x ;  dx <- dx1997$dxF$dx ;
    p2 <- c( 0.1028357, 1.41335e-05) ;  p5 <- c(0.115, 4.825569e-06, 0.008, 0.35, 91   ) ; compPlot5p(x,dx,p2,p5,yr=1997) ;
                                        
    dx2007 <- prepareYearMFDx(sT,"2007") ;
 # M 2007
    x <- dx2007$dxM$x ;  dx <- dx2007$dxM$dx ;
    p2 <- c(0.1028357, 1.804533e-05); p5 <- c( 0.115, 6.914726e-06, 0.008, 0.35, 91) ; compPlot5p(x,dx,p2,p5,clr=4,drawLeg=0,yr=2007) ;
 # F 2007
    x <- dx2007$dxF$x ;  dx <- dx2007$dxF$dx ;
    p2 <- c( 0.1119258, 5.854791e-06) ;p5 <- c(0.1004395, 1.521759e-05, 0.01882413, 0.215882, 92.5  ) ; compPlot5p(x,dx,p2,p5,yr=2007) ;

    dx2017 <- prepareYearMFDx(sT,"2017") ;
 # M 2017
    x <- dx2017$dxM$x ;  dx <- dx2017$dxM$dx ;                                     
    p2 <- c(0.109486, 8.507181e-06) ; p5 <- c(0.1002597,1.847507e-05,0.01430287,0.2028746,92.5);  compPlot5p(x,dx,p2,p5,clr=4,drawLeg=0,yr=2017) ;
 # F 2017   
    x <- dx2017$dxF$x ;  dx <- dx2017$dxF$dx ; 
    p2 <- c(0.118,3.19e-06) ; p5 <- c(0.115,3.913839e-06,0.03849455,0.3461503,92.5) ; compPlot5p(x,dx,p2,p5,yr=2017) ;
  
    if (savePng==1)  dev.off() ;

}

  


movAve3 <- function(y,w=c(1,1,1)) {
    n <- length(y) ;
    y0 <- c(y[1],y,y[n]) ;
    yp <- c(y[1],y[1],y) ;
    yf <- c(y,y[n],y[n]) ;
    yma <- (yp*w[1]+y0*w[2]+yf*w[3])/sum(w)
    yma[2:(n+1)]
}



#p2 <- c(0.118,3.19e-06) ; p3 <- c(0.0118,6.392e-06,1.437) ; p3m <- c(0.125,1.665e-06,0.000234) ;


#compPlot(x,dx,p2,p5)

#w2x         p1       p2    value fevals gevals niter convcode kkt1  kkt2
#nlminb 0.01127964 9.132166 386771.3     25     34    12        0 TRUE FALSE
#g2x         p1           p2    value fevals gevals niter convcode kkt1 kkt2
#nlminb 0.1176342 3.189598e-06 373976.7     50     44    17        0   NA   NA
#gMx         p1           p2           p3    value fevals gevals niter convcode
#nlminb 0.125079 1.664866e-06 0.0002339014 372483.2    136    198    55        0
#g3x         p1           p2       p3    value fevals gevals niter convcode
#nlminb 0.01176981 6.392421e-06 1.436836 373191.8    200    389   118        1
#g5x         p1           p2         p3        p4   p5         value fevals
#nlminb   0.109981 5.994283e-06 0.03692791 0.3438071 92.5  373312    199
# g6x        p1           p2         p3        p4       p5           p6      value fevals gevals niter convcode 
#nlminb   0.1093708 6.048863e-06 0.03983736 0.3460905 91.57632 0.0004999988  373140  28     24     4        1   
# g6x        p1           p2         p3        p4       p5           p6    value fevals gevals niter convcode kkt1 kkt2 xtimes
#nlminb   0.1190457 2.757254e-06 0.03136957 0.4043552 92.50002 0.0001964004 372204.5    200    547    71        1   NA   NA  0.076
# g5x        p1           p2         p3        p4   p5         value fevals
#nlminb   0.108987 6.419392e-06 0.03849455 0.3461503 91.5  373235.1    134
# g5x        p1           p2         p3        p4   p5         value fevals
#nlminb   0.115 3.913839e-06 0.03849455 0.3461503 92.5  373461.5     29


#p2 <- c(0.118,3.19e-06) ; p3 <- c(0.0118,6.392e-06,1.437) ; p3m <- c(0.125,1.665e-06,0.000234) ;


                                     #p5 <- c(0.1093708, 6.215833e-06, 0.03983736, 0.3460905, 91.57632)
                                        # Very good fit... dxF2017
   


                                        # M2017
#    p5 <- c( 0.115, 3.913839e-06, 0.03849455, 0.3461503, 92.5   );
#    g5x<-optimx(p5,fn=minusLogLikGompAugm5p,method=c("Nelder-Mead","BFGS","nlm","L-BFGS-B","nlminb"),lower=c(0,0,0,0,0)) ;
#  


# dx2016 <- prepareYearMFDx(sT,"2016") ; x <- dx2016$dxF$x ;  dx <- dx2016$dxF$dx ;
# F2016
# p2 <- c( 0.1169112, 3.406272e-06); p5 <- c(0.1067405,7.906248e-06,0.02289794,0.2607802,92.5) ; compPlot5p(x,dx,p2,p5) ;
# M2016 
# x <- dx2016$dxM$x ;  dx <- dx2016$dxM$dx ;
# p2 <- c(0.107513,1.010393e-05);  p5 <- c(0.09925484,1.958454e-05,0.02686547,0.2784424,90.5); compPlot5p(x,dx,p2,p5) ;


# source("gompertz_fitting_1.R")
#
#> load("../data/sT-2018.RData")
#
#  g2x<-optimx(startG,fn=minusLogLikGomp,method=c("Nelder-Mead","BFGS","nlm","L-BFGS-B","nlminb"),lower=c(0,0))


#  startG <- c(0.11 3e-06)
#  startGM <- c( 0.1, 1.8e-10,0.0005)
#  startG3 <- c(0.11,0.00051,1)
#  startGM4 <- c(0.11,0.0005,0.005,1)
#


#  g2x<-optimx(startG,fn=minusLogLikGomp,method=c("Nelder-Mead","BFGS","nlm","L-BFGS-B","nlminb"),lower=c(0,0))
#  gMx<-optimx(startGM,fn=minusLogLikGompMak,method=c("Nelder-Mead","BFGS","nlm","L-BFGS-B","nlminb"),lower=c(0,0,0))
#  g3x<-optimx(startG3,fn=minusLogLikGomp3p,method=c("Nelder-Mead","BFGS","nlm","L-BFGS-B","nlminb"),lower=c(0,0,0))
#  g4Mx<-optimx(startGM4,fn=minusLogLikGompMak4p,method=c("Nelder-Mead","BFGS","nlm","L-BFGS-B","nlminb"),lower=c(0,0,0,0))
#  g5x<-optimx(startG5,fn=minusLogLikGompAugm5p,method=c("Nelder-Mead","BFGS","nlm","L-BFGS-B","nlminb"),lower=c(0,0,0,0,0))
#  g6x<-optimx(p6,fn=minusLogLikGompMakAugm6p,method=c("Nelder-Mead","BFGS","nlm","L-BFGS-B","nlminb"),lower=c(0,0,0,0,0,0))

# F2017
# 0.1176342 3.189598e-06
#p2 <- c(0.118,3.19e-06) ;




#minusLogLikGompAugm5p 

#F(x,a,b,lambda) <- 1 - exp(-lambda*x - a/b*(exp(b*x)-1))
#f(x,a,b,lambda) <- (a*exp(b*x)+lambda) * exp(-lambda*x - a/b*(exp(b*x)-1))
#lf(x,a,b,lambda) <- b*x + log(a) + log(1+lambda/a*exp(-b*x)) - lambda*x -a/b - a/b*exp(b*x)

#lfGompMak(a,b,l):=log(b*exp(a*x)+l)-l*x+((-b)/a)*(exp(a*x)-1);
#lfGomp3p(a,b,k):=log(b)+a*x^k+log(k)+(k-1)*log(x)+((-b)/a)*exp(a*x^k)+b/a;
#lfGompMak4p(a,b,k,l):=log(b*k*exp(a*x^k)*x^(k-1)+l)-l*x+((-b)/a)*exp(a*x^k)+b/a;

#gfGompMak:matrix([(b*x*exp(a*x))/(b*exp(a*x)+l)-(b*x*exp(a*x))/a+(b*exp(a*x))/a^2-b/a^2,exp(a*x)/(b*exp(a*x)+l)-exp(a*x)/a+1/a,1/(b*exp(a*x)+l)-x]);

#hfGompMak:matrix([(-(b^2*x^2*exp(2*a*x))/(b^2*exp(2*a*x)+2*b*l*exp(a*x)+l^2))+(b*x^2*exp(a*x))/(b*exp(a*x)+l)-(b*x^2*exp(a*x))/a+(2*b*x*exp(a*x))/a^2-(2*b*exp(a*x))/a^3+(2*b)/a^3,(-(b*x*exp(2*a*x))/(b^2*exp(2*a*x)+2*b*l*exp(a*x)+l^2))+(x*exp(a*x))/(b*exp(a*x)+l)-(x*exp(a*x))/a+exp(a*x)/a^2-1/a^2,-(b*x*exp(a*x))/(b^2*exp(2*a*x)+2*b*l*exp(a*x)+l^2)],[(-(b*x*exp(2*a*x))/(b^2*exp(2*a*x)+2*b*l*exp(a*x)+l^2))+(x*exp(a*x))/(b*exp(a*x)+l)-(x*exp(a*x))/a+exp(a*x)/a^2-1/a^2,-exp(2*a*x)/(b^2*exp(2*a*x)+2*b*l*exp(a*x)+l^2),-exp(a*x)/(b^2*exp(2*a*x)+2*b*l*exp(a*x)+l^2)],[-(b*x*exp(a*x))/(b^2*exp(2*a*x)+2*b*l*exp(a*x)+l^2),-exp(a*x)/(b^2*exp(2*a*x)+2*b*l*exp(a*x)+l^2),-1/(b^2*exp(2*a*x)+2*b*l*exp(a*x)+l^2)]);

#gfGomp3p:matrix([(-(b*x^k*exp(a*x^k))/a)+(b*exp(a*x^k))/a^2+x^k-b/a^2,(-exp(a*x^k)/a)+1/b+1/a,(-b*x^k*exp(a*x^k)*log(x))+a*x^k*log(x)+log(x)+1/k]);

#hfGomp3p:matrix([(-(b*x^(2*k)*exp(a*x^k))/a)+(2*b*x^k*exp(a*x^k))/a^2-(2*b*exp(a*x^k))/a^3+(2*b)/a^3,(-(x^k*exp(a*x^k))/a)+exp(a*x^k)/a^2-1/a^2,x^k*log(x)-b*x^(2*k)*exp(a*x^k)*log(x)],[(-(x^k*exp(a*x^k))/a)+exp(a*x^k)/a^2-1/a^2,-1/b^2,-x^k*exp(a*x^k)*log(x)],[x^k*log(x)-b*x^(2*k)*exp(a*x^k)*log(x),-x^k*exp(a*x^k)*log(x),(-a*b*x^(2*k)*exp(a*x^k)*log(x)^2)-b*x^k*exp(a*x^k)*log(x)^2+a*x^k*log(x)^2-1/k^2]);

#gfGompMak4p:matrix([(b*k*x^(2*k-1)*exp(a*x^k))/(b*k*x^(k-1)*exp(a*x^k)+l)-(b*x^k*exp(a*x^k))/a+(b*exp(a*x^k))/a^2-b/a^2,(k*x^(k-1)*exp(a*x^k))/(b*k*x^(k-1)*exp(a*x^k)+l)-exp(a*x^k)/a+1/a,(a*b*k*x^(2*k-1)*exp(a*x^k)*log(x))/(b*k*x^(k-1)*exp(a*x^k)+l)+(b*k*x^(k-1)*exp(a*x^k)*log(x))/(b*k*x^(k-1)*exp(a*x^k)+l)-b*x^k*exp(a*x^k)*log(x)+(b*x^(k-1)*exp(a*x^k))/(b*k*x^(k-1)*exp(a*x^k)+l),1/(b*k*x^(k-1)*exp(a*x^k)+l)-x]);

#hfGompMak4p:matrix([(-(b^2*k^2*x^(4*k-2)*exp(2*a*x^k))/(b^2*k^2*x^(2*k-2)*exp(2*a*x^k)+2*b*k*l*x^(k-1)*exp(a*x^k)+l^2))+(b*k*x^(3*k-1)*exp(a*x^k))/(b*k*x^(k-1)*exp(a*x^k)+l)-(b*x^(2*k)*exp(a*x^k))/a+(2*b*x^k*exp(a*x^k))/a^2-(2*b*exp(a*x^k))/a^3+(2*b)/a^3,(-(b*k^2*x^(3*k-2)*exp(2*a*x^k))/(b^2*k^2*x^(2*k-2)*exp(2*a*x^k)+2*b*k*l*x^(k-1)*exp(a*x^k)+l^2))+(k*x^(2*k-1)*exp(a*x^k))/(b*k*x^(k-1)*exp(a*x^k)+l)-(x^k*exp(a*x^k))/a+exp(a*x^k)/a^2-1/a^2,(-(a*b^2*k^2*x^(4*k-2)*exp(2*a*x^k)*log(x))/(b^2*k^2*x^(2*k-2)*exp(2*a*x^k)+2*b*k*l*x^(k-1)*exp(a*x^k)+l^2))-(b^2*k^2*x^(3*k-2)*exp(2*a*x^k)*log(x))/(b^2*k^2*x^(2*k-2)*exp(2*a*x^k)+2*b*k*l*x^(k-1)*exp(a*x^k)+l^2)+(a*b*k*x^(3*k-1)*exp(a*x^k)*log(x))/(b*k*x^(k-1)*exp(a*x^k)+l)+(2*b*k*x^(2*k-1)*exp(a*x^k)*log(x))/(b*k*x^(k-1)*exp(a*x^k)+l)-b*x^(2*k)*exp(a*x^k)*log(x)-(b^2*k*x^(3*k-2)*exp(2*a*x^k))/(b^2*k^2*x^(2*k-2)*exp(2*a*x^k)+2*b*k*l*x^(k-1)*exp(a*x^k)+l^2)+(b*x^(2*k-1)*exp(a*x^k))/(b*k*x^(k-1)*exp(a*x^k)+l),-(b*k*x^(2*k-1)*exp(a*x^k))/(b^2*k^2*x^(2*k-2)*exp(2*a*x^k)+2*b*k*l*x^(k-1)*exp(a*x^k)+l^2)],[(-(b*k^2*x^(3*k-2)*exp(2*a*x^k))/(b^2*k^2*x^(2*k-2)*exp(2*a*x^k)+2*b*k*l*x^(k-1)*exp(a*x^k)+l^2))+(k*x^(2*k-1)*exp(a*x^k))/(b*k*x^(k-1)*exp(a*x^k)+l)-(x^k*exp(a*x^k))/a+exp(a*x^k)/a^2-1/a^2,-(k^2*x^(2*k-2)*exp(2*a*x^k))/(b^2*k^2*x^(2*k-2)*exp(2*a*x^k)+2*b*k*l*x^(k-1)*exp(a*x^k)+l^2),(-(a*b*k^2*x^(3*k-2)*exp(2*a*x^k)*log(x))/(b^2*k^2*x^(2*k-2)*exp(2*a*x^k)+2*b*k*l*x^(k-1)*exp(a*x^k)+l^2))-(b*k^2*x^(2*k-2)*exp(2*a*x^k)*log(x))/(b^2*k^2*x^(2*k-2)*exp(2*a*x^k)+2*b*k*l*x^(k-1)*exp(a*x^k)+l^2)+(a*k*x^(2*k-1)*exp(a*x^k)*log(x))/(b*k*x^(k-1)*exp(a*x^k)+l)+(k*x^(k-1)*exp(a*x^k)*log(x))/(b*k*x^(k-1)*exp(a*x^k)+l)-x^k*exp(a*x^k)*log(x)-(b*k*x^(2*k-2)*exp(2*a*x^k))/(b^2*k^2*x^(2*k-2)*exp(2*a*x^k)+2*b*k*l*x^(k-1)*exp(a*x^k)+l^2)+(x^(k-1)*exp(a*x^k))/(b*k*x^(k-1)*exp(a*x^k)+l),-(k*x^(k-1)*exp(a*x^k))/(b^2*k^2*x^(2*k-2)*exp(2*a*x^k)+2*b*k*l*x^(k-1)*exp(a*x^k)+l^2)],[(-(a*b^2*k^2*x^(4*k-2)*exp(2*a*x^k)*log(x))/(b^2*k^2*x^(2*k-2)*exp(2*a*x^k)+2*b*k*l*x^(k-1)*exp(a*x^k)+l^2))-(b^2*k^2*x^(3*k-2)*exp(2*a*x^k)*log(x))/(b^2*k^2*x^(2*k-2)*exp(2*a*x^k)+2*b*k*l*x^(k-1)*exp(a*x^k)+l^2)+(a*b*k*x^(3*k-1)*exp(a*x^k)*log(x))/(b*k*x^(k-1)*exp(a*x^k)+l)+(2*b*k*x^(2*k-1)*exp(a*x^k)*log(x))/(b*k*x^(k-1)*exp(a*x^k)+l)-b*x^(2*k)*exp(a*x^k)*log(x)-(b^2*k*x^(3*k-2)*exp(2*a*x^k))/(b^2*k^2*x^(2*k-2)*exp(2*a*x^k)+2*b*k*l*x^(k-1)*exp(a*x^k)+l^2)+(b*x^(2*k-1)*exp(a*x^k))/(b*k*x^(k-1)*exp(a*x^k)+l),(-(a*b*k^2*x^(3*k-2)*exp(2*a*x^k)*log(x))/(b^2*k^2*x^(2*k-2)*exp(2*a*x^k)+2*b*k*l*x^(k-1)*exp(a*x^k)+l^2))-(b*k^2*x^(2*k-2)*exp(2*a*x^k)*log(x))/(b^2*k^2*x^(2*k-2)*exp(2*a*x^k)+2*b*k*l*x^(k-1)*exp(a*x^k)+l^2)+(a*k*x^(2*k-1)*exp(a*x^k)*log(x))/(b*k*x^(k-1)*exp(a*x^k)+l)+(k*x^(k-1)*exp(a*x^k)*log(x))/(b*k*x^(k-1)*exp(a*x^k)+l)-x^k*exp(a*x^k)*log(x)-(b*k*x^(2*k-2)*exp(2*a*x^k))/(b^2*k^2*x^(2*k-2)*exp(2*a*x^k)+2*b*k*l*x^(k-1)*exp(a*x^k)+l^2)+(x^(k-1)*exp(a*x^k))/(b*k*x^(k-1)*exp(a*x^k)+l),(-(a^2*b^2*k^2*x^(4*k-2)*exp(2*a*x^k)*log(x)^2)/(b^2*k^2*x^(2*k-2)*exp(2*a*x^k)+2*b*k*l*x^(k-1)*exp(a*x^k)+l^2))-(2*a*b^2*k^2*x^(3*k-2)*exp(2*a*x^k)*log(x)^2)/(b^2*k^2*x^(2*k-2)*exp(2*a*x^k)+2*b*k*l*x^(k-1)*exp(a*x^k)+l^2)-(b^2*k^2*x^(2*k-2)*exp(2*a*x^k)*log(x)^2)/(b^2*k^2*x^(2*k-2)*exp(2*a*x^k)+2*b*k*l*x^(k-1)*exp(a*x^k)+l^2)+(a^2*b*k*x^(3*k-1)*exp(a*x^k)*log(x)^2)/(b*k*x^(k-1)*exp(a*x^k)+l)+(3*a*b*k*x^(2*k-1)*exp(a*x^k)*log(x)^2)/(b*k*x^(k-1)*exp(a*x^k)+l)+(b*k*x^(k-1)*exp(a*x^k)*log(x)^2)/(b*k*x^(k-1)*exp(a*x^k)+l)-a*b*x^(2*k)*exp(a*x^k)*log(x)^2-b*x^k*exp(a*x^k)*log(x)^2-(2*a*b^2*k*x^(3*k-2)*exp(2*a*x^k)*log(x))/(b^2*k^2*x^(2*k-2)*exp(2*a*x^k)+2*b*k*l*x^(k-1)*exp(a*x^k)+l^2)-(2*b^2*k*x^(2*k-2)*exp(2*a*x^k)*log(x))/(b^2*k^2*x^(2*k-2)*exp(2*a*x^k)+2*b*k*l*x^(k-1)*exp(a*x^k)+l^2)+(2*a*b*x^(2*k-1)*exp(a*x^k)*log(x))/(b*k*x^(k-1)*exp(a*x^k)+l)+(2*b*x^(k-1)*exp(a*x^k)*log(x))/(b*k*x^(k-1)*exp(a*x^k)+l)-(b^2*x^(2*k-2)*exp(2*a*x^k))/(b^2*k^2*x^(2*k-2)*exp(2*a*x^k)+2*b*k*l*x^(k-1)*exp(a*x^k)+l^2),(-(a*b*k*x^(2*k-1)*exp(a*x^k)*log(x))/(b^2*k^2*x^(2*k-2)*exp(2*a*x^k)+2*b*k*l*x^(k-1)*exp(a*x^k)+l^2))-(b*k*x^(k-1)*exp(a*x^k)*log(x))/(b^2*k^2*x^(2*k-2)*exp(2*a*x^k)+2*b*k*l*x^(k-1)*exp(a*x^k)+l^2)-(b*x^(k-1)*exp(a*x^k))/(b^2*k^2*x^(2*k-2)*exp(2*a*x^k)+2*b*k*l*x^(k-1)*exp(a*x^k)+l^2)],[-(b*k*x^(2*k-1)*exp(a*x^k))/(b^2*k^2*x^(2*k-2)*exp(2*a*x^k)+2*b*k*l*x^(k-1)*exp(a*x^k)+l^2),-(k*x^(k-1)*exp(a*x^k))/(b^2*k^2*x^(2*k-2)*exp(2*a*x^k)+2*b*k*l*x^(k-1)*exp(a*x^k)+l^2),(-(a*b*k*x^(2*k-1)*exp(a*x^k)*log(x))/(b^2*k^2*x^(2*k-2)*exp(2*a*x^k)+2*b*k*l*x^(k-1)*exp(a*x^k)+l^2))-(b*k*x^(k-1)*exp(a*x^k)*log(x))/(b^2*k^2*x^(2*k-2)*exp(2*a*x^k)+2*b*k*l*x^(k-1)*exp(a*x^k)+l^2)-(b*x^(k-1)*exp(a*x^k))/(b^2*k^2*x^(2*k-2)*exp(2*a*x^k)+2*b*k*l*x^(k-1)*exp(a*x^k)+l^2),-1/(b^2*k^2*x^(2*k-2)*exp(2*a*x^k)+2*b*k*l*x^(k-1)*exp(a*x^k)+l^2)]);



#    plane1: 5*x + 4*y + 5*z =0 ;
#    draw3d(implicit(plane1,x,-4,4,y,-4,4,z,-6,6)) ;
#    draw3d(enhanced3d=true,implicit(plane1,x,-4,4,y,-4,4,z,-6,6)) ;

#    ellips1: x^2/3 + y^2 + z^2 = 3 ; 
#    draw3d(enhanced3d=true,implicit(ellips1,x,-4,4,y,-4,4,z,-6,6)) ;
    
