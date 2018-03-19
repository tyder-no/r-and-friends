#
# source("prime_testing.R")
#
#
#

source("sieve.R")



#
graphics.off()

#
#
#
#

testSieve <- function(MAXN) {

    print(t1 <- system.time(p1 <- sieveLoop(MAXN)))
    print(t2 <- system.time(p1 <- sieveVector(MAXN)))
    c(MAXN,t1[3],t2[3])
}



timeSieve <- function(MAXN=10000) {

  t1 <- system.time(p1 <- sieveVector(MAXN))
  list(time=t1,primes=p1)
}


comparePrimeSieves <- function(firstN=10000,nSteps=10,stepSize=2){

    currSize <- firstN ;
    resM <- matrix(0,nrow=nSteps,ncol=3,byrow=T)
    for (i in 1:nSteps) {
        t1 <- testSieve(currSize) ;
        resM[i,] <- c(currSize,t1[2],t1[3])    
        currSize <- currSize * stepSize ;
        if (currSize>500000) print(paste("CurrSize: ",currSize)) ;
    }

    resM
}


plotComparison <- function(lc0) {

    plotter <- function()  {
        plot(lc0[,1],lc0[,2],col=2,xlab="log(N)",ylab="log(Time,sec)",ylim=c(-7,6))
        points(lc0[,1],lc0[,3],col=3)
        abline(lmloop$coef[1],lmloop$coef[2],col=2,lty=3)
        abline(lmvect$coef[1],lmvect$coef[2],col=3,lty=3)
        legend(13,-4,lty=c(3,3),col=c(2,3),
               legend=c(paste("Loop, T= ",round(exp(lmloop$coef[1])*1000000,4),"*","N^",round(lmloop$coef[2],3)),
                        paste("Vector, T= ",round(exp(lmvect$coef[1])*1000000,4),"*","N^",round(lmvect$coef[2],3))))
    }   

    lmloop <- lm(lc0[,2]~lc0[,1])
    lmvect <- lm(lc0[,3]~lc0[,1])
     
    X11() ;  plotter() ;
    png(file="sieve_loop_vector.png") ;  plotter() ;  dev.off() ;
        
}


compareAndPlotComparisons <- function(nSteps=15) {
    c0 <-comparePrimeSieves(nSteps=15)
    lc0 <- log(c0)  
    plotComparison(lc0) 
}



timeAndCountPrimes <- function(firstN=100,nSteps=10,stepSize=2){

    currSize <- firstN ;
    resM <- matrix(0,nrow=nSteps,ncol=3,byrow=T)
    for (i in 1:nSteps) {
        t1 <- system.time(p1 <- sieveVector(currSize)) ;
        resM[i,] <- c(currSize,t1[3],length(p1))    
        currSize <- currSize * stepSize ;
        if (currSize>500000) print(paste("CurrSize: ",currSize)) ;
    }

    resM
}

primeProperties <- function(tc1) {

    plotterPNT <- function() {
        plot(log(tc1[,1]),tc1[,1]/tc1[,3],col=2,ylim=c(4,18),xlab="log(N)",ylab="N/pi(N)")
        points(log(tc1[,1]),log(tc1[,1])-1,type="l",lty=3,col=2)
        legend(12,8,lty=c(3),col=c(2),legend=c("Y=log(N)-1"))

    }

    
    
 # 1 Prime number theorem illustrated x/pi(x) = log(x)-1
    X11() ;  plotterPNT() ;
    png(file="sieve_PNT.png") ;  plotterPNT() ;  dev.off() ;
    
    

 # 2 Timing in R vector implementation
 #   lf <- lm(log(tc1[6:20,2]) ~ log(tc1[6:20,1]))
 # line y = -17.1 + 1.15x    
 # ie in this case approx time = 1/exp(17.1)*N^1.15
    
  #  X11()
  #  plot(log(tc1[,1]),log(tc1[,2]),col=2)  
  #  abline(coef(lf)[1],coef(lf)[2]) 
    
    
}


#> tc1
#          [,1]   [,2]    [,3]
# [1,]      100  0.000      26
# [2,]      200  0.000      47
# [3,]      400  0.001      79
# [4,]      800  0.000     140
# [5,]     1600  0.000     252
# [6,]     3200  0.001     453
# [7,]     6400  0.001     835
# [8,]    12800  0.002    1527
# [9,]    25600  0.004    2819
#[10,]    51200  0.007    5240
#[11,]   102400  0.017    9806
#[12,]   204800  0.034   18368
#[13,]   409600  0.080   34580
#[14,]   819200  0.231   65349
#[15,]  1638400  0.540  123842
#[16,]  3276800  1.109  235371
#[17,]  6553600  2.638  448183
#[18,] 13107200  6.628  855733
#[19,] 26214400 16.144 1637055
#[20,] 52428800 41.951 3137913


# tc1 <- timeAndCountPrimes(nSteps=20)    
# primeProperties(tc1)
