options(encoding="UTF-8")

# source("life_expectation_examples.R")
library(httr)
# rjstat is used for converting SSB JSON -> Data frame
library(rjstat)
# jsonlite is used mainly for converting metadata 
library(jsonlite)
# Reshape is used for filtering/transforming/grouping 
library(reshape)

# JSON handling functions and utility 

source("ssb_mortality_table_testing.R")
source("ssb-json-functions.R")


graphics.off()


createAndSaveELT <- function() {
    eLTTable <- getExpectedLifeTimeTable()
    save(eLTTable,file="../data/expectedLT_1966-2017.RData")
}


mkplot0 <- function(createPng=0) {

    pointsAndRegr <- function(eL,yrs=years,ages=c(67,eAges),col=1,newPlot=0,main="All") {

       rCoef <- matrix(0,nrow=length(ages),ncol=2) ; 
       for (i in 1:length(ages)) {   
        tS <-unlist(ages[i]+eL[(ages[i]+1),3:54]) ;
        rL  <- lm(tS~yrs) ;  regrL <- rL$coef ;       
        if (ages[i]==67)  plot(yrs,tS,type="l",ylim=c(75,100),lty=1,lwd=2,xlab="Year",ylab="Life expect",main=main,col=col)  
        else  points(yrs,tS,type="l",lty=2,lwd=2,col=col) ;
        abline(regrL[1],regrL[2],lty=4,col=col )
        rCoef[i,] <- regrL ;
       }
      rCoef  
    }

    computeCrossings <- function(rCoef) {
        nCat <- length(rCoef[,2]) ;
        for (i in 1:(nCat-1)) {
            xi <- (rCoef[i,1]-rCoef[i+1,1])/(rCoef[i+1,2]-rCoef[i,2])
            yi <- rCoef[i,1] + rCoef[i,2]*xi
            print(paste("xi: ",xi,"  yi: ",yi)) ;
        }
    }
    
    pointsAndCurves <- function(rCoef,ages=c(67,eAges),col=1,newPlot=0,main="Life expectancy increases") {
        if (newPlot==1) {
            plot(ages,100*rCoef[,2],ylim=c(0,10),xlab="Age",ylab="E(Increase 100 yrs)",main=main,col=col) ;
            points(ages,100*rCoef[,2],type="l",lty=4,col=col) ;
        }
        else  {
            points(ages,100*rCoef[,2],col=col) ;
            points(ages,100*rCoef[,2],type="l",lty=4,col=col) ;  
        }
    }

    relativePensions <- function(eL,yrs=years,d0=eLA[68,42],col=1,newPlot=0,main="Pension level with constant savings"  )  {
        if (newPlot==1) {
            plot(yrs,d0/eL[68,3:54],type="l",ylim=c(0.8,1.4),lty=1,lwd=2,xlab="Year",ylab="Relative pension level",main=main,col=col)
            abline(1,0)
        }
        else  points(yrs,d0/eL[68,3:54],type="l",lty=1,col=col,lwd=2)
    }
    
    X11(width=4,height=12)  ;
    if (createPng>0) png(file="life_expect_1.png",width=480,height=1200) ;
        
    par(mfrow=c(3,1)) ;  years <- 1966:2017 ;  eAges <- c(75,80,85,90,95) ;
  # Plot life  expectations for all   
    eLA <- eLTTable[eLTTable[,1]==0,] ;  rCoefA <- pointsAndRegr(eLA) ;
  # For men   
    eLM <- eLTTable[eLTTable[,1]==1,] ;    rCoefM <- pointsAndRegr(eLM,col=4,main="Men") ;
  # For women    
    eLF <- eLTTable[eLTTable[,1]==2,] ;     rCoefF <- pointsAndRegr(eLF,col=2,main="Women") ;
 
    if (createPng>0)  dev.off() ;

    X11(width=10,height=6) ;
    if (createPng>0) png(file="life_expect_2.png",width=800,height=480) ;
   
    par(mfrow=c(1,2)) ; 
  # 100 yr increases life expectancy  
    pointsAndCurves(rCoefA,col=1,newPlot=1) ;  pointsAndCurves(rCoefM,col=4) ;  pointsAndCurves(rCoefF,col=2) ;
  # Relative pensions 1966-2017
    relativePensions(eLA,col=1,newPlot=1) ; relativePensions(eLM,col=4) ; relativePensions(eLF,col=2) ;
    
    computeCrossings(rCoefA) ;

    if (createPng>0)  dev.off() ;

}



