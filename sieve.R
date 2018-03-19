# tay 20180307
# sieve.R - Eratosthenes again
#
# source("sieve.R")
#


sieveLoop <- function(MAXN=1000){
  numbers <- seq(1:MAXN) ; 
  maxSmallFactor <- round(sqrt(MAXN)) ;
  currPrime <- 2 ;   numbers[1] <- 0 ;
  while(currPrime <= maxSmallFactor) {
      checkNum <- currPrime*currPrime ;
      while(checkNum<=MAXN){
          numbers[checkNum] <- 0 ;
          checkNum <- checkNum + currPrime ;
      }
      currPrime <- currPrime + 1 ;
      while (numbers[currPrime]==0)  currPrime <- currPrime + 1 ;
  }

  numbers[numbers!=0]
}



sieveVector <- function(MAXN=1000){
    numbers <- 3:MAXN
    maxSmallFactor <- round(sqrt(MAXN)) ;
    currPrime <- 2 ; primes <- c(2) 
    while(currPrime <= maxSmallFactor) {
        numbers <- numbers[numbers %% currPrime!=0]
        currPrime <- numbers[1] ;  
        primes <- c(primes,currPrime)
    }
    primes <- c(primes,numbers)
    primes
}


sieveVectorPreAlloc <- function(MAXN=1000){
    numbers <- 3:MAXN
    maxSmallFactor <- round(sqrt(MAXN)) ;
    currPrime <- 2 ; primes <- numeric(round(MAXN/log(MAXN)+100)) ;
    primes[1] <- 2 ; nPrimes <- 1 ;
                                       
    while(currPrime <= maxSmallFactor) {
        numbers <- numbers[numbers %% currPrime!=0]
        currPrime <- numbers[1] ; nPrimes <- nPrimes + 1 ; 
        primes[nPrimes] <-  currPrime
    }
    primes <- primes[1:nPrimes] ;
    primes <- c(primes,numbers)
    primes
}


sieve <-  function(MAXN){
    sieveVector(MAXN)
}












