options(encoding="UTF-8")

# source("ssb_mortality_table_testing.R")
library(httr)
# rjstat is used for converting SSB JSON -> Data frame
library(rjstat)
# jsonlite is used mainly for converting metadata 
library(jsonlite)
# Reshape is used for filtering/transforming/grouping 
library(reshape)

source("ssb-json-functions.R")
#graphics.off()




# This query will fetch the whole table 07902

getQueryData07902 <- function() {
'{
  "query": [
  {
    "code":"Kjonn",
     "selection": {
    "filter": "all",
    "values":["*"] }
  },{
    "code":"AlderX",
     "selection": {
    "filter": "all",
    "values":["*"] }
  },{
    "code":"ContentsCode",
     "selection": {
    "filter": "all",
    "values":["*"] }
  },{
    "code":"Tid",
    "selection": {
    "filter": "all",
    "values":["*"] }
   }   ],
  "response": {
  "format": "json-stat"
  } 
}'
}




getMDValuesLabels07902 <- function(){

    getValuesAndLabels("07902")
    
}


getAllMortalityData07902 <- function(){

    getJSONData("07902",getQueryData07902())
    
}


#
#
#

testMetaDataAndData07902 <- function(){
    metaData <- getMDValuesLabels07902() ; 
    tableData <- getAllMortalityData07902()
    list(metaData=metaData,tableData=tableData)
}

#
# Example of "manual" setup for data extraction 
#

pickMortalityYears <- function(mdVL,yearsPicked=c(47,48,49,50,51)){
    mdVL$Kjonn[1,3] <- 10 ; # All
    mdVL$AlderX[1,3] <- 10 ; # All
    mdVL$ContentsCode[1,3] <- 10 ; # All
    for (i in 1:length(yearsPicked))  mdVL$Tid[yearsPicked[i],3] <- 1 ; # Selected
    
    mdVL
}

#
# Fetching expected life time data
#

getMortalityYearsData <- function(yearsPicked){
    
    metaData <- getMDValuesLabels07902() ;     
    mdVL <- pickMortalityYears(metaData,yearsPicked=yearsPicked) ;
    eQuery <- createQueryFromDF(mdVL) ;
    eData <- getJSONData("07902",eQuery)
  
    eData
}


#
#  Using melt&cast from reshape to transform to columns
#

transformToColumnsYearsData <- function(eData,yearsPicked){

    eData$expLT <- eData$value # To avoid default naming collision on "value"
    eData$value <- NULL ;  # Drop value column

    if (length(yearsPicked)==1) {
        eData$Tid <- NULL # Drop trivial column
        meltYears <- melt(eData,id=c("Kjonn","AlderX","ContentsCode")) ;
        meltYears$expLT <- NULL ;  # Drop trivial column
        colExp <- cast(meltYears,Kjonn+AlderX ~ContentsCode)
    }
    else {
        meltYears <- melt(eData,id=c("Kjonn","AlderX","ContentsCode","Tid")) ;
        meltYears$expLT <- NULL ;
        colExp <- cast(meltYears,Kjonn+AlderX ~Tid+ContentsCode)
    }
    colExp 
}


#
# Putting it all together - this function will give expected life times in columns
# Grouped by sex in first column  
#

getMortalityYearsDataTable <- function(yearsPicked=c(51)){

    eData <- getMortalityYearsData(yearsPicked)
    transformToColumnsYearsData(eData,yearsPicked)

}

#
# Example of "manual" setup for data extraction 
#

pickMortalityVars <- function(mdVL,contCode=3){
    mdVL$Kjonn[1,3] <- 10 ; # All
    mdVL$AlderX[1,3] <- 10 ; # All
    mdVL$ContentsCode[contCode,3] <- 1 ; # Expected life years
    mdVL$Tid[1,3] <- 10 ; # All
    
    mdVL
}


#
# Fetching expected life time data
#

getMortalityData <- function(contCode){
    
    metaData <- getMDValuesLabels07902() ;     
    mdVL <- pickMortalityVars(metaData,contCode=contCode) ;
    eQuery <- createQueryFromDF(mdVL) ;
    eData <- getJSONData("07902",eQuery)
    eData$ContentsCode <- NULL # Drop trivial column
    eData
}


#
#  Using melt&cast from reshape to transform to columns
#

transformToColumnsData <- function(eData){
    eData$expLT <- eData$value # To avoid default naming collision on "value"
    eData$value <- NULL ;  # Drop value column
    meltExpectation <- melt(eData,id=c("Kjonn","AlderX","Tid")) ;
    meltExpectation$expLT <- NULL ;  # Drop trivial column
    colExp <- cast(meltExpectation,Kjonn+AlderX~Tid)

    colExp 
}


#
# Putting it all together - this function will give expected life times in columns
# Grouped by sex in first column  
#

getMortalityDataTable <- function(contCode){

    eData <- getMortalityData(contCode)
    transformToColumnsData(eData)

}



getSurvivalTable <- function(){

    getMortalityDataTable(1)
  
}

getDeathsTable <- function(){

    getMortalityDataTable(2)
  
}



getExpectedLifeTimeTable <- function(){

    getMortalityDataTable(3)
  
}


getProbabilityOfDeathTable <- function(){

    getMortalityDataTable(4)
  
}


#survTable <- getSurvivalTable()

#numDeaths <- getDeathsTable()
#eLTTable <- getExpectedLifeTimeTable()
#probDeath <- getProbabilityOfDeathTable()

#mort2016 <- getMortalityYearsDataTable(yearsPicked=c(51))
#mort20142016 <- getMortalityYearsDataTable(yearsPicked=c(49,50,51))
#tableData <- getAllMortalityData07902()

#> t2017 <-  t07902[t07902$Tid==2017,]
#> table(t2017$ContentsCode)
#             Dode Dodssannsynlighet   ForvGjenLevetid   LevendePerTusen 
#              321               321               321               321 
#> lx <-  t2017[t2017$ContentsCode=="LevendePerTusen",5]
#> ex <-  t2017[t2017$ContentsCode=="ForvGjenLevetid",5]
#> qx <-  t2017[t2017$ContentsCode=="Dodssannsynlighet",5]
#> dx <-  t2017[t2017$ContentsCode=="Dode",5]
#> t2017d <-  t07902[t07902$Tid==2017&t07902$ContentsCode=="Dode",]
#> b2017 <- t2017d ;
#> b2017$value <- NULL ; b2017$ContentsCode <- NULL ;
#> dt2017 <- data.frame(b2017,lx,qx,dx,ex)

mortTableYear <- function(year,df=t07902) {
    tYr <-  df[df$Tid==year,] ; tYrD <-  tYr[tYr$ContentsCode=="Dode",] ;
    bYr <- tYrD ; bYr$value <- NULL ; bYr$ContentsCode <- NULL ; bYr$Tid <- NULL ;
    cCodes <- c("Dode","Dodssannsynlighet","ForvGjenLevetid","LevendePerTusen") 
    mData <- matrix(0,nrow=length(bYr[,1]),ncol=4)
    for (i in 1:4) mData[,i] <- tYr[tYr$ContentsCode==cCodes[i],5] 
    dfYr <- data.frame(bYr,mData)
    names(dfYr) = c("Kjonn","AlderX","dx","qx","ex","lx") ;
    dfYr 
}


#> dt2017All <- dt2017[dt2017$Kjonn==0,]
#> dt2017M <- dt2017[dt2017$Kjonn==1,]
#> dt2017F <- dt2017[dt2017$Kjonn==2,]
#> plot(as.numeric(dt2017All$AlderX),dt2017All$ex)
#> points(as.numeric(dt2017All$AlderX),dt2017M$ex,type="l",col=4)
#> points(as.numeric(dt2017All$AlderX),dt2017F$ex,type="l",col=2)

testPlotEX <- function(df,savePng=0) {
    plotter <- function(y,newP=0,col=1) {
        if (newP==1) plot(as.numeric(dfA$AlderX),y,type="l",xlab="Age",ylab="Expected years left",col=col)
        else points(as.numeric(dfA$AlderX),y,type="l",col=col)
        legend(65,75,col=c(1,2,4),lty=c(1,1,1),legend=c("All","Women","Men"))
    }

    X11()
    if (savePng>0) png(file="life_expect_0.png") ;
    dfA <- df[df$Kjonn==0,] ;  dfM <- df[df$Kjonn==1,] ; dfF <- df[df$Kjonn==2,] ;
    plotter(dfA$ex,newP=1) ;  plotter(dfM$ex,col=4) ;  plotter(dfF$ex,col=2) ;
    if (savePng>0) dev.off() ;
}
