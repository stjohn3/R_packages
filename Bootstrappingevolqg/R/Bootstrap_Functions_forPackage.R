#' Bootstrapping the MeanMatrixStatistic from the evolqg package
#'
#' This function takes a single dataframe and calculates one of the 10 
#' statistics that the MeanMatrixStatistic function can produce. The user must 
#' provide the data frame and the specific statistics they would like to caculate
#' 
#'
#' @param d dataframe
#' @param i iteration number that will be provided by the boot package (you do not enter anything)
#' @param measurement one of the following statistics: "MeanSquaredCorrelation", "pc1.percent", "ICV","EigenSd","respondability","evolvability","conditional.evolvability", "autonomy", "flexibility", or "constraints"
#' @return boot data frame from which you can make further calculations
#' @export
boot.meanmatrixstats<- function(d, i, measurement){
  
  if(measurement=="MeanSquaredCorrelation"){
    d[i,2:19]%>%
      cov(.)%>%
      as.matrix(.)%>%
      MeanMatrixStatistics(.)%>%
      as.data.frame()%>%
      slice(1)%>%
      as.numeric()%>%
      return()  
  }
  else if(measurement=="pc1.percent"){
    d[i,2:19]%>%
      cov(.)%>%
      as.matrix(.)%>%
      MeanMatrixStatistics(.)%>%
      as.data.frame()%>%
      slice(2)%>%
      as.numeric()%>%
      return()
  }
  else if(measurement=="ICV"){
    d[i,2:19]%>%
      cov(.)%>%
      as.matrix(.)%>%
      MeanMatrixStatistics(.)%>%
      as.data.frame()%>%
      slice(3)%>%
      as.numeric()%>%
      return()
  }
  else if(measurement=="EigenSd"){
    d[i,2:19]%>%
      cov(.)%>%
      as.matrix(.)%>%
      MeanMatrixStatistics(.)%>%
      as.data.frame()%>%
      slice(4)%>%
      as.numeric()%>%
      return()
  }
  else if(measurement=="respondability"){
    d[i,2:19]%>%
      cov(.)%>%
      as.matrix(.)%>%
      MeanMatrixStatistics(.)%>%
      as.data.frame()%>%
      slice(5)%>%
      as.numeric()%>%
      return()
  }
  else if(measurement=="evolvability"){
    d[i,2:19]%>%
      cov(.)%>%
      as.matrix(.)%>%
      MeanMatrixStatistics(.)%>%
      as.data.frame()%>%
      slice(6)%>%
      as.numeric()%>%
      return()
  }
  else if(measurement=="conditional.evolvability"){
    d[i,2:19]%>%
      cov(.)%>%
      as.matrix(.)%>%
      MeanMatrixStatistics(.)%>%
      as.data.frame()%>%
      slice(7)%>%
      as.numeric()%>%
      return()
  }
  else if(measurement=="autonomy"){
    d[i,2:19]%>%
      cov(.)%>%
      as.matrix(.)%>%
      MeanMatrixStatistics(.)%>%
      as.data.frame()%>%
      slice(8)%>%
      as.numeric()%>%
      return()
  }
  else if(measurement=="flexibility"){d[i,2:19]%>%
      cov(.)%>%
      as.matrix(.)%>%
      MeanMatrixStatistics(.)%>%
      as.data.frame()%>%
      slice(9)%>%
      as.numeric()%>%
      return()}
  else if(measurement=="constraints"){
    d[i,2:19]%>%
      cov(.)%>%
      as.matrix(.)%>%
      MeanMatrixStatistics(.)%>%
      as.data.frame()%>%
      slice(10)%>%
      as.numeric()%>%
      return()
  }
  
}

#' Bootstrapping the RandomSkewer function from the evolqg package
#'
#' This function takes a single dataframe that is the combination of the two data frames of interest (example: rbind(pop1.df, pop2,df)).
#' It is important that the rbind'd single dataframe contains a column that identifies which rows belong to which population, so that the
#' function can split up the data based on the user entered population name. 
#' This function takes the single dataframe, splits it into two data frames (one per population), bootstraps covariance matrices for each,
#' and finally calculates stats from the RandomSkewers function
#' 

#'
#' @param d dataframe
#' @param i iteration number that will be provided by the boot package (you do not enter anything)
#' @param measurement one of the following statistics: "correlation", "probability", or "correlation_sd"
#' @param matrix.1.pop a character that gives the name of the first population (will be used for filtering)
#' @param matrix.2.pop a character that gives the name of the second population (will be used for filtering)
#' @return boot data frame from which you can make further calculations
#' @export
boot.randomskewer<-function(d, i, measurement, matrix.1.pop, matrix.2.pop){
  d[i,]%>%
    filter(dataframe==matrix.1.pop)%>%
    nrow()->length.matrix1
  #print(length.matrix1)
  
  d[i,]%>%
    filter(dataframe==matrix.2.pop)%>%
    nrow()->length.matrix2
  #print(length.matrix2)
  
  if(measurement=="correlation"){
    
    d[i,]%>%
      filter(dataframe==matrix.1.pop)%>%
      dplyr::select(2:19)%>%
      slice(1:length.matrix1)%>%
      cov(.)%>%
      as.matrix(.)->Matrix.1
    
    
    d[i,]%>%
      filter(dataframe==matrix.2.pop)%>%
      dplyr::select(2:19)%>%
      slice(1:length.matrix2)%>%
      cov(.)%>%
      as.matrix(.)->Matrix.2    
    
    RandomSkewers(Matrix.1, Matrix.2)%>%
      as.data.frame()%>%
      slice(1)%>%
      as.numeric()%>%
      return()  
  }
  else if(measurement=="probability"){
    
    d[i,]%>%
      filter(dataframe==matrix.1.pop)%>%
      dplyr::select(2:19)%>%
      slice(1:length.matrix1)%>%
      cov(.)%>%
      as.matrix(.)->Matrix.1
    
    d[i,]%>%
      filter(dataframe==matrix.2.pop)%>%
      dplyr::select(2:19)%>%
      slice(1:length.matrix2)%>%
      cov(.)%>%
      as.matrix(.)->Matrix.2    
    
    RandomSkewers(Matrix.1, Matrix.2)%>%
      as.data.frame()%>%
      slice(2)%>%
      as.numeric()%>%
      return()  
  }
  else if(measurement=="correlation_sd"){
    
    d[i,]%>%
      filter(dataframe==matrix.1.pop)%>%
      dplyr::select(2:19)%>%
      slice(1:length.matrix1)%>%
      cov(.)%>%
      as.matrix(.)->Matrix.1
    
    d[i,]%>%
      filter(dataframe==matrix.2.pop)%>%
      dplyr::select(2:19)%>%
      slice(1:length.matrix2)%>%
      cov(.)%>%
      as.matrix(.)->Matrix.2    
    
    RandomSkewers(Matrix.1, Matrix.2)%>%
      as.data.frame()%>%
      slice(3)%>%
      as.numeric()%>%
      return()  
  }
}

#' Bootstrapping the PCAsimilarity function from the evolqg package
#'
#' This function takes a single dataframe that is the combination of the two data frames of interest (example: rbind(pop1.df, pop2,df)).
#' It is important that the rbind'd single dataframe contains a column that identifies which rows belong to which population, so that the
#' function can split up the data based on the user entered population name. 
#' This function takes the single dataframe, splits it into two data frames (one per population), bootstraps covariance matrices for each,
#' and finally calculates the PCAsimilarity score
#' 

#'
#' @param d dataframe
#' @param i iteration number that will be provided by the boot package (you do not enter anything)
#' @param matrix.1.pop a character that gives the name of the first population (will be used for filtering)
#' @param matrix.2.pop a character that gives the name of the second population (will be used for filtering)
#' @return boot data frame from which you can make further calculations
#' @export
boot.pca.similarity<-function(d, i, matrix.1.pop, matrix.2.pop){
  
  #Get row lengths for each data frame
  d[i,]%>%
    filter(dataframe==matrix.1.pop)%>%
    nrow(.)->length.matrix1
  
  d[i,]%>%
    filter(dataframe==matrix.2.pop)%>%
    nrow(.)->length.matrix2
  
  #print("these are the number of rows in each matrix:")
  #print(length.matrix1)
  #print(length.matrix2)
  
  
  if(length.matrix1<2 | length.matrix2<2){
    return("NaN")
  } 
  
  else{
    #Filter for the first population, select columns 2-19, and then select rows 1:length of matrix
    # print("calculate covariance matrices for pop 1")
    d[i,]%>%
      filter(dataframe==matrix.1.pop)%>%
      dplyr::select(2:19)%>%
      dplyr::slice(1:length.matrix1)%>%
      cov(.)%>%
      as.matrix(.)->Matrix.1
    
    #print("calculate covariance matrices for pop 2")
    d[i,]%>%
      filter(dataframe==matrix.2.pop)%>%
      dplyr::select(2:19)%>%
      dplyr::slice(1:length.matrix2)%>%
      cov(.)%>%
      as.matrix(.)->Matrix.2    
    
    #print("running actual comparison")
    
    PCAsimilarity(Matrix.1, Matrix.2)%>%
      return()  
  }    
  
}

#' Bootstrapping the MatrixDistance function from the evolqg package
#'
#' This function takes a single dataframe that is the combination of the two data frames of interest (example: rbind(pop1.df, pop2,df)).
#' It is important that the rbind'd single dataframe contains a column that identifies which rows belong to which population, so that the
#' function can split up the data based on the user entered population name. 
#' This function takes the single dataframe, splits it into two data frames (one per population), bootstraps covariance matrices for each,
#' and finally calculates the PCAsimilarity score
#' 

#'
#' @param d dataframe
#' @param i iteration number that will be provided by the boot package (you do not enter anything)
#' @param matrix.1.pop a character that gives the name of the first population (will be used for filtering)
#' @param matrix.2.pop a character that gives the name of the second population (will be used for filtering)
#' @return boot data frame from which you can make further calculations
#' @export
boot.matrix.distance<-function(d, i, matrix.1.pop, matrix.2.pop){
  d[i,]%>%
    filter(dataframe==matrix.1.pop)%>%
    nrow()->length.matrix1
  #print(length.matrix1)
  
  d[i,]%>%
    filter(dataframe==matrix.2.pop)%>%
    nrow()->length.matrix2
  #print(length.matrix2)
  
  d[i,]%>%
    filter(dataframe==matrix.1.pop)%>%
    dplyr::select(2:19)%>%
    slice(1:length.matrix1)%>%
    cov(.)%>%
    as.matrix(.)->Matrix.1
  
  d[i,]%>%
    filter(dataframe==matrix.2.pop)%>%
    dplyr::select(2:19)%>%
    slice(1:length.matrix2)%>%
    cov(.)%>%
    as.matrix(.)->Matrix.2    
  
  MatrixDistance(Matrix.1, Matrix.2)%>%
    return()  
}

#' Prints the mean and confidence intervals for any of the above functions
#'
#' This function can be used to print the mean and upper & lower confidence intervals for any of the above functions in one go
#'
#' @param enter.dataframe dataframe to be bootstrapped. The type of dataframe is dependent upon the function you will use
#' @param iterations iteration number that will be provided to the boot package (User enters these and this values is used as i in other functions)
#' @param which.function name of the function you wish to bootstrap: boot.meanmatrixstats, boot.randomskewer, boot.pca.similarity, or boot.matrix.distance (no quotes need here)
#' @param stat DEFAULT=NULL, for function with only 1 calculation (PCAsimilarity or MatrixDistance) this parameter can be left blank. For functions that calculate many
#' different statistics (MeanMatrixStatistic or RandomSkewers), you should enter the name of the stat you'd like to calculate here in quotes (see other function help pages for more info).
#' @param matrix.1.pop DEFAULT=NULL, only enter population names (in quotes) for functions that require two populations for calculations. 
#' @param matrix.2.pop DEFAULT=NULL, only enter population names (in quotes) for functions that require two populations for calculations. 
#' @return boot data frame from which you can make further calculations
#' @export
meanandci<-function(enter.dataframe, iterations, which.function, stat=NULL, matrix.1.pop=NULL,matrix.2.pop=NULL){
  
  #uses this print for MeanMatrixStatistics
  if(is.null(matrix.1.pop)==TRUE){
    boot.mean<-boot(d=enter.dataframe,
                    statistic=which.function,
                    measurement=stat,
                    R=iterations)
    boot.ci<-boot.ci(boot.mean,conf = 0.95, type = "basic")
    
    print(paste("mean is:", boot.mean$t0, "confidence intervals are:", boot.ci$basic[,4],boot.ci$basic[,5]))  }
  else{
    #uses this pring for PCAsimilarity and MatrixDistance
    if(is.null(stat)==TRUE){
      boot.mean<-boot(d=enter.dataframe,
                      statistic=which.function,
                      R=iterations,
                      matrix.1.pop=matrix.1.pop,
                      matrix.2.pop=matrix.2.pop)
      boot.ci<-boot.ci(boot.mean,conf = 0.95, type = "basic")
      
      print(paste("mean is:", boot.mean$t0, "confidence intervals are:", boot.ci$basic[,4],boot.ci$basic[,5]))
    }
    else{
      #use this print for random skewers 
      boot.mean<-boot(d=enter.dataframe,
                      statistic=which.function,
                      measurement=stat,
                      R=iterations,
                      matrix.1.pop=matrix.1.pop,
                      matrix.2.pop=matrix.2.pop)
      boot.ci<-boot.ci(boot.mean,conf = 0.95, type = "basic")
      
      print(paste("mean is:", boot.mean$t0, "confidence intervals are:", boot.ci$basic[,4],boot.ci$basic[,5]))
    }
    
  }
}

#' Prints the mean and confidence intervals for any of the above functions
#'
#' This function can be used to print the mean and upper & lower confidence intervals for any of the above functions in one go, but it is specifcally
#' meant to be used with the small datasets for the PCAsimilarity comparison.
#'
#' @param enter.dataframe dataframe to be bootstrapped. The type of dataframe is dependent upon the function you will use
#' @param iterations iteration number that will be provided to the boot package (User enters these and this values is used as i in other functions)
#' @param which.function name of the function you wish to bootstrap: boot.meanmatrixstats, boot.randomskewer, boot.pca.similarity, or boot.matrix.distance (no quotes need here)
#' @param stat DEFAULT=NULL, for function with only 1 calculation (PCAsimilarity or MatrixDistance) this parameter can be left blank. For functions that calculate many
#' different statistics (MeanMatrixStatistic or RandomSkewers), you should enter the name of the stat you'd like to calculate here in quotes (see other function help pages for more info).
#' @param matrix.1.pop DEFAULT=NULL, only enter population names (in quotes) for functions that require two populations for calculations. 
#' @param matrix.2.pop DEFAULT=NULL, only enter population names (in quotes) for functions that require two populations for calculations. 
#' @return boot data frame from which you can make further calculations
#' @export
meanandci.werid.PCAdfs<-function(enter.dataframe, iterations, which.function, stat=NULL, matrix.1.pop=NULL,matrix.2.pop=NULL){
  #uses this part for PCAsimilarity and MatrixDistance
  boot.mean<-boot(d=enter.dataframe,
                  statistic=which.function,
                  R=iterations,
                  matrix.1.pop=matrix.1.pop,
                  matrix.2.pop=matrix.2.pop)
  
  boot.mean$t%>%
    as.data.frame()->dataset
  
  names(dataset)<-"boot.mean.estimates"
  
  dataset%<>%
    filter(boot.mean.estimates!="NaN")%>%
    mutate(boot.mean.estimates=as.numeric(boot.mean.estimates))
  
  n<-nrow(dataset)
  calculated.mean<-boot.mean$t0
  s<-sd(as.numeric(dataset$boot.mean.estimates))
  
  margin <- qt(0.975,df=n-1)*s/sqrt(n)
  
  lowerinterval <- calculated.mean - margin
  
  upperinterval <- calculated.mean + margin
  
  print(paste("the mean is:",boot.mean$t0))
  
  print(paste("confidence intervals are:", lowerinterval,upperinterval))
  
  
}
