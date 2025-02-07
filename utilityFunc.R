#utility functions required 

#read the markers from 
extract_markers <- function(ff, markers){
  ff@exprs[, markers]
}

#estimate L1-norm between channels
manhattan_dist <- function(a, b){
  dist <- abs(a-b)
  dist <- sum(dist)
  return(dist)
}

###############################################################################################################################
loadBaselineData <- function(myMetaData,marker_positions){
  
  #find the baseline sample
  a <- which(myMetaData$Baseline==1)
  msg<-paste('Processing baseline sample: ', myMetaData[a,'Filename'],' ID: ', myMetaData[a,'ID'] ,sep='')
  cat(msg)
  cat('\n')
  baseline <- myMetaData[a,'Filename']
  baseline <- as.character(baseline)
  fcs <- read.flowSet(baseline, transformation = FALSE, truncate_max_range = FALSE)
  
  panel_fcs <- pData(parameters(fcs[[1]]))
  panel_fcs <- as.data.frame(panel_fcs)
  panel_fcs$Position <- 1:nrow(panel_fcs)
  
  
  fcs <- fsApply(fcs, function(x, cofactor = 5){
    #careful here
    colnames(x) <- panel_fcs$name
    expr <- exprs(x)
    expr <- asinh(expr/cofactor)
    exprs(x) <- expr
    x
  })
  marker_names <- as.data.frame(panel_fcs$desc)
  myMarkers <- as.character(panel_fcs$desc[marker_positions$Positions])
  
  
  dataBaseline    <- flowFrame(fsApply(fcs, extract_markers, marker_positions$Positions))
  colnames(dataBaseline) <- myMarkers
  #make it easy to handle as data frame
  dataBaseline <- data.frame(dataBaseline@exprs)
  dataBaseline$ID <- 'Baseline'
  dataBaseline$TimePoint <- as.character(myMetaData[a,'TimePoint'])
  return(dataBaseline)
}

###############################################################################################################################
loadConditionData <- function(currentFile,currentID,TimePoint,marker_positions){
  
  #msg<-paste('Processing baseline sample: ', currentFile,' ID: ', currentID ,sep='')
  #cat(msg)
  #cat('\n')
  fcs <- read.flowSet(currentFile, transformation = FALSE, truncate_max_range = FALSE)
  
  panel_fcs <- pData(parameters(fcs[[1]]))
  panel_fcs <- as.data.frame(panel_fcs)
  panel_fcs$Position <- 1:nrow(panel_fcs)
  
  
  fcs <- fsApply(fcs, function(x, cofactor = 5){
    #careful here
    colnames(x) <- panel_fcs$name
    expr <- exprs(x)
    expr <- asinh(expr/cofactor)
    exprs(x) <- expr
    x
  })
  marker_names <- as.data.frame(panel_fcs$desc)
  myMarkers <- as.character(panel_fcs$desc[marker_positions$Positions])
  
  conditionData    <- flowFrame(fsApply(fcs, extract_markers, marker_positions$Positions))
  colnames(conditionData) <- myMarkers
  #make it easy to handle as data frame
  conditionData <- data.frame(conditionData@exprs)
  conditionData$ID <- currentID
  conditionData$TimePoint <- TimePoint
  return(conditionData)
}


################################################################################################

estimateSARA <- function(dataBaseline,conditionData,scoreMarkers,N){
  
  #here make the vectors of equal size
  #if(nrow(dataBaseline)>nrow(conditionData)){
  #  fixSize <- nrow(conditionData)
  #  a <- sample(nrow(dataBaseline))
  #  aSel <- a[1:fixSize]
  #  dataBaseline <- dataBaseline[aSel,]
  #}else{
    
  #  fixSize <- nrow(dataBaseline)
  #  a <- sample(nrow(conditionData))
  #  aSel <- a[1:fixSize]
  #  conditionData <- conditionData[aSel,]
  #}
  
  SARAscores <- matrix(-1,nrow = 1,ncol=ncol(dataBaseline)-2)
  IDX <- 1
  for(idx in scoreMarkers){
      
    str <- paste('  Loading channel:',idx,sep='')
    cat(str)
    cat('\n')
    channelA <- dataBaseline[,idx]
    channelB <- conditionData[,idx]
    
    #TODO: add density plots ?
    
    #generate the ecdf
    ecdfA<- ecdf(channelA)
    ecdfB<- ecdf(channelB)
    A <- ecdfA(channelA)
    B <- ecdfB(channelB)
    myEMD <- manhattan_dist(A,B)
    #perform permutations
    myCounter <- 0
    for(myPerm in 1:N){
      
      #shuffle the condition channel and generate random ecdf
      C <- sample(channelB)
      ecdfC <- ecdf(C)
      myC <- ecdfC(C)
      #permuted 
      permEMD <- manhattan_dist(A,myC)
      if(permEMD>myEMD){
        myCounter <- myCounter + 1
      }
    }
    pvalue <- myCounter/N
    #myEMD is a distance --> larger distance more different the samples
    #sign(mean(channelB)-mean(channelA) --> +1 if mean at the condition is greater than mean at baseline
    #1-pvalue is a score --> small p value means more significant changes so 1-p as score gets higher
    #overall higher SARA score means greater differences between the two experimental conditions
    myScore <- myEMD*sign(mean(channelB)-mean(channelA))*(1-pvalue)
    SARAscores[1,IDX] <- myScore
    IDX <- IDX +1
  
  }
  return(SARAscores)
}

#modify version of SARA
################################################################################################

estimateSARA_v2 <- function(dataBaseline,conditionData,scoreMarkers,N){
  
  #here make the vectors of equal size
  if(nrow(dataBaseline)>nrow(conditionData)){
    fixSize <- nrow(conditionData)
    #a <- sample(nrow(dataBaseline))
    #aSel <- a[1:fixSize]
    dataBaseline <- dataBaseline[fixSize,]
  }else{
    fixSize <- nrow(dataBaseline)
    #a <- sample(nrow(conditionData))
    #aSel <- a[1:fixSize]
    conditionData <- conditionData[fixSize,]
  }
  
  SARAscores <- matrix(-1,nrow = 1,ncol=ncol(dataBaseline)-2)
  IDX <- 1
  for(idx in scoreMarkers){
    
    str <- paste('  Loading channel:',idx,sep='')
    cat(str)
    cat('\n')
    channelA <- dataBaseline[,idx]
    channelB <- conditionData[,idx]
    
    #TODO: add density plots ?
    
    #generate the ecdf
    ecdfA<- ecdf(channelA)
    ecdfB<- ecdf(channelB)
    A <- ecdfA(channelA)
    B <- ecdfB(channelB)
    myEMD <- manhattan_dist(A,B)
    #perform permutations
    myCounter <- 0
    for(myPerm in 1:N){
      
      #permute the baseline
      #shuffle the condition channel and generate random ecdf
      D <- sample(channelA)
      ecdfD <- ecdf(D)
      myD <- ecdfD(D)
      
      #permute the condition channel and generate random ecdf
      C <- sample(channelB)
      ecdfC <- ecdf(C)
      myC <- ecdfC(C)
      #permuted 
      permEMD <- manhattan_dist(myD,myC)
      if(permEMD>myEMD){
        myCounter <- myCounter + 1
      }
    }
    pvalue <- myCounter/N
    #myEMD is a distance --> larger distance more different the samples
    #sign(mean(channelB)-mean(channelA) --> +1 if mean at the condition is greater than mean at baseline
    #1-pvalue is a score --> small p value means more significant changes so 1-p as score gets higher
    #overall higher SARA score means greater differences between the two experimental conditions
    myScore <- myEMD*sign(mean(channelB)-mean(channelA))*(1-pvalue)
    SARAscores[1,IDX] <- myScore
    IDX <- IDX +1
    
  }
  return(SARAscores)
}


estimateSARA_v3 <- function(dataBaseline,conditionData,scoreMarkers,N){
#here make the vectors of equal size
if(nrow(dataBaseline)>nrow(conditionData)){
  fixSize <- nrow(conditionData)
  #a <- sample(nrow(dataBaseline))
  #aSel <- a[1:fixSize]
  #dataBaseline <- dataBaseline[fixSize,]
}else{
  fixSize <- nrow(dataBaseline)
  #a <- sample(nrow(conditionData))
  #aSel <- a[1:fixSize]
  #conditionData <- conditionData[fixSize,]
}

SARAscores <- matrix(-1,nrow = 1,ncol=ncol(dataBaseline)-2)
IDX <- 1
for(idx in scoreMarkers){
  
  str <- paste('  Loading channel:',idx,sep='')
  cat(str)
  cat('\n')
  channelA <- dataBaseline[,idx]
  channelB <- conditionData[,idx]
  
  #TODO: add density plots ?
  
  #generate the ecdf
  ecdfA<- ecdf(channelA)
  ecdfB<- ecdf(channelB)
  A <- ecdfA(channelA)
  B <- ecdfB(channelB)
  A <- A[1:fixSize]
  B <- B[1:fixSize]
  myEMD <- manhattan_dist(A,B)
  #perform permutations
  myCounter <- 0
  for(myPerm in 1:N){
    
    #permute the baseline
    #shuffle the condition channel and generate random ecdf
    D <- sample(channelA)
    ecdfD <- ecdf(D)
    myD <- ecdfD(D)
    myD <- myD[1:fixSize]
    #permute the condition channel and generate random ecdf
    C <- sample(channelB)
    ecdfC <- ecdf(C)
    myC <- ecdfC(C)
    myC <- myC[1:fixSize]
    #permuted 
    permEMD <- manhattan_dist(myD,myC)
    if(permEMD>myEMD){
      myCounter <- myCounter + 1
    }
  }
  pvalue <- myCounter/N
  #myEMD is a distance --> larger distance more different the samples
  #sign(mean(channelB)-mean(channelA) --> +1 if mean at the condition is greater than mean at baseline
  #1-pvalue is a score --> small p value means more significant changes so 1-p as score gets higher
  #overall higher SARA score means greater differences between the two experimental conditions
  myScore <- myEMD*sign(mean(channelB)-mean(channelA))*(1-pvalue)
  SARAscores[1,IDX] <- myScore
  IDX <- IDX +1
  
}
return(SARAscores)
}

#estimate the value on the same "grid"
estimateSARA_v4 <- function(dataBaseline,conditionData,scoreMarkers,N){
  
  fixSize <- 10000
  
  SARAscores <- matrix(-1,nrow = 1,ncol=ncol(dataBaseline)-2)
  IDX <- 1
  for(idx in scoreMarkers){
    
    str <- paste('  Loading channel:',idx,sep='')
    cat(str)
    cat('\n')
    channelA <- dataBaseline[,idx]
    channelB <- conditionData[,idx]
    
    myGrid <- sample(c(channelA,channelB),size=fixSize)
    
    #TODO: add density plots ?
    
    #generate the ecdf
    ecdfA<- ecdf(channelA)
    ecdfB<- ecdf(channelB)
    A <- ecdfA(myGrid)
    B <- ecdfB(myGrid)

    myEMD <- manhattan_dist(A,B)
    #perform permutations
    myCounter <- 0
    for(myPerm in 1:N){
      
      #permute the baseline
      #shuffle the condition channel and generate random ecdf
      D <- sample(channelA)
      C <- sample(channelB)
      myGrid <- sample(c(channelA,channelB),size=fixSize)
      
      ecdfD <- ecdf(D)
      ecdfC <- ecdf(C)
      
      myD <- ecdfD(myGrid)
      myC <- ecdfC(myGrid)
      #permuted 
      permEMD <- manhattan_dist(myD,myC)
      if(permEMD>myEMD){
        myCounter <- myCounter + 1
      }
    }
    pvalue <- myCounter/N
    #myEMD is a distance --> larger distance more different the samples
    #sign(mean(channelB)-mean(channelA) --> +1 if mean at the condition is greater than mean at baseline
    #1-pvalue is a score --> small p value means more significant changes so 1-p as score gets higher
    #overall higher SARA score means greater differences between the two experimental conditions
    myScore <- myEMD*sign(mean(channelB)-mean(channelA))*(1-pvalue)
    SARAscores[1,IDX] <- myScore
    IDX <- IDX +1
    
  }
  return(SARAscores)
}




#extra function to see how the L1 norm is affected by sampling
studyPermutations <- function(dataB,dataC,myMarker,N){
  
  #version 1 - fix the small vector to the dimensions of the big
  scores_v1 <- matrix(-1,nrow = N,ncol = 1)
  for(reps in 1:N){
      dataBaseline <- dataB
      conditionData <- dataC
      #here make the vectors of equal size
      if(nrow(dataBaseline)>nrow(conditionData)){
        fixSize <- nrow(conditionData)
        a <- sample(nrow(dataBaseline))
        aSel <- a[1:fixSize]
        dataBaseline <- dataBaseline[aSel,]
      }else{
        fixSize <- nrow(dataBaseline)
        a <- sample(nrow(conditionData))
        aSel <- a[1:fixSize]
        conditionData <- conditionData[aSel,]
      }
      channelA <- dataBaseline[,myMarker]
      channelB <- conditionData[,myMarker]
      #generate the ecdf
      ecdfA<- ecdf(channelA)
      ecdfB<- ecdf(channelB)
      A <- ecdfA(channelA)
      B <- ecdfB(channelB)
      myEMD <- manhattan_dist(A,B)
      scores_v1[reps,1] <- myEMD
  }
  
  #version 2 - select the same number of cells
  scores_v2 <- matrix(-1,nrow = N,ncol = 1)
  fixSize <- 50000
  for(reps in 1:N){
    dataBaseline <- dataB
    conditionData <- dataC
    
    a <- sample(nrow(dataBaseline))
    aSel <- a[1:fixSize]
    dataBaseline <- dataBaseline[aSel,]
    a <- sample(nrow(conditionData))
    aSel <- a[1:fixSize]
    conditionData <- conditionData[aSel,]
    
    channelA <- dataBaseline[,myMarker]
    channelB <- conditionData[,myMarker]
    #generate the ecdf
    ecdfA<- ecdf(channelA)
    ecdfB<- ecdf(channelB)
    A <- ecdfA(channelA)
    B <- ecdfB(channelB)
    myEMD <- manhattan_dist(A,B)
    scores_v2[reps,1] <- myEMD
  }
  
  #version 3 - fix the small vector to the dimensions of the big
  scores_v3 <- matrix(-1,nrow = N,ncol = 1)
  for(reps in 1:N){
    dataBaseline <- dataB
    conditionData <- dataC
      
    a <- sample(nrow(dataBaseline))
    aSel <- a[1:fixSize]
    dataBaseline <- dataBaseline[aSel,]
      
    a <- sample(nrow(conditionData))
    aSel <- a[1:fixSize]
    conditionData <- conditionData[aSel,]
    
    channelA <- dataBaseline[,myMarker]
    channelB <- conditionData[,myMarker]
    #generate the ecdf
    ecdfA<- ecdf(channelA)
    ecdfB<- ecdf(channelB)
    A <- ecdfA(channelA)
    B <- ecdfB(channelB)
    myEMD <- manhattan_dist(A,B)
    scores_v3[reps,1] <- myEMD
  }
  
  channelA <- dataB[,myMarker]
  channelB <- dataC[,myMarker]
  #generate the ecdf
  ecdfA<- ecdf(channelA)
  ecdfB<- ecdf(channelB)
  A <- ecdfA(channelA)
  B <- ecdfB(channelB)
  myEMD <- manhattan_dist(A,B)
}





#############################
#original code from MATLAB SARA.m

# Statistical Analysis of Response Amplitude
# -----------------------------------------------------------------------
  # Score response distribution y from unperturbed distribution x,
#  using niter randomized values in null distribution
# 
# [score,pval,obs,mdiff,null] = SARA( x, y, niter )

# the (c) belong to the original publication, here is just a translation in R using the original code
#############################

emd_calc <- function( x, y, xedges ){
  # Earth mover distance
  # bin the data
  x1 <- histc( x, xedges )
  x1 <- x1$cnt
  y1 <- histc( y, xedges )
  y1 <- y1$cnt
  # CDF
  x1 <- cumsum( x1 / sum(x1) )
  y1 <- cumsum( y1 / sum(y1) )
  # emd
  emd <- sum( abs( x1 - y1 ) )
  return(emd)
  
}


check_significance <- function(obs,null,alpha){
  
  # Estimate confidence bounds on p-value using normal approximation to
  # binomial distribution
  # Only valid if #{null>=obs} > ~10
  # stop = true if 95% CI includes alpha
  # For details, see Knijnenberg et al., Bioinformatics 2009
  
  null <- na.omit(null)
  N <- nrow(null)
  a <- which(null>=obs)
  p_ecdf <- length(a)/N
  sigma <- p_ecdf*(1-p_ecdf)/N;
  lb <- p_ecdf - 1.96*sigma;
  # stop if bottom of CI is greater than alpha
  if( lb > alpha){
    stop <- TRUE
  }else{
    stop <- FALSE
  }
  return(stop)
}

originalSARA <- function(x,y,niter){
  
  null <- data.frame()
  #pool the data together
  xy <- c(x,y)
  n <- length(xy)
  nn <- length(x)
  
  # use fixed equipartition of space - default 100 but we could make more 
  xedges <- linspace(min(xy),max(xy),1000)
  
  # observed statistic
  obs <- emd_calc( x, y, xedges )
  
  # build null
  stop <- FALSE
  M <- 0
  ii <- 1
  while( stop==FALSE && ii < niter+1){
    
    a <- sample(xy)
    x1 <- a[1:nn]
    y1 <- a[(nn+1):length(xy)]
    #x1 <- sample(x)
    #y1 <- sample(y)
    myObs <- emd_calc(x1,y1,xedges)
    null <- rbind(null,myObs)
    if(myObs>=obs){
      M <- M + 1
    }
    
    # every 50 permutations, check significance
    # check requires #{null>obs} > ~10
    if( (ii%%50==0) && M >= 10){
      stop <- check_significance( obs, null, 0.05 )
      if(stop==TRUE){
        str <- paste('Stopping, the trick is working: ',ii,sep='')
        print(str)
       
      }
    }
    ii <- ii + 1
    
  }
  # # perms executed
  N <- ii-1;
  
  # compute p-value
  if(M>0){
    pseudocount <- 0
  }else{
    pseudocount <- 1;
  }
  a <- which(null>=obs)
  pval <- (pseudocount + length(a)) / N;
  mdiff <- median(y) - median(x)
  score <- obs * sign(mdiff) * (1 - pval )
  return(score)
}



