

tailUQ <- function(x)
{
	if( length(x) > 1 )
	{
		return(tail(x, 2)[1] == 'UQ')
	}
	else
	{
		return(FALSE)
	}
}
	
combinedNorm <- function(combinedDataFile, outputFile, numCols, keepUqTag = TRUE, isQuantile = TRUE,  sep = '	', scaleTo = 5000000)
{
  comData <- read.table(combinedDataFile, sep = sep, row.names = 1, header = TRUE)

  countData <- comData[, 1:numCols]
  mapInfo <- comData[ , numCols + 1]

  isUqTag <- NULL

  if(keepUqTag)
  {
  	print("removing the non-uniquely mapped tags...")
  
  	mapInfoItems <- strsplit(as.character(mapInfo), ' ')

	filter <- sapply(mapInfoItems, tailUQ)

	countData <- countData[filter, ]
	mapInfo <- mapInfo[filter]
  }
  
  normData <- scaleByEmpiricalBayes(countData, scaleTo = scaleTo)

  if(isQuantile)
  {
      print("robust quantile normalization by Bolstad")
	  library(preprocessCore)
      normData.quantile <- normalize.quantiles.robust(normData)
	  rownames(normData.quantile) <- rownames(normData)
	  colnames(normData.quantile) <- colnames(normData)
	  normData <- normData.quantile
  }

  out <- cbind(data.frame(normData), mapInfo)

  write.table( t(c('tag', colnames(comData))), outputFile, sep = '	', quote = FALSE, row.names = FALSE, col.names = FALSE)

  write.table( out, outputFile, sep = '	', quote = FALSE, row.names = TRUE, col.names = FALSE, append = TRUE)
}

simpleGoodTurningEstimation <- function( freqTable , zeroCompensation = TRUE)
{
  rList <- as.integer(rownames(freqTable))
  nr <- as.integer(freqTable)

  if(rList[1] != 1)
    {
      print(paste("The smallest tag frequency observed is: ", rList[1]))
    }

  lenR <- length(rList)

  NOrig <- sum(rList * nr)

  P0 <- nr[1] / NOrig

  print(c("P0: ", toString(P0)))
  
  rNearRight <- rList[2:lenR]
  rNearLeft  <- c(0, rList[1:(lenR - 2)])  # only (lenR - 1) data points to fit
  
  zr <-  2 * nr[1:(lenR - 1)] /(rNearRight - rNearLeft)

  fit <- glm( log(zr) ~ log( rList[1:(lenR - 1)] ) )
  
  intercept <- fit$coefficients[1]
  exponent  <- fit$coefficients[2]

  print( c( "fitting results (intercept, exponent): ", toString(intercept), toString(exponent)) )

  ACoeff <- exp(intercept)

  rIndex0<- which(rNearRight - rList[1:(lenR - 1)] > 1)[1]

  if( is.na(rIndex0) )
    rIndex0 <- lenR

  #rIndex0 is the index of r in the original rList, you have to start to use smoothed curve to calculate the rStar

  if( rIndex0 > 1 )
    {
      
      rSub  <- rList[1:(rIndex0 - 1)]
      nrSub <- nr[1:(rIndex0 - 1)]

      #rSub and nrSub are the sub region of r, where they have next continuous r value;

      lenNrSub <- length(nrSub)
      
      x <- (rSub + 1) * (rList[2:(lenNrSub + 1)] / nrSub )

      y <- (rSub + 1) * (rSub + 1)^exponent /  rSub ^exponent

      varCmp <- 1.96 * sqrt( (rSub + 1)^2 * nr[2:(lenNrSub + 1)] / nrSub^2 * (1 + nr[2:(lenNrSub + 1)] / nrSub)  )

      rIndex1 <- which(abs(x -y) < varCmp)[1]

      isRStarDone <- FALSE
      
      if( is.na(rIndex1) ) #all of the sub should use x (original data)
        {
          rIndexStar <- rIndex0
        }
      else if( rIndex1 == 1) # all of list should use y (smoothed curve)
        {
          rStar <- (rList + 1) * (rList + 1)^exponent /  rList ^exponent
          isRStarDone = TRUE
        }
      else # 1 : (rIndex1 - 1) using x; rIndex1 : end using y;
        {
          rIndexStar <- rIndex1
        }
    }
  else
    # everything using y (smoothing curve)
    {
      rStar <- (rList + 1) * (rList + 1)^exponent /  rList ^exponent
      isRStarDone = TRUE
    }

  if(!isRStarDone)
    {
      rSub <- rList[1:(rIndexStar - 1)]
      nrSub <- nr[1:(rIndexStar - 1)]

      lenNrSub <- length(nrSub)
      
      rStar <- (rSub + 1) * (nr[2:(lenNrSub + 1)] / nrSub)

      rSub <- rList[rIndexStar:lenR]
      
      rStar <- c(rStar, (rSub + 1) * (rSub + 1)^exponent /  rSub ^exponent)
    }

  stopifnot(length(rStar) == lenR)
  
  # now renormalize r* to find the pStar
  
  NStar <- sum(rStar * nr)
  
  pStar <- (1 - P0) * rStar / NStar

#  numUnseen <- NStar * P0

#  print("p(0): ")
#  print(P0)
#  print(numUnseen)
  
  if( zeroCompensation )
    {
      pStar <- as.matrix(c(pStar[1], pStar))
    }
  else
    {
      pStar <- as.matrix(c(0, pStar))
    }

  rownames(pStar) <- as.character(c(0, rList))
      
#  print(pStar[1 , ])

  return(pStar)
}

scaleByEmpiricalBayes <- function(rawData, zeroCompensation = FALSE, scaleTo = 5000000)
{
  rowNum <- dim(rawData)[1]
  colNum <- dim(rawData)[2]

  pEst <- NULL

  T <- colSums(rawData)

  Tmax <- max(T)
  
  for( col in 1:colNum )
    {
      x <- rawData[, col]
      pStar <- simpleGoodTurningEstimation(table(x[x != 0]), zeroCompensation)
      
      pEst <- cbind(pEst, pStar[ as.character(x), ])
    }

  countEstPerM <- pEst * scaleTo

  countEstPerM <- as.matrix(countEstPerM)
  rownames(countEstPerM) <- rownames(rawData)
  
  return(countEstPerM)
}
combinedNorm("/Users/joshgruber/Documents/Data/GEO_datasets/SAGE-seq-GSE32017/GSM1115019-5262.txt", "combine-norm.dat", 1, keepUqTag = FALSE, isQuantile = FALSE, scaleTo=5000000)
