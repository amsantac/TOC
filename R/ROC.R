ROC <- function(index, boolean, mask=NULL, nthres=NULL, thres=NULL, NAval=0, P=NA, Q=NA, progress=FALSE) {

if(!is.null(nthres) & !is.null(thres)) stop("Enter nthres OR thres as input, not both at the same time")

# create the mask if not provided
if(is.null(mask)) 
    mask <- boolean*0 + 1
mask[mask == NAval] <- NA

# mask out nodata cells in the index and boolean maps
index <- index*mask
boolean <- boolean*mask

# extract cell values from the boolean and index maps 
boolval <- getValues(boolean)
indval <- getValues(index)

# extract total number of cells with ones and zeros in the boolean map vector
boolvals <- boolval[!is.na(boolval)]
ones.bool <- sum(as.bit(boolvals))
zeros.bool <- length(boolvals) - ones.bool

# generic function for crosstabing two boolean vectors 
func_logical3   <- function(v1,v2){
    r1  <- sum(v1 & v2)
    r2  <- sum(v1 & !v2)
    return(c(r1, r2))
}

# extract only no NA values from the boolean and index vectors
# length of (not NA) boolean and index vectors must be equal as the mask was applied previously
boolval <- boolval[!is.na(boolval)]
indval <- indval[!is.na(indval)]

if (length(boolval)!=length(indval)) stop('different NA values in input maps')


zeroIndVal <- indval*0
maxInd <- max(indval, na.rm=TRUE)

ifelse(!is.null(thres), minInd <- min(thres), minInd <- min(indval, na.rm=TRUE))

ifelse(!is.null(nthres), newThres <- (maxInd - minInd)/(nthres-2)*(0:(nthres-2)) + minInd, ifelse(!is.null(thres), 
                                                                                                  newThres <- thres, newThres <- unique(indval)))
newThres <- sort(newThres, decreasing=TRUE)

res <- cbind(newThres, "Hits"=0, "HitsRate"=0, "falseAlarms"=0, "falseAlarmsRate"=0)

# loop for reclassifying the index vector for each new threshold
# and then perform the crosstab with the boolean vector

for (j in 2:(nrow(res))){
  i <- newThres[j]
  zeroIndVal[which(indval > i)]  <- 1
  
  xb <- as.bit(zeroIndVal)
  yb <- as.bit(boolval)
  crsstb <- func_logical3(xb,yb)
  
  res[j,"Hits"] <- crsstb[1]
  res[j,"HitsRate"] <- crsstb[1]/ones.bool*100
  res[j,"falseAlarms"] <- crsstb[2]
  res[j,"falseAlarmsRate"] <- crsstb[2]/zeros.bool*100
  zeroIndVal <- indval*0
  
  if(progress){
    pb <- txtProgressBar(min = 0, max = 100, style = 3)
    setTxtProgressBar(pb, round(j/nrow(res)*100, 0))
  }
}


tocd <- as.data.frame(rbind(res, c(NA, ones.bool, 100, zeros.bool, 100)))

names(tocd) <- c("Threshold", "A", "HitsRate", "B", "falseAlarmsRate")

validPixels <- ones.bool + zeros.bool

population <- validPixels * res(index)[1] * res(index)[2]
if(!is.na(P) & !is.na(Q)){
population <- P + Q
}

tocd$Model1 <- tocd$HitsRate/100
tocd$falseAlarms1 <- tocd$falseAlarmsRate/100
tocd$Uniform <- tocd$falseAlarms1

tocd1 <- tocd
tocd1[nrow(tocd1), "Threshold"] <- paste("<= ", minInd)

id <- order(tocd$falseAlarms1)
AUC <- sum(tocd$Model1[-length(tocd$Model1)] * diff(tocd$falseAlarms1)) + sum(diff(tocd$falseAlarms1[id])*diff(tocd$Model1[id]))/2 

# calculate uncertainty
thist <- hist(unique(index), breaks=c(sort(tocd$Threshold)), plot=FALSE)
tocd$counts <- c(0, thist$counts[length(thist$counts):1], 0)

uncertain <- 0
for (i in 2:(nrow(tocd)-2)){
if(tocd[i,"counts"] > 1) {
area <- (tocd[i,"falseAlarms1"] - tocd[i-1,"falseAlarms1"])*(tocd[i,"Model1"] - tocd[i-1,"Model1"])
uncertain <- uncertain + area
}
}

i <- i+1
if(tocd[i,"counts"] > 2) {
area <- (tocd[i,"falseAlarms1"] - tocd[i-1,"falseAlarms1"])*(tocd[i,"Model1"] - tocd[i-1,"Model1"])
uncertain <- uncertain + area
}

# output
output <- new("Roc", table=tocd1[, c(1,3,5:7)], AUC=AUC, maxAUC = AUC + uncertain/2, minAUC = AUC - uncertain/2)

return(output) 

}
