TOC <- function(index, boolean, mask=NULL, nthres=NULL, thres=NULL, NAval=0, P=NA, Q=NA, progress=FALSE) {

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

# define the thresholds vector
ifelse(!is.null(thres), minInd <- min(thres), minInd <- min(indval, na.rm=TRUE))

ifelse(!is.null(nthres), newThres <- (maxInd - minInd)/(nthres-2)*(0:(nthres-2)) + minInd, ifelse(!is.null(thres), 
                                                                                                  newThres <- thres, newThres <- unique(indval)))
newThres <- sort(newThres, decreasing=TRUE)

# create results data.frame
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

# calculate population 
validPixels <- ones.bool + zeros.bool

population <- validPixels * res(index)[1] * res(index)[2]
if(!is.na(P) & !is.na(Q)){
population <- P + Q
}

# calculate TOC table as shown in TOCfigure1.xlsx created by Pontius
tocd$Model1 <- tocd$HitsRate/100
tocd$falseAlarms1 <- tocd$falseAlarmsRate/100
tocd$Uniform <- tocd$falseAlarms1
maxA <- tocd$A[nrow(tocd)]
maxB <- tocd$B[nrow(tocd)]
tocd$m <- (maxA - tocd$A)/(maxA + maxB) 
tocd$h <- tocd$A/(maxA + maxB)
tocd$f <- tocd$B/(maxA + maxB) 
tocd$c <- 1 - tocd$m - tocd$h -tocd$f
prevalence <- tocd$h[nrow(tocd)]

tocd$Hits <- tocd$Model1 * prevalence * population
tocd$hitsFalseAlarms <- tocd$Hits + tocd$falseAlarms1*(1-prevalence)*population
tocd$hitsMisses <- prevalence*population
tocd$maximum <- pmin(tocd$hitsMisses, tocd$hitsFalseAlarms)
tocd$minimum <- pmax(0, tocd$hitsFalseAlarms + tocd$hitsMisses - population)
tocd$Uniform1 <- tocd$hitsMisses * tocd$hitsFalseAlarms / population

tocd1 <- tocd
tocd1[nrow(tocd1), "Threshold"] <- paste("<= ", minInd)

tocd2 <- tocd1[, c("Threshold", "hitsFalseAlarms", "Hits")]

# adjustment to population if user provides P and Q
if(!is.na(P) & !is.na(Q)){
tocd1$hitsFalseAlarmsP <- P * tocd1$Model1 + Q * tocd1$falseAlarms1
tocd1$HitsP <- P * tocd1$Model1
tocd2 <- tocd1[, c("Threshold", "hitsFalseAlarms", "Hits", "hitsFalseAlarmsP", "HitsP")]
}

# extract map units for plotting purposes with plot.TOC
units <- paste("square", strsplit(strsplit(CRSargs(crs(index)), "+units=")[[1]][2], " ")[[1]][1])

# calculate totalAUC in squared map units and AUC as a proportion 
id <- order(tocd2$hitsFalseAlarms)
totalAUC <- sum(tocd2$Hits[-length(tocd2$Hits)] * diff(tocd2$hitsFalseAlarms)) + sum(diff(tocd2$hitsFalseAlarms[id])*diff(tocd2$Hits[id]))/2 - ((prevalence * population)^2)/2 
AUC <- totalAUC/(population * prevalence * population - (prevalence * population)^2)

colnames(tocd2)[2] <- "Hits+FalseAlarms"
if (any(colnames(tocd2) == "hitsFalseAlarmsP")) colnames(tocd2)[4] <- "Hits+FalseAlarmsP"

thist <- hist(unique(index), breaks=c(sort(tocd$Threshold)), plot=FALSE)
tocd$counts <- c(0, thist$counts[length(thist$counts):1], 0)

# calculate uncertainty
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
output <- new("Toc", table=tocd2, prevalence=prevalence*population, population=population, units=units, AUC=AUC,  
               maxAUC = AUC + uncertain/2, minAUC = AUC - uncertain/2)
return(output) 

}
