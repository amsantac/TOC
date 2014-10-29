TOC <- function(input, ref, mask=NULL, nthres, NAvalue=0) {

if(is.null(mask)) 
    mask <- ref - ref + 1
mask[mask == NAvalue] <- NA

input <- input*mask
ref <- ref*mask

refval <- getValues(ref)
inpval <- getValues(input)

srted.ref <- sort(refval, index.return=TRUE, na.last=NA)
srted.inp <- sort(inpval, index.return=TRUE, na.last=NA)

maxRank <- max(srted.inp$ix)
srted.inp$inputRank <- 1:maxRank

mrg.ir <- merge(srted.ref, srted.inp, by.x="ix", by.y="ix")
mrg.ir$thresB <- 0

newThres <- ceiling(maxRank/nthres*seq(1:nthres))
newThres <- sort(newThres, decreasing=TRUE)

ones.ref <- table(refval)[[2]]
zeros.ref <- table(refval)[[1]]

res <- cbind(newThres, "truePositive"=0, "truePosRate"=0, "falsePositive"=0, "falsePosRate"=0)

func_logical3   <- function(v1,v2){
    r1  <- sum(v1 & v2)
    r2  <- sum(v1 & !v2)
    return(c(r1, r2))
}

for (j in 2:(nrow(res))){
i <- newThres[j]
mrg.ir[which(mrg.ir$inputRank > i), "thresB"]  <- 1

xb <- as.bit(mrg.ir$thresB)
yb <- as.bit(mrg.ir$x.x)
crsstb <- func_logical3(xb,yb)

res[j,"truePositive"] <- crsstb[1]
res[j,"truePosRate"] <- crsstb[1]/ones.ref*100
res[j,"falsePositive"] <- crsstb[2]
res[j,"falsePosRate"] <- crsstb[2]/zeros.ref*100
mrg.ir$thresB <- 0
}

res <- as.data.frame(rbind(res, c(0, ones.ref, 100, zeros.ref, 100)))
return(res)
}
