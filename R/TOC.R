TOC <- function(input, ref, mask=NULL, nthres, NAvalue=0) {

if(is.null(mask)) mask <- ref-ref+1
mask[mask == NAvalue] <- NA

input <- input*mask
ref <- ref*mask

refval <- getValues(ref)
inpval <- getValues(input)

srted.ref <- sort(refval, index.return=TRUE, na.last=NA)
srted.inp <- sort(inpval, index.return=TRUE, na.last=NA)

maxInp <- max(srted.inp$ix)
newThres <- ceiling(maxInp/nthres*seq(1:nthres))
newThres <- sort(newThres, decreasing=TRUE)
srted.inp$ordInp <- 1:maxInp

mrg.ir <- merge(srted.ref, srted.inp, by.x="ix", by.y="ix")
mrg.ir$thresB <- 0

ones.ref <- table(refval)[[2]]
zeros.ref <- table(refval)[[1]]

res <- cbind(newThres, "truePositive"=0, "truePosRate"=0, "falsePositive"=0, "falsePosRate"=0)

for (j in 2:(nrow(res))){
i <- newThres[j]
mrg.ir[which(mrg.ir$ordInp > i), "thresB"]  <- 1

crsstb <- table(mrg.ir$thresB, mrg.ir$x.x)
mrg.ir$thresB <- 0
res[j,"truePositive"] <- crsstb[2,2]
res[j,"truePosRate"] <- crsstb[2,2]/ones.ref*100
res[j,"falsePositive"] <- crsstb[2,1]
res[j,"falsePosRate"] <- crsstb[2,1]/zeros.ref*100

}

res <- as.data.frame(rbind(res, c(0, ones.ref, 100, zeros.ref, 100)))
return(res)
}
