TOCplot <- function(truePositive, truePositRate, falsePositive, falsePositRate, population){
tocd <- as.data.frame(cbind(truePositive, falsePositive, truePositRate, falsePositRate))
names(tocd) <- c("A",  "B", "truePositive", "falsePositive")
tocd$Model1 <- tocd$truePositive/100
tocd$falsePositive1 <- tocd$falsePositive/100
tocd$Uniform <- tocd$falsePositive1
maxA <- tocd$A[nrow(tocd)]
maxB <- tocd$B[nrow(tocd)]
tocd$m <- (maxA - tocd$A)/(maxA + maxB) 
tocd$h <- tocd$A/(maxA + maxB)
tocd$f <- tocd$B/(maxA + maxB) 
tocd$c <- 1 - tocd$m - tocd$h -tocd$f
prevalence <- tocd$h[nrow(tocd)]

tocd$Model11 <- tocd$Model1 * prevalence * population
tocd$hitsFalseAlarms <- tocd$Model11 + tocd$falsePositive1*(1-prevalence)*population
tocd$hitsMisses <- prevalence*population
tocd$maximum <- pmin(tocd$hitsMisses, tocd$hitsFalseAlarms)
tocd$minimum <- pmax(0, tocd$hitsFalseAlarms + tocd$hitsMisses - population)
tocd$Uniform1 <- tocd$hitsMisses * tocd$hitsFalseAlarms / population

plot(tocd$hitsFalseAlarms, tocd$minimum, type="l", lty="dashed", xlab="Hits + False Alarms (square kilometers)", ylab="Hits (square kilometers)", lwd=2, col=rgb(128,100,162, maxColorValue=255))
lines(tocd$hitsFalseAlarms, tocd$hitsMisses, lwd=3, col=rgb(146,208,80, maxColorValue=255))
lines(tocd$hitsFalseAlarms, tocd$Uniform1, lty="dotted", lwd=2, col=rgb(0,0,255, maxColorValue=255))
lines(tocd$hitsFalseAlarms, tocd$maximum, lty="dotdash", lwd=2, col=rgb(79,129,189, maxColorValue=255))
lines(tocd$hitsFalseAlarms, tocd$Model11, lwd=2, col=rgb(255,0,0, maxColorValue=255))
points(tocd$hitsFalseAlarms, tocd$Model11, pch=17, col=rgb(255,0,0, maxColorValue=255))

legend("right", c("Hits + Misses", "Maximum", "Model 1", "Uniform" , "Minimum"), col = c(rgb(146,208,80, maxColorValue=255), rgb(79,129,189, maxColorValue=255), rgb(255,0,0, maxColorValue=255), rgb(0,0,255, maxColorValue=255), rgb(128,100,162, maxColorValue=255)), lty = c(1, 4, 1, 3, 2), pch = c(NA, NA, 17, NA, NA),
       merge = TRUE, bty="n", lwd=c(3, 2, 2, 2, 2))
return(tocd)
}
