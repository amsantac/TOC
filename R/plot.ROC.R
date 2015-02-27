plot.ROC <- function(roc, labelThres=FALSE, digits=3){

# population <- toc$population
# prevalence <- toc$prevalence/population
# units <- roc$units

rocd <- roc$ROCtable
# if((!is.null(tocd$HitsP) & !is.null(tocd$"Hits+FalseAlarmsP"))==TRUE){
# tocd$Hits <- tocd$HitsP
# tocd$"Hits+FalseAlarms" <- tocd$"Hits+FalseAlarmsP"
# }


#units <- strsplit(strsplit(CRSargs(crs(index)), "+units=")[[1]][2], " ")[[1]][1]
old.par <- par(no.readonly = TRUE)
par(oma = c(0, 0, 0, 4))
par(mgp = c(1.5, 1, 0))

#plot(c(0, population*(1-prevalence), population), c(0, 0, prevalence * population), type="l", lty="dashed", 
# xlab=paste("False Alarms/(False Alarms + Correct Rejections)"), ylab=paste("Hits/(Hits+Misses)"), lwd=2, 
# col=rgb(128,100,162, maxColorValue=255), bty="n", xaxt="n", yaxt="n", xlim=c(0, 1.05*population), 
# ylim=c(0, 1.05*prevalence * population), asp=1/prevalence)

plot(rocd$falseAlarms1, rocd$Model1, type="l", lty=1, 
     xlab=paste("False Alarms/(False Alarms + Correct Rejections)"), ylab=paste("Hits/(Hits+Misses)"), lwd=2, 
     col=rgb(255,0,0, maxColorValue=255), bty="n", xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), asp=1)
axis(1, pos=0)
axis(2, pos=0)
points(rocd$falseAlarms1, rocd$Model1, pch=17, col=rgb(255,0,0, maxColorValue=255))

# maximum
lines(c(0, 1, 1), c(1, 1, 0), lwd=1, col="black")
#lines(rocd$falseAlarms1, rocd$Model1, type="l", lty="dashed", col=rgb(255,0,0, maxColorValue=255))

# hits+misses
#lines(c(0, population), rep(prevalence*population, 2), lwd=3, col=rgb(146,208,80, maxColorValue=255))

# uniform
lines(c(0, 1), c(0, 1), lty="dotted", lwd=2, col=rgb(0,0,255, maxColorValue=255))

#lines(tocd$"Hits+FalseAlarms", tocd$maximum, lty="dotdash", lwd=2, col=rgb(79,129,189, maxColorValue=255))

# model
# lines(tocd$"Hits+FalseAlarms", tocd$Hits, lwd=2, col=rgb(255,0,0, maxColorValue=255))
# points(tocd$"Hits+FalseAlarms", tocd$Hits, pch=17, col=rgb(255,0,0, maxColorValue=255))
if(labelThres == TRUE) text(rocd$falseAlarms1, rocd$Model1, round(as.numeric(rocd$Threshold), digits))

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

legend("right", c("Model", "Uniform"), 
       col = c(rgb(255,0,0, maxColorValue=255), rgb(0,0,255, maxColorValue=255)), lty = c(1, 3), pch = c(17, NA),
       merge = TRUE, bty="n", lwd=c(2, 2))
par(old.par)

}


#setMethod("plot", signature("TOC"), plot.TOC)

#}
