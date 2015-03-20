plot.TOC <- function(toc, labelThres=FALSE, digits=3, modelLeg="Model", ...){

population <- toc$population
prevalence <- toc$prevalence/population
units <- toc$units

tocd <- toc$TOCtable
if((!is.null(tocd$HitsP) & !is.null(tocd$"Hits+FalseAlarmsP"))==TRUE){
tocd$Hits <- tocd$HitsP
tocd$"Hits+FalseAlarms" <- tocd$"Hits+FalseAlarmsP"
}


#units <- strsplit(strsplit(CRSargs(crs(index)), "+units=")[[1]][2], " ")[[1]][1]
old.par <- par(no.readonly = TRUE)
par(oma = c(0, 0, 0, 4))
par(mgp = c(1.5, 1, 0))

plot(c(0, population*(1-prevalence), population), c(0, 0, prevalence * population), type="l", lty="dashed", 
     xlab=paste0("Hits+False Alarms (", units, ")"), ylab=paste0("Hits (", units, ")"), 
     lwd=2, col=rgb(128,100,162, maxColorValue=255), bty="n", xaxt="n", yaxt="n", xlim=c(0, 1.05*population), 
     ylim=c(0, 1.05*prevalence * population), asp=1/prevalence, ...)

axis(1, pos = 0, xaxp = c(0, population, 5))
axis(2, pos = 0, yaxp = c(0, prevalence * population, 5))

# maximum
lines(c(0, prevalence * population, population), c(0, prevalence * population, prevalence * population), 
      lty="dotdash", lwd=2, col=rgb(79,129,189, maxColorValue=255)) 

# hits+misses
lines(c(0, population), rep(prevalence*population, 2), lwd=3, col=rgb(146,208,80, maxColorValue=255))

# uniform
lines(c(0, population), c(0, prevalence*population), lty="dotted", lwd=2, col=rgb(0,0,255, maxColorValue=255))

#lines(tocd$"Hits+FalseAlarms", tocd$maximum, lty="dotdash", lwd=2, col=rgb(79,129,189, maxColorValue=255))

# model
lines(tocd$"Hits+FalseAlarms", tocd$Hits, lwd=2, col=rgb(255,0,0, maxColorValue=255))
points(tocd$"Hits+FalseAlarms", tocd$Hits, pch=17, col=rgb(255,0,0, maxColorValue=255))
if(labelThres == TRUE) text(tocd$"Hits+FalseAlarms", tocd$Hits, round(as.numeric(tocd$Threshold), digits))

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

legend("right", c("Hits+Misses", "Maximum", modelLeg, "Uniform", "Minimum"), 
       col = c(rgb(146,208,80, maxColorValue=255), rgb(79,129,189, maxColorValue=255), rgb(255,0,0, maxColorValue=255), rgb(0,0,255, maxColorValue=255), rgb(128,100,162, maxColorValue=255)), 
       lty = c(1, 4, 1, 3, 2), pch = c(NA, NA, 17, NA, NA),
       merge = TRUE, bty="n", lwd=c(3, 2, 2, 2, 2))
par(old.par)

}


#setMethod("plot", signature("TOC"), plot.TOC)

#}
