.plotROC <- function(object, labelThres=FALSE, digits=3, modelLeg="Model", ...){

rocd <- object@table

old.par <- par(no.readonly = TRUE)
par(oma = c(0, 0, 0, 4))
par(mgp = c(1.5, 1, 0))

plot(rocd$falseAlarms1, rocd$Model1, type="l", lty=1, 
     xlab=paste("False Alarms/(False Alarms + Correct Rejections)"), ylab=paste("Hits/(Hits+Misses)"), lwd=2, 
     col=rgb(255,0,0, maxColorValue=255), bty="n", xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), asp=1, ...)
axis(1, pos=0, xaxp = c(0, 1, 5))
axis(2, pos=0, xaxp = c(0, 1, 5))
points(rocd$falseAlarms1, rocd$Model1, pch=17, col=rgb(255,0,0, maxColorValue=255))

# maximum
lines(c(0, 1, 1), c(1, 1, 0), lwd=1, col="black")

# uniform
lines(c(0, 1), c(0, 1), lty="dotted", lwd=2, col=rgb(0,0,255, maxColorValue=255))

if(labelThres == TRUE) text(rocd$falseAlarms1, rocd$Model1, round(as.numeric(rocd$Threshold), digits))

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

legend("right", c(modelLeg, "Uniform"), 
       col = c(rgb(255,0,0, maxColorValue=255), rgb(0,0,255, maxColorValue=255)), lty = c(1, 3), pch = c(17, NA),
       merge = TRUE, bty="n", lwd=c(2, 2))
par(old.par)

}

