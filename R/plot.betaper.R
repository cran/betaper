`plot.betaper` <-
function(x, xlab = "Distances (km)", 
     ylab="Sorensen's similarity index", pch = c(18,3), ...){

plot(x$perm[,1] ~ x$geodist, xlab = xlab, ylab = ylab, type = "n", 
...)
for (i in 1:length(x$perm[1,])) {
par(new=T)
abline(lm(x$perm[,i] ~ x$geodist), lty=1, col = "darkgrey")
}
abline(lm(apply(x$perm, 1, mean) ~ x$geodist))
}

