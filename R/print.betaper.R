`print.betaper` <-
function(x, ...) {
m <- c(min(x$cordis), mean(x$cordis)-2*sd(x$cordis), mean(x$cordis), mean(x$cordis)+2*sd(x$cordis), max(x$cordis))
names(m) <- c("minimum", "-2 SD", "mean", "+2 SD", "maximum")
names(n) <-ifelse(names(n)=="0", "Fully identified", names(n))
names(n) <- ifelse(names(n)=="1", "Identified to genus", names(n))
names(n) <- ifelse(names(n)=="2", "Identified to family", names(n))
names(n) <- ifelse(names(n)=="3", "Fully undetermined", names(n))
cat("Correlation values", "\n")
print(m)
cat("\n")
print(n)
}

