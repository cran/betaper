`betaper` <-
function(data, geodist, index = NULL, nsim = 100, vegdist.method = "bray", binary = FALSE, cor.method = "pearson"){
require(vegan)

# Index that determines the matching condition

  index0 <- c("Indet", "indet", "", " ", as.character(c(1:100)), 
             "sp", paste("sp", as.character(c(1:100)), sep=""),  
             paste("sp", as.character(c(1:100)), sep=" "))
  index <- c(index0, index)

# Create a species vector by pasting 'Genus' and 'Specific' columns

colnames(data)[1:3] <- c("Family", "Genus", "Specific")
data <- cbind(Species = paste(data$Genus, data$Specific), data)

# Conditioning sentences that establish whether the species have been 
# correctly determined at any taxonomic level (0 = correctly classified; 
# 1 = classified at genus resolution; 2 = classified at family resolution;
# 3 = family unknown)

cond <- rep(0, length(data$Species))
cond <- as.numeric(data$Specific%in%index)
cond <- cond + as.numeric(data$Genus%in%index)
cond <- cond + as.numeric(data$Family%in%index)


# Estimate the floristic distance between each pair of plots by 
# randomising the identification of unknown species

results <- geodist

for (p in 1:nsim){
print(paste("Esta es la simulaci<f3>n n<fa>mero", p))

data.A <- data

# STEP 1. Randomise non-identified species at any taxonomic level

for (i in 1:length(data.A[cond==3, 1])){

# Create an index to identify the plots where species 'i' is found

pr <- as.numeric(c(rep(0,4), data.A[cond==3,-c(1:4)][i,]))
pr <- ifelse(is.na(pr) == "TRUE", 0, pr)

# Select all species found in the same plots than species 'i' 

temp1 <- apply(as.data.frame(data.A[,pr>0]), 1, sum)

# Re-assign the non-identified species 'i' to any of the
# new possible identities

data.A$Species[cond==3][i] <-  
sample(rep(c(as.character(data.A$Species[cond==3][i]), as.character(data.A$Species[temp1==0])), 2), 1)
}

# STEP 2. Randomise non-identified species to family level

for (i in 1:length(data.A[cond==2, 1])){
pr2 <- as.numeric(c(rep(0,4), data.A[cond==2,-c(1:4)][i,]))
pr2 <- ifelse(is.na(pr2) == "TRUE", 0, pr2)

# Create a new data frame containing only the species 
# belonging to the same family than the species 'i'

data.f <- data.A[data.A$Family == data.A$Family[cond==2][i], ]

# Select all species found in the same plots than species 'i' 

temp1 <- apply(as.data.frame(data.f[,pr2>0]), 1, sum)

# Re-assign the non-identified species 'i' to any of the
# new possible identities

data.A$Species[cond==2][i] <- ifelse(dim(data.f)[1]==0, NA, 
sample(rep(c(as.character(data.A$Species[cond==2][i]), as.character(data.f$Species[temp1==0])), 2), 1))
}

# STEP 3. Randomise non-identified species to genus level

for (i in 1:length(data.A[cond==1, 1])){
pr3 <- as.numeric(c(rep(0,4), data.A[cond==1,-c(1:4)][i,]))
pr3 <- ifelse(is.na(pr3) == "TRUE", 0, pr3)

# Create a new data frame containing only the species 
# belonging to the same genus than the species 'i'

data.g <- data.A[data.A$Genus == data.A$Genus[cond==1][i], ]

# Select all species found in the same plots than species 'i' 

temp1 <- apply(as.data.frame(data.g[,pr3>0]), 1, sum)

# Re-assign the non-identified species 'i' to any of the
# new possible identities

data.A$Species[cond==1][i] <- ifelse(dim(data.g)[1]==0, NA, 
sample(rep(c(as.character(data.A$Species[cond==1][i]), as.character(data.g$Species[temp1==0])), 2), 1))
}


# As a result of re-assigning species identities there might be repeated 
# species names. We need to re-arrange the species matrix to merge species 
# with the same identities.

fun1 <- function(arg1) {apply(arg1, 2, sum)}

data.A2 <- do.call(rbind, by(data.A[,-c(1:4)], data.A$Species, fun1))
results <- cbind(results, vegdist(t(data.A2), method = vegdist.method, binary =  binary))

}

names(results) <- c("geodist", paste(rep("vegdis", nsim), c(1:nsim), 
         sep=""))

# 'cordis' is an object that contains Pearson's correlations between 
# geographical distances and each array of floristic distances resulting
# from simulation (as many correlation values as nsim)
 cordis <- apply(results[,-1], 2, function(x) 
    cor.test(x, 1 - results[, 1], method = cor.method)$estimate)

n <- table(cond)
names(n) <-ifelse(names(n)=="0", "Fully identified", names(n))
names(n) <- ifelse(names(n)=="1", "Identified to genus", names(n))
names(n) <- ifelse(names(n)=="2", "Identified to family", names(n))
names(n) <- ifelse(names(n)=="3", "Fully undetermined", names(n))
attributes(dimnames(n))$names <- "Taxonomic uncertainty"

final <- list(perm = 1 - results[, -1], cordis = cordis, 
    taxunc = n, geodist = geodist, cond = cond)
class(final) <- "betaper"
final

}

