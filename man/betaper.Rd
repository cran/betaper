\name{betaper}
\title{ Distance decay of similarity among floristic inventories in the face of taxonomic uncertainty }
\alias{betaper}
\alias{print.betaper}
\description{
  Function betaper computes the uncertainty range of the estimated correlation between floristic and geographical or environmental 
  distances in the face of taxonomic uncertainty.
}
\usage{
betaper(data, geodist, index = NULL, nsim = 100, vegdist.method = "bray", binary = FALSE, cor.method = "pearson")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ Community data matrix. The three first columns are factors referring to the family, genus and species specific names. 
  The remaining columns are numeric vectors indicating species abundances at each site. }
  \item{geodist}{ Dissimilarity matrix or a dist object of geographical or environmental distances. }
  \item{index}{ List of additional parameters to determine the level at which species have been identified. Default values include 'Indet', 
  'indet', 'sp', 'sp1' to 'sp100', 'sp 1' to 'sp 100', '', and ' '. }
  \item{nsim}{ Number of simulations of species' identities in assessing bounds to the estimated correlation. }
  \item{vegdist.method}{ Dissimilarity index, as in \code{\link[vegan]{vegdist}}. Partial match to "manhattan", "euclidean", "canberra", "bray", 
  "kulczynski", "jaccard", "gower", "morisita", "horn", "mountford", "raup" , "binomial" or "chao". }
  \item{binary}{ Perform presence/absence standardization before calculation of vegetation distances using \code{\link[vegan]{vegdist}}
  function. }
  \item{cor.method}{ Correlation method, as accepted by \code{\link[stats]{cor.test}}: "pearson", "spearman" or "kendall".  }
}
\details{
The method describes the manner to incorporate taxonomic uncertainty to the calculation of the correlation coefficient that measures 
the degree of correspondence between floristic and geographical or environmental distances, by estimating some credible bounds under 
plausible scenarios of re-assigned species identities. The approach involves the following steps:

Step 1. Codification of the species according to their degree of identification. Identified species are codified as 0, morphospecies identified 
only to genus (e.g. Eugenia sp1) as 1, only to family (e.g. Myrtaceae sp1) as 2, and fully unidentified morphospecies (e.g. Undetermined sp1) 
are codified as 3.

Step 2. Morphospecies identified only to genus are then randomly re-assigned in the group of species and morphospecies that share the same 
genus provided they are not found in the same sites. In the re-assignment of the species identity, the species can also get its own identity. 

Step 3. Once the re-assignment of morphospecies has been accomplished, we calculate a floristic distance matrix. A distance matrix is a 
square, symmetric matrix in which rows and columns are sites and cells are the floristic distance between pairs of sites. A single distance 
matrix thus represents the floristic pattern among sites for a single set of inventories. The Bray-Curtis' or Sorensen's coefficient is used by 
default to calculate floristic distance between sites. Note that although dissimilarity indexes are calculated with the vegdist function, function 
betaper calculates the complementary similarity index (i.e. 1 - dissimilarity index).

Step 4. The Pearson's correlation coefficient is computed between the floristic distance matrix and the geographical or environmental distance
matrix.

Step 5. Steps 2 to 4 are iterated n times. The process can be time-consuming if the floristic matrix is large. The mean and standard deviation, 
and the minimum and maximum values of the Pearson correlation coefficient, r, are then estimated. These measures provide credible bounds 
to the true correlation coefficient in the face of taxonomic uncertainty.
}
\value{
  The function return a list of class \code{betaper} with the following components
  \item{perm }{Permutated floristic distance matrices coerced to columns in a data frame. As many permutations (columns) as \code{nsim}.}
  \item{cordis }{Pearson's correlation values (\code{r}) between simulated floristic distances and geographical or environmental distances.}
  \item{taxunc }{Summary of the number of species fully identified (0), identified to genus (1), identified to family (2), or fully undetermined (3).}
  \item{geodist }{Dissimilarity matrix of geographical or environmental distances.}
  \item{cond }{Numeric vector indicating the assignment of species to any of the \code{taxunc} categories.}
  ...
}
\author{ Luis Cayuela \email{luis.cayuela@uah.es} and Marcelino de la Cruz \email{marcelino.delacruz@upm.es}}
\note{ The method uses floristic similarities instead of floristic dissimilarities, as calculated by \code{vegdist} function. The method simply 
computes floristic similarities by substracting the value output by \code{vegdist} to one
}
\seealso{ \code{\link[vegan]{mantel}} for computation of the Mantel statistic, \code{\link[vegan]{vegdist}} for methods to calculate vegetation
distances. See also \code{\link[vegan]{bioenv}} for methods to select environmental variables. 
}
\examples{
\dontrun{
require(vegan)
################################################
# Tree data from 16 0.1-ha plots from the Highlands of Chiapas, southern Mexico

data(HCP)
data(HCP.coords)

geodist.HCP <- dist(HCP.coords)/1000	# Geographical distances in kms

# Define a new index that includes the terms used in the HCP dataset to define undetermined taxa at genus and family levels
index.new <- c(paste("Genus", as.character(c(1:20)), sep=""), paste("Family", as.character(c(1:20)), sep=""))

beta.HCP <- betaper(HCP, geodist = geodist.HCP, index = index.new)
beta.HCP	# This results in a relatively low uncertainty range of the estimated correlation
beta.HCP$taxunc	# Number of identified and unidentified species at different taxonomic levels

################################################
# Tree data from nine lowland sites in Western Amazonia

data(Amazonia)
data(soils)

# Define a new index that includes the terms used in the Amazonia dataset to define undetermined taxa at different taxonomic levels
index.Amazon <- c(paste("sp.", rep(1:20), sep=""), "Indet.", "indet.")

beta.Amazonia <- betaper(Amazonia, geodist = dist(log(soils)), index = index.Amazon, nsim = 50)	# This can take a few minutes
beta.Amazonia	# The uncertainty range of the estimated correlation is much larger than in the previous example
beta.Amazonia$taxunc	# Number of identified and unidentified species at different taxonomic levels
}
}
\keyword{ multivariate }
\keyword{ htest }
