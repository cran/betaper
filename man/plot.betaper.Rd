\name{plot.betaper}
\alias{plot.betaper}
\title{ Plot of distance decay of similarity among floristic inventories in the face of taxonomic uncertainty }
\description{
  Plot of maximum and minimum correlation bounds between floristic and geographical or environmental distances under simulated scenarios 
  of species identities in the face of taxonomic uncertainty
}
\usage{
plot.betaper(x, xlab = "Distances (km)", ylab = "Sorensen's similarity index", pch = c(18, 3), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a \code{betaper} object. }
  \item{xlab}{ x-axis label. }
  \item{ylab}{ y-axis label. }
  \item{pch}{ A vector of the form c(x1, x2) specifying a symbol or a single character to represent the points that minimize and maximize 
  the correlation between floristic and geographical or environmental distances. }
 \item{...}{ any other graphical parameter. }

}
\details{
  This plot function produces a quick, standard plot with estimated minimum and maximum correlation bounds between floristic and geographical
  or environmental distances.
}
\author{ Luis Cayuela \email{luis.cayuela@uah.es} and Marcelino de la Cruz \email{marcelino.delacruz@upm.es} }
\seealso{ \code{\link{betaper}} }
\examples{
\dontrun{
require(vegan)
################################################
# Tree data from 16 0.1-ha plots from the Highlands of Chiapas, southern Mexico

data(HCP)
data(HCP.coords)

geodist.HCP <- dist(HCP.coords)/1000	# Geographical distances in kms

# Define a new index that includes the terms used in the HCPdataset to define undetermined taxa at genus and family levels
index.new <- c(paste("Genus", as.character(c(1:20)), sep=""), paste("Family", as.character(c(1:20)), sep=""))

beta.HCP <- betaper(HCP, geodist = geodist.HCP, index = index.new)
beta.HCP	# This results in a relatively low uncertainty range of the estimated correlation
beta.HCP$taxunc	# Number of identified and unidentified species at different taxonomic levels

plot(beta.HCP)

################################################
# Tree data from nine lowland sites in Western Amazonia

data(Amazonia)
data(soils)

# Define a new index that includes the terms used in the Amazonia dataset to define undetermined taxa at different taxonomic levels
index.Amazon <- c(paste("sp.", rep(1:20), sep=""), "Indet.", "indet.")

beta.Amazonia <- betaper(Amazonia, geodist = dist(log(soils)), index = index.Amazon, nsim = 50)	# This can take a few minutes
beta.Amazonia	# The uncertainty range of the estimated correlation is much larger than in the previous example
beta.Amazonia$taxunc	# Number of identified and unidentified species at different taxonomic levels

plot(beta.Amazonia)
  }
}
\keyword{ hplot }