\name{mantel.pertables}
\Rdversion{1.1}
\alias{mantel.pertables}
\alias{plot.mantel.pertables}
\alias{print.mantel.pertables}

\title{
Function to assess the efects of taxonomic uncertainty on Mantel tests
}
\description{
This function asses the effects of taxonomic uncertainty on the coefficient of correlation and the p-values of a Mantel test.
}
\usage{
mantel.pertables(pertab, env, dist.method = "bray", binary = FALSE, cor.method = "pearson", permutations = 100)
\method{plot}{mantel.pertables} (x, xlab = "Environmental distance", ylab = "Sorensen's similarity index", pch = 19, \dots)
}

\arguments{
  \item{pertab}{
A pertables object (i.e. a list of simulated community data matrices obtained with \code{\link{pertables}}.
}
  \item{env}{
Data frame with the environmental variables.
}
  \item{dist.method}{
Method to compute the dissimilarity matrices from the biological and environmental data tables. One of the methods described in  function \code{\link[vegan]{vegdist}} of the package \pkg{vegan}.
}
  \item{binary}{
Value for the argument \code{binary} in the function \code{\link[vegan]{vegdist}} of the package \pkg{vegan}.
}
  \item{cor.method}{
Correlation method, as accepted by \code{\link[stats]{cor}}: "pearson", "spearman" or "kendall". 
}
  \item{permutations}{
Number of permutations in assessing significance.
}
\item{x}{
\code{mantel.pertables} object to plot.
}
  \item{xlab}{
Label to name x-axis
}
  \item{ylab}{
Label to name y-axis
}

\item{pch}{Plotting 'character', i.e., symbol to use in the distance decay plot. See  \code{\link[graphics]{points}} for examples of use of this graphical argument.}
 \item{\dots}{Additional graphical parameters passed to plot.}
}


\value{\code{mantel.pertables} returns an object of class\code{mantel.pertables}, basically a list with the following components:

\item{mantel }{A list with two components: \code{mantel.raw}, an object of class \code{'mantel'}, i.e. the results of applying \code{\link[vegan]{mantel}} to the original biological data table without the unidentified species, and \code{ptax}, a p-value showing the probability of obtaining the same mantel statistic under different scenarios of taxonomic uncertainty.}
\item{simulation }{A list with the results of the simulation: \code{results}, i.e. a data.frame with all the simulated mantel statistics and p-values; \code{mantel.quant}, i.e. a data.frame with the summary of \code{results} by quantiles; \code{vegdist}, i.e. a list with all the dissimilarity matrices employed.}

The objects of class \code{mantel.pertables} have \code{print} and \code{plot} S3 methods for a simple access to results. See the examples. 
}

\references{
Cayuela, L., De la Cruz, M. and Ruokolainen, K. 
A method to incorporate the effect of taxonomic uncertainty on multivariate analyses of ecological data (\emph{Ecography, submitted})
}
\author{
Luis Cayuela and Marcelino de la Cruz
}

\seealso{\code{pertables}, \code{\link[vegan]{mantel}}}

\examples{
  \dontrun{
require(vegan)

data(Amazonia)
data(soils)

# Define a new index that includes the terms used in the \code{Amazonia} dataset to define undetermined taxa at different taxonomic levels
index.Amazon <- c(paste("sp.", rep(1:20), sep=""), "Indet.", "indet.")

#Generate a pertables object (i.e. a list of biological data tables simulated from taxonomic uncertainty)
Amazonia100 <- pertables(Amazonia, index=index.Amazon, nsim=100)

# Assess the effects of taxonomic uncertainty on a Mantel test of biological dissimilarity  correlated to soil dissimilarity among sites:

Amazonia.mantel <- mantel.pertables(pertab=Amazonia100, env=soils, dist.method = "bray")

Amazonia.mantel

plot(Amazonia.mantel)
 
}
}
\keyword{ multivariate }

