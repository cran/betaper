\name{adonis.pertables}
\Rdversion{1.1}
\alias{adonis.pertables}
\alias{print.adonis.pertables}
\alias{plot.adonis.pertables}

\title{ 
Function to assess the efects of taxonomic uncertainty on permutational multivariate analysis of variance using distance matrices
}
\description{
This function asses the effects of taxonomic uncertainty on the R2 coefficients and the p-values of a permutational multivariate analysis of variance using distance matrices.
}
\usage{
adonis.pertables(formula = X ~ ., data, permutations = 5, method = "bray")
\method{plot}{adonis.pertables}(x, \dots)
}

\arguments{
  \item{formula}{A typical model formula such as 'Y ~ A + B*C', but where 'Y' is a pertables object (i.e. a list of simulated community data matrices obtained with \code{\link{pertables}}; 'A', 'B', and 'C' may be factors or continuous variables.}
  \item{data}{The data frame from which 'A', 'B', and 'C' would be drawn.}
  \item{permutations}{Number of replicate permutations used for the hypothesis tests (F tests) for each simulated community data matrices obtained with \code{pertables}.}
  \item{method}{The name of any method used in 'vegdist' to calculate pairwise distances.}
  \item{x}{\code{adonis.pertables} object to plot.}
  \item{\dots}{Additional graphical parameters passed to plot.}
}

\value{
\code{adonis.pertables} returns an object of class\code{adonis.pertables}, basically a list with the following components:
  \item{raw }{An object of class \code{adonis}, i.e. the results of applying \code{\link[vegan]{mantel}} to the original biological data table without the unidentified species. This includes p-values for each explanatory variable showing the probability of obtaining the same F statistic under different scenarios of taxonomic uncertainty.}
  \item{simulation }{A list with the results of the simulation: \code{F}, i.e. a data.frame with all the simulated pseudo-F (columns) for each explanatory variable (rows); \code{R2}, i.e. a data.frame with all the simulated R2 coefficients (columns) for each explanatory variable (rows); \code{pvalue}, i.e. a data.frame with all the simulated p-values (columns) for each explanatory variable (rows); \code{R2.quant}, i.e. a data.frame with the summary of \code{R2} by quantiles; \code{p.quant},, i.e. a data.frame with the summary of \code{pvalue} by quantiles.}

The objects of class \code{adonis.pertables} have  \code{print} and \code{plot} S3 methods for a simple access to results. See the examples. 
}
\references{
Cayuela, L., De la Cruz, M. and Ruokolainen, K. 
A method to incorporate the effect of taxonomic uncertainty on multivariate analyses of ecological data (\emph{Ecography, submitted})
}
\author{
Luis Cayuela and Marcelino de la Cruz
}

\seealso{\code{pertables}, \code{\link[vegan]{adonis}}}

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

Amazonia.adonis <- adonis.pertables(Amazonia100 ~ Ca + K + Mg + Na, data=soils)

Amazonia.adonis

plot(Amazonia.adonis)
}
}
\keyword{ multivariate }
