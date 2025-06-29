\name{rda_pertables}
\alias{rda_pertables}
\alias{plot.rda_pertables}
\alias{print.rda_pertables}

\title{
Function to assess the efects of taxonomic uncertainty on [Partial] Redundance Analysis
}
\description{
This function asses the effects of taxonomic uncertainty on two widely used parameters of a [Partial] Redundance Analysis, i.e. the 'percentage explained variance' (sometimes referred to as \emph{R-squared}) and the \emph{'pseudo-F' }. 
}
\usage{
rda_pertables(fml, data, scale=FALSE,...)
\method{plot}{rda_pertables} (x, pch = 18, ...)

}

\arguments{
  \item{fml}{
Model formula, where the left hand side gives a pertables object (i.e. a list of simulated community data matrices obtained with \code{\link{pertables}},  right hand side gives the constraining variables, and conditioning variables can be given within a special function \code{Condition}.
}
  \item{data}{
Data frame containing the variables on the right hand side of the model formula.
}

\item{scale}{
Scale species to unit variance (like correlations).
}

\item{x}{rda_pertables object to plot.
}
 \item{pch}{Plotting 'character', i.e., symbol to use in the RDA plot. See  \code{\link[graphics]{points}} for examples of use of this graphical argument.}
 \item{\dots}{Additional graphical parameters passed to plot.}

}
\details{
This function is a wrapper to submit a \code{\link{pertables}} object to \code{rda} function of the \pkg{vegan} package. See the documentation of \code{\link[vegan]{cca}} for details about \code{formula} and \code{Condition} use. } 

\value{\code{rda_pertables} returns  an object of class\code{'rda_pertables'}, basically a list with the following components:

\item{raw }{An object of class class\code{'rda'}. The results of applying \code{rda} to the original biological data table without the unidentified species.}
\item{simulation }{A list with the results of the simulation: \code{'results'}, i.e. a data.frame with all the simulated R-squared and pseudo-F values; \code{'rda.quant'}, i.e. a data.frame with the summary of \code{'results'} by quantiles; \code{'sites'} i.e. a list with the scores of  the sites of all the simulated data tables and  \code{'biplot'}, i.e.  a list with the scores of the environmental data in all the analyses}


The objects of class \code{'rda_pertables'} have  \code{print}  and \code{plot} S3 methods for a simple access to results. See the examples. 

}
\references{
Cayuela, L., De la Cruz, M. and Ruokolainen, K. (2011). A method to incorporate the effect of taxonomic uncertainty on multivariate analyses of ecological data. \emph{Ecography}, 34: 94-102. \doi{10.1111/j.1600-0587.2009.05899.x}.   
}
\author{
Luis Cayuela and Marcelino de la Cruz
}

\seealso{\code{pertables}, \code{\link[vegan]{cca}}}

\examples{
 


data(Amazonia)
data(soils)

# Define a new index that includes the terms used in the \code{Amazonia} dataset to define
# undetermined taxa at different taxonomic levels

index.Amazon <- c(paste("sp.", rep(1:20), sep=""), "Indet.", "indet.")


#Generate a pertables object (i.e. a list of biological data tables simulated from taxonomic
# uncertainty)
 \dontrun{
Amazonia100 <- pertables(Amazonia, index=index.Amazon, nsim=100)

# Assess the effects of taxonomic uncertainty on a RDA analysis of biological data explained
# by all the environmental variables of the soil data:

Amazonia.rda <- rda_pertables(Amazonia100 ~., data=soils)

Amazonia.rda

plot(Amazonia.rda)
 
}

# Fast example for Rcheck

Amazonia4.p2 <- pertables.p2(Amazonia[1:50,], index=index.Amazon, nsim=4, ncl=2, iseed=4)
set.seed(2)
Amazonia.rda <- rda_pertables(Amazonia4.p2 ~., data=soils)

Amazonia.rda

plot(Amazonia.rda)
 

}
\keyword{ multivariate }

