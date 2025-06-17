\name{pertables}
\alias{pertables}
\alias{pertables.p2}
\alias{print.pertables}

\title{
Function to incorporate the effect of taxonomic uncertainty on multivariate analyses of ecological data.
}
\description{
 This function implements a permutational method to incorporate taxonomic uncertainty on multivariate analyses typically used in the analysis of ecological data. The procedure is based on iterative randomizations that randomly re-assign non identified species in each site to any of the other species found in the remaining sites.
}
\usage{
pertables(data, index = NULL, nsim = 100)
pertables.p2(data, index = NULL, nsim = 100, ncl=2, iseed = NULL)
}

\arguments{
  \item{data}{
 Community data matrix. The three first columns are factors referring to the family, genus and species specific names. 
  The remaining columns are numeric vectors indicating species abundances at each site. 
}
  \item{index}{
List of additional parameters to determine the level at which species have been identified. Default values include 'Indet', 
  'indet', 'sp', 'sp1' to 'sp100', 'sp 1' to 'sp 100', '', and ' '.
}
  \item{nsim}{
Number of simulations of species' identities, i.e., number of data tables to simulate.
}
  \item{ncl}{
Number of clusters for parallel simulation.
}
\item{iseed}{
An integer to be supplied to clusterSetRNGStream, or NULL not to set reproducible seeds.
}
}
\details{
The procedure is implemented in two sequential steps:

Step 1. Morphospecies identified only to genus are randomly re-assigned with the same probability within the group of species and morphospecies that share the same genus, provided they are not found in the same sites. In the re-assignment of the species identity, the species considered can also receive its own identity. For instance, let's assume we have three floristic inventories. In site A we have \emph{Eugenia sp1} and \emph{E. nesiotica}. In site B we have \emph{Eugenia nesiotica}, \emph{E. principium} and \emph{E. salamensis}. In site C we have \emph{Eugenia sp2} and \emph{E. salamensis}. \emph{Eugenia sp1} can be thus re-identified with equal probability as \emph{Eugenia sp2, E. principium, E. salamensis} or just maintain its own identity (\emph{Eugenia sp1}). In the latter case, this means that we assume that \emph{E. sp1} is a completely different species, although we do not know its true identity. On the contrary, we cannot re-identify \emph{E. sp1} as \emph{E. nesiotica} because they were found in the same site, so we are quite certain that \emph{E. sp1} is different from \emph{E. nesiotica}. The same is applied to species identified only to family and fully unidentified species. Note that when collating inventories from different researchers, we must rename all unidentified species. This is because two researchers can use the same label, e.g. \emph{Eugenia sp1}, even though this name does not necessarily refer to the same species. For a verification of the biological identity of \emph{Eugenia sp1} one would need to cross-check the vouchers bearing the same name.


Step 2. Step 1 is iterated \emph{nsim} times. As a result, \emph{nsim} matrices are obtained, all of which contain the same number of sites but variable number of species depending on the resulting re-assignment of morphospecies, The process can be time-consuming if community data matrices are large.

Function \code{pertables.p2} implements a parallelized version which considerably reduces computation time.

}
\value{
The function return a list of class \code{pertables} with the following components

\item{taxunc }{Summary of the number of species fully identified (0), identified to genus (1), identified to family (2), or fully undetermined (3).}
\item{pertables }{A list with all the simulated data matrices.}
\item{raw}{The raw data matrix, without the unidentified especies.}

}
\references{
Cayuela, L., De la Cruz, M. and Ruokolainen, K. (2011). A method to incorporate the effect of taxonomic uncertainty on multivariate analyses of ecological data. \emph{Ecography}, 34: 94-102. \doi{10.1111/j.1600-0587.2009.05899.x}.   
}
\author{
Luis Cayuela and Marcelino de la Cruz
}

\examples{

data(Amazonia)
data(soils)

# Define a new index that includes the terms used in the \code{Amazonia} dataset to define
# undetermined taxa at different taxonomic levels

index.Amazon <- c(paste("sp.", rep(1:20), sep=""), "Indet.", "indet.")

#Generate a pertables object (i.e. a list of biological data tables simulated from taxonomic
# uncertainty)
 
 \dontrun{
# compare prformance of pertables and pertables.p2
nsim <-100
ncl <-2
gc()
t0<- Sys.time()
 Amazonia100<- pertables(Amazonia, index=index.Amazon, nsim=nsim)
 Sys.time()-t0
gc()
t0<- Sys.time()
 Amazonia100.p2<- pertables.p2(Amazonia, index=index.Amazon, nsim=nsim, ncl=ncl)
 Sys.time()-t0
}
# Example for Rcheck

Amazonia4.p2<- pertables.p2(Amazonia, index=index.Amazon, nsim=4, ncl=2)

}
\keyword{ multivariate }

