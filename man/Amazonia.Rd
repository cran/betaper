\name{Amazonia}
\alias{Amazonia}
\alias{soils}
\docType{data}
\title{ Tree abundance and soil data in Western Amazonia}
\description{
  The \code{Amazonia} data frame has tree counts in nine 0.16-hectare inventory plots in Western Amazonia. \code{soils} contains data on soil
  cations at each location.
}
\usage{
data(Amazonia)
data(soils)
}
\format{
  \code{Amazonia} is a data frame with 1188 observations (species) and 12 columns (taxonomic descripcion and sites). The three first columns 
  refer to family, genus and specific species Latin names. Columns 4 to 12 have tree abundance data for nine inventory plots.
  
  \code{soils} is a data frame with 9 observations (inventory plots) and 4 columns (variables). Soil variables (Ca, K, Mg, Na) are given in cmol/kg.
}
\details{
  Data from Western Amazonia includes tree inventories at nine lowland sites (approximately 100-150 m above sea level) near Iquitos, Peru. 
  The sites were selected to represent regional variations in geology and were distributed along a soil nutrient gradient ranging from poor 
  loamy soils to richer clayey soils. Each inventory consisted of 20 x 20 m plots (0.16 ha total area) distributed along 1.3-km transects. 
  At each site, K. Ruokolainen and colleagues identified to species or morphospecies all woody, free-standing stems of > 2.5 cm dbh. The 
  full inventories sampled 3980 individuals from 1188 species or morphospecies.
}
\references{
  Higgins, M.A. & Ruokolainen, K. 2004. Rapid tropical forest inventory: a comparison of techniques based on inventory data from western 
  Amazonia. \emph{Conservation Biology} 18(3): 799-811.
  
  Ruokolainen, K., Tuomisto, H., Macia, M.J., Higgins, M.A. & Yli-Halla, M. 2007. Are floristic and edaphic patterns in Amazonian rain forests 
  congruent for trees, pteridophytes and Melastomataceae? \emph{Journal of Tropical Ecology} 23: 13-25.
}
\examples{
data(Amazonia)
data(soils)
}
\keyword{datasets}
