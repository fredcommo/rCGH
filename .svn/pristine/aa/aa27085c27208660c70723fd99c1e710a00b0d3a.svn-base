\name{plotLOH}
\alias{plotLOH}
\alias{plotLOH-methods}
\alias{plotLOH,rCGH-method}

\title{Allelic Differences Visualization}

\description{
This function display a static view of the allele differences, when available.
}

\usage{
\S4method{plotLOH}{rCGH}(object, Title=NULL)
}

\arguments{
    \item{object}{
    : An object of class \code{"\linkS4class{rCGH}"}
    }
    \item{Title}{
    : character string. A title for the density plot. If \code{NULL} (default),
    the sample name (when exists) is used.
    }
}

\value{
None.
}

\author{
Frederic Commo
}

\seealso{
\code{\link{plotDensity}}, \code{\link{plotProfile}}, \code{\link{multiplot}},
\code{\link{view}}
}

\examples{

filePath <- system.file("extdata", "Affy_cytoScan.cyhd.CN5.CNCHP.txt.bz2",
    package = "rCGH")
cgh <- readAffyCytoScan(filePath, sampleName = "AffyScHD")
cgh <- adjustSignal(cgh, nCores=1)
cgh <- EMnormalize(cgh)
cgh <- segmentCGH(cgh, nCores=1)

# Static visalizations
plotLOH(cgh)

}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
