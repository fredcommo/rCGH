\name{plotDensity}
\alias{plotDensity}
\alias{plotDensity-methods}
\alias{plotDensity,rCGH-method}

\title{Visualizing the Log2Ratios Density and Centralization Decision}

\description{
This function display the distribution of the Log2Ratios, as well as how the
\code{"\link{EMnormalize}"} step estimates the mixture of gaussian
populations, and choose a centralization value.
}

\usage{
\S4method{plotDensity}{rCGH}(object, breaks=NULL, Title=NULL,...)
}

\arguments{
    \item{object}{
    : An object of class \code{"\linkS4class{rCGH}"}
    }
    \item{breaks}{
    : The number of breaks to use. See \code{\link{hist}}.
    When \code{NULL} (default), \code{breaks} is arbitrarily defined from the
    number of values to draw.
    }
    \item{Title}{
    : character string. A title for the density plot. If \code{NULL} (default),
    the sample name (when exists) will be used.
    }
    \item{...}{
    : Other graphical parameters supported by \code{\link{par}}.
    }
}

\value{
None.
}

\author{
Frederic Commo
}

\seealso{
\code{\link{plotProfile}}, \code{\link{plotLOH}}, \code{\link{multiplot}},
\code{\link{view}}
}

\examples{

filePath <- system.file("extdata", "Affy_cytoScan.cyhd.CN5.CNCHP.txt.bz2",
    package = "rCGH")
cgh <- readAffyCytoScan(filePath, sampleName = "AffyScHD")
cgh <- adjustSignal(cgh, nCores=1)
cgh <- EMnormalize(cgh)
plotDensity(cgh)

}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
