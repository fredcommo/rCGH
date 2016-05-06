\name{multiplot}
\alias{multiplot}
\alias{multiplot-methods}
\alias{multiplot,rCGH-method}

\title{Static Genomic Profile and LOH Visualization}

\description{
This function display a static view of the genomic profile and the allelic
difference stored in an object of class \code{"\linkS4class{rCGH}"}.\cr
If no allelic difference is available, the genomic profile only is displayed.
}

\usage{
\S4method{multiplot}{rCGH}(object, symbol=NULL, gain=.5, loss=(-.5), minLen=10,
L=matrix(seq(1, 12)), p=c(2/3, 1/3), Title=NULL, ylim=NULL)
}

\arguments{
    \item{object}{
    : An object of class \code{"\linkS4class{rCGH}"}
    }
    \item{symbol}{
    : character. A valid HUGO symbol (case insensitive).
    }
    \item{gain}{
    : numeric. A gain threshold value (in \eqn{Log2(Ratio)}) from where
    gained segments will be shown, in blue.
    }
    \item{loss}{
    : numeric.  A loss threshold value (in \eqn{Log2(Ratio)}) from where
    lossed segments will be shown, in red.
    }
    \item{minLen}{
    : numeric. The mininal length for a segment, in Kb.
    When \code{NULL} (default), segments are reported as they have been
    computed with \code{\link{segmentCGH}}, segments shorter than the
    specified value are re-merged otherwise.
    }
    \item{L}{
    : matrix. A matrix defining how the layout is built. Default is 12 lines.
    }
    \item{p}{
    : numeric. The proportion of each plot in the plot window.
    Default is \code{2/3, 1/3}, which corresponds to 8 and 4 lines for the 
    genomic profile and the LOH plot, respectiviely.
    }
    \item{Title}{
    : character string. A title for the density plot. If \code{NULL} (default),
    the sample name (when exists) will be used.
    }
    \item{ylim}{
    : numeric. A vector of two values specifying the y-axis range.
    See \code{\link{plotProfile}}.
    }
}

\value{
None.
}

\note{
If no allelic difference is available, the genomic profile only is displayed.
}

\author{
Frederic Commo
}

\seealso{
\code{\link{plotDensity}}, \code{\link{plotProfile}},
\code{\link{plotLOH}}, \code{\link{view}}
}

\examples{

filePath <- system.file("extdata", "Affy_cytoScan.cyhd.CN5.CNCHP.txt.bz2",
    package = "rCGH")
cgh <- readAffyCytoScan(filePath, sampleName = "AffyScHD")
cgh <- adjustSignal(cgh, nCores=1)
cgh <- EMnormalize(cgh)
cgh <- segmentCGH(cgh, nCores=1)

# Static visalizations
multiplot(cgh, "erbb2")

}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
