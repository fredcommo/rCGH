\name{plotProfile}
\alias{plotProfile}
\alias{plotProfile-methods}
\alias{plotProfile,rCGH-method}

\title{Static Genomic Profile Visualization}

\description{
This function display a static view of the genomic profile stored in an
object of class \code{"\linkS4class{rCGH}"}.
}

\usage{
\S4method{plotProfile}{rCGH}(object, showCopy=FALSE, symbol=NULL,
    gain=.5, loss=(-.5), minLen = 10, pCol = "grey50",
    GLcol = c("blue", "red3"), Title=NULL, ylim=NULL)
}

\arguments{
    \item{object}{
    : An object of class \code{"\linkS4class{rCGH}"}
    }
    \item{showCopy}{
    : logical. To show the estimated copy numbers instead of the Log2Ratios.
    default is \code{FALSE}.
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
    \item{pCol}{
    : string. The probe points color. DEfault is \code{"grey50"}.
    }
    \item{GLcol}{
    : vector. A vector of 2 colors: the gained and lost segments colors,
    respectively. Default is \code{"blue"} for gains and \code{"red3"} for
    losses.
    }
    \item{Title}{
    : character string. A title for the density plot. If \code{NULL} (default),
    the sample name (when exists) is used.
    }
    \item{ylim}{
    : numeric. A vector of two values specifying a range for the y-axis.
    If \code{NULL} (default), the range of Log2Ratio is used.
    }
}

\value{
None.
}

\author{
Frederic Commo
}

\seealso{
\code{\link{plotDensity}}, \code{\link{plotLOH}}, \code{\link{multiplot}},
\code{\link{view}}
}

\examples{

filePath <- system.file("extdata", "Affy_cytoScan.cyhd.CN5.CNCHP.txt.bz2",
    package = "rCGH")
cgh <- readAffyCytoScan(filePath, sampleName = "AffyScHD")
cgh <- adjustSignal(cgh, nCores=1)
cgh <- segmentCGH(cgh, nCores=1)
cgh <- EMnormalize(cgh)

# Static visalization using Log2Ratios
plotProfile(cgh, symbol = "erbb2")

# Static visalization using estimated copy numbers
plotProfile(cgh, showCopy = TRUE, symbol = "erbb2")

}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
