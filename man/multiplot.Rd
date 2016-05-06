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
\S4method{multiplot}{rCGH}(object, symbol=NULL, gain=.5,
    loss=(-.5), minLen=10, pCol = "grey50", GLcol = c("blue", "red3"),
    L=matrix(seq(1, 12)), p=c(1/2, 1/4, 1/4), Title=NULL, ylim=NULL)
}

\arguments{
    \item{object}{
    : An object of class \code{"\linkS4class{rCGH}"}
    }
    \item{symbol}{
    : character. A valid HUGO symbol (case insensitive).
    }
    \item{gain}{
    : numeric. A gain threshold value (in \eqn{Log2(Ratio)}). Segments greater,
    or equal to, this value will be colored, as specified by \code{GLcol}.
    }
    \item{loss}{
    : numeric.  A loss threshold value (in \eqn{Log2(Ratio)}). Segments lower,
    or equal to, this value will be colored, as specified by \code{GLcol}.
    }
    \item{minLen}{
    : numeric. The mininal length for a segment, expressed in Kb.
    When \code{NULL} (default), segments are reported as they have been
    computed by \code{\link{segmentCGH}}. Segments shorter than the
    specified value are re-merged otherwise.
    }
    \item{pCol}{
    : string. The probe points color. Default is \code{"grey50"}.
    }
    \item{GLcol}{
    : vector. A vector of 2 colors: the gained and lost segments colors,
    respectively. Default is \code{"blue"} for gains and \code{"red3"} for
    losses.
    }
    \item{L}{
    : matrix. A matrix defining the layout. Default is 12 lines.
    }
    \item{p}{
    : numeric. The proportion of each plot within the plot window.
    Default is \code{1/2, 1/4, 1/4}, which corresponds to 6-4-4 lines for the 
    genomic profile in Log2R, in copy number, and the LOH plot, respectiviely,
    and given a 12-line layout.
    }
    \item{Title}{
    : character string. A title for the plot. If \code{NULL} (default),
    the sample name (when exists) is used.
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
cgh <- segmentCGH(cgh, nCores=1)
cgh <- EMnormalize(cgh)

# Static visalizations
multiplot(cgh, symbol = "erbb2")

}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
