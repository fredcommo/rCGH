\name{recenter<--methods}
\docType{methods}
\alias{recenter}
\alias{recenter<-}
\alias{recenter<--methods}
\alias{recenter<-,rCGH-method}

\title{ Recentering a Genomic Profile }

\description{
This function allows the user to recenter a genomic profile stored in
an object of class \code{"\linkS4class{rCGH}"}.\cr
Peaks are indexed from 1 to k, from left to right, as they appear on the 
\code{\link{plotDensity}} after the \code{\link{EMnormalize}} step.
}

\usage{
\S4method{recenter}{rCGH}(object) <- value
}

\arguments{
    \item{object}{
    : An object of class \code{"\linkS4class{rCGH}"}
    }
    \item{value}{
    : numeric. What peak number to choose to recenter the genomic profile.
    }
}

\value{
An object of class \code{"\linkS4class{rCGH}"}
}

\author{
Frederic Commo
}

\note{
When a profile is recentered, the stored workflow parameters are updated.
see \code{\link{getParam}}.
}

\seealso{
\code{\link{EMnormalize}}, \code{\link{plotDensity}}
}

\examples{

filePath <- system.file("extdata", "Affy_cytoScan.cyhd.CN5.CNCHP.txt.bz2",
    package = "rCGH")
cgh <- readAffyCytoScan(filePath, sampleName = "AffyScHD")
cgh <- adjustSignal(cgh, nCores=1)
cgh <- EMnormalize(cgh)
cgh <- segmentCGH(cgh, nCores=1)

# Default peak choice center the profile on the 2nd peak
plotDensity(cgh)

# Recentering on the 3rd density peak
recenter(cgh) <- 3
plotDensity(cgh)

}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
