\name{readAgilent}
\alias{readAgilent}

\title{
Agilent Dual-Color Hybridization \code{"\linkS4class{rCGH-Agilent}"}
    Constructor.
}

\description{
A constructor function taking as input an Agilent FE .txt file, exported from
Feature Extraction, possibly in a compressed format (.bz2 or .gz).\cr
Agilent from 44 to 400K are supported.
}

\usage{
readAgilent(filePath, sampleName = NA,
    labName = NA, supFlags = TRUE, verbose = TRUE)
}

\arguments{
    \item{filePath}{
    : string. A path to an Agilent FE (.txt) file.
    }
    \item{sampleName}{
    : string. A sample Id. Optional.
    }
    \item{labName}{
    : string. A lab Id. Optional.
    }
    \item{supFlags}{
    : should the flagged probes be suppressed. Default is \code{TRUE}.
    }
    \item{verbose}{
    : logical. if \code{TRUE} (default), progress is printed.
    }
}

\details{
When available in the file preambule, several array information will be stored
in \code{Object@info}: scanning date, grid version,...\cr
Any other useful item can be stored using \code{\link{setInfo}}.
}

\value{
An object of class \code{"\linkS4class{rCGH}"}
}

\author{
Frederic Commo
}

\seealso{
\code{\link{readAffyCytoScan}}, \code{\link{readAffySNP6}},
\code{\link{setInfo}}, \code{\link{getInfo}}
}

\examples{

filePath <- system.file("extdata", "Agilent4x180K.txt.bz2", package = "rCGH")
cgh <- readAgilent(filePath, sampleName = "Agilent4x180K", labName = "myLab")
cgh

}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
