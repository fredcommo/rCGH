\name{readGeneric}
\alias{readGeneric}

\title{Generic rCGH object \code{"\linkS4class{rCGH-generic}"} Constructor}

\description{
A constructor function which takes a tabulated .txt file as input, possibly in a
compressed format (.bz2 or .gz).
Notice that precise column names are mandatory, see the \code{details} section.
}

\usage{
readGeneric(filePath, sampleName=NA, labName=NA,
    genome = c("hg19", "hg18", "hg38"), ploidy = 2, verbose=TRUE)
}

\arguments{
    \item{filePath}{
    : string. A path to an Generic .txt file.
    }
    \item{sampleName}{
    : string. A sample Id. Optional.
    }
    \item{labName}{
    : string. A lab Id. Optional.
    }
    \item{genome}{
    : string. The genome build to use. Supported genomes are \code{hg18},
    \code{hg19} (default), and \code{hg38}.
    }
    \item{ploidy}{
    : numeric. A priori ploidy value, when known, to adjust the estimation
    of copy numbers. Default is 2.
    }
    \item{verbose}{
    : logical. When \code{TRUE} (default), progress is printed.
    }
}

\details{
This generic constructor does not expect any preamble.
Mandatory columns are:
    \describe{
        \item{\code{ProbeName}:}{
        Character strings. Typicaly the probe ids.
        }
        \item{\code{ChrNum}:}{
        numeric. The chromosome numbers. In case Chr X and Y are used and named
        as \code{"X"} and \code{"Y"}, these notations will be converted into
        23 and 24, respectively.
        }
        \item{\code{ChrStart}:}{
        numeric. The chromosomal probes locations.
        }
        \item{\code{Log2Ratio}:}{
        numeric. The corresponding Log2Ratios.
        }
    }

}

\value{
An object of class \code{"\linkS4class{rCGH}"}
}


\author{
Frederic Commo
}

\seealso{
\code{\link{readAgilent}}, \code{\link{readAffySNP6}},
\code{\link{readAffyCytoScan}}, \code{\link{readAffyOncoScan}},
\code{\link{setInfo}}, \code{\link{getInfo}}
}

\examples{

filePath <- system.file("extdata", "generic.txt.bz2", package = "rCGH")
cgh <- readGeneric(filePath, sampleName = "demo")
cgh

}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }

