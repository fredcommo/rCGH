\name{segmentCGH}
\alias{segmentCGH}
\alias{segmentCGH-methods}
\alias{segmentCGH,rCGH-method}

\title{Genomic Profile Segmentation}

\description{
A function for performing the Log2Ratio segmentation on an object of class
\code{"\linkS4class{rCGH}"}. See the \code{details} section below.
}

\usage{
\S4method{segmentCGH}{rCGH}(object, Smooth=TRUE, UndoSD = NULL,
minLen = 10, nCores=NULL, verbose = TRUE)
}

\arguments{
    \item{object}{
    : An object of class \code{"\linkS4class{rCGH}"}
    }
    \item{Smooth}{
    : logical. Should the LRR be smoothed before being segmented.
    See \code{DNAcopy} for details.
    }
    \item{UndoSD}{
    : numeric. When not specified (default is \code{NULL}), the UndoSD value is
    estimated from the Log2Ratios. See \code{DNAcopy} for details.
    }
    \item{minLen}{
    : numeric. The minimal length for a segment, in Kb. Shorter segments
    will be merged to the closest adjacent one. Default value is 10(Kb).
    }
    \item{nCores}{
    : numeric. The number of cores to use in order to speed up the computation.
    When \code{NULL} (default), half of the available cores are used.
    See \code{\link[parallel]{mclapply}}.
    }
    \item{verbose}{
    : logical. if \code{TRUE} (default), progress is printed.
    }
}

\details{
This function is a wrapper for the \code{DNAcopy}, \code{\link[DNAcopy]{CNA}}
and \code{\link[DNAcopy]{segment}} functions, which allows parallelization and
data-driven parameterization.\cr
In addition to the usual \code{DNAcopy} output, the segmentation table contains
the probes Log2Ratio standard deviation for each segment, as well as there
length, in Kb.
}

\value{
An object of class \code{"\linkS4class{rCGH}"}
}

\references{
\href{http://www.ncbi.nlm.nih.gov/pubmed/17234643}{Venkatraman ES1, Olshen AB.
A faster circular binary segmentation algorithm for the analysis of array CGH
data.Bioinformatics. 2007 Mar 15;23(6):657-63.}
}

\author{
Frederic Commo
}

\seealso{
\code{\link[DNAcopy]{CNA}}, \code{\link[DNAcopy]{segment}},
\code{\link[parallel]{mclapply}}
}

\examples{

filePath <- system.file("extdata", "Affy_cytoScan.cyhd.CN5.CNCHP.txt.bz2",
    package = "rCGH")
cgh <- readAffyCytoScan(filePath, sampleName = "AffyScHD")
cgh <- adjustSignal(cgh, nCores=1)
cgh <- EMnormalize(cgh)
cgh <- segmentCGH(cgh, nCores=1)
st <- getSegTable(cgh)
head(st)

}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
