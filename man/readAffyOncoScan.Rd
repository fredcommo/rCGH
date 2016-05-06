\name{readAffyOncoScan}
\alias{readAffyOncoScan}

\title{Affymetrix OncoScan \code{"\linkS4class{rCGH-oncoScan}"} Constructor}

\description{
A constructor function which takes an Affymetrix oncoScan .txt
file as input, possibly in a compressed format (.bz2 or .gz).\cr
This can be either a \code{'ProbeSets,CopyNumber.tsv'} alone, or merged with its
corresponding \code{'ProbeSets,AllelicData.tsv'} file. See the \code{details} 
section.
}

\usage{
readAffyOncoScan(filePath, sampleName=NA, labName=NA,
    genome = c("hg19", "hg18", "hg38"),
    ploidy = 2, verbose=TRUE)
}

\arguments{
    \item{filePath}{
    : string. A path to an Affymetrix cychp.txt (or cnchp.txt) file.
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
    : numeric. A priori plody value, when known, to adjust the estimation
    of copy numbers. Default is 2.
    }
    \item{verbose}{
    : logical. When \code{TRUE} (default), progress is printed.
    }
}

\details{
Merging the \code{'ProbeSets,CopyNumber.tsv'} with the
\code{'ProbeSets,AllelicData.tsv'} may produce a file containing both probes
\code{Log2Ratio} and \code{AllelicDifference}. 
}

\value{
An object of class \code{"\linkS4class{rCGH}"}
}

\references{
\href{http://www.affymetrix.com/estore/partners_programs/programs/developer/tools/powertools.affx}{Affymetrix Power Tools}
}

\author{
Frederic Commo
}

\seealso{
\code{\link{readAgilent}}, \code{\link{readAffySNP6}},
\code{\link{readGeneric}}, \code{\link{readAffyCytoScan}},
\code{\link{setInfo}}, \code{\link{getInfo}}
}

\examples{

# Just a toy file
filePath <- system.file("extdata", "oncoscan.tsv.bz2", package = "rCGH")
cgh <- readAffyOncoScan(filePath, sampleName = "AffyOncoScan")
cgh

}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }

