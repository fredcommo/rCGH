\name{readIllumina}
\alias{readIllumina}

\title{
Illumina HumanOmniExpress \code{"\linkS4class{rCGH-Illumina}"} Constructor
}

\description{
A constructor function which takes an Illumina HumanOmniExpress .txt
file as input, possibly in a compressed format (.bz2 or .gz).\cr
These files are exported from GenomeStudio.
}

\usage{
readIllumina(filePath, sampleName=NA, labName=NA,
    genome = c("hg19", "hg18", "hg38"), verbose=TRUE)
}

\arguments{
    \item{filePath}{
    : string. A path to an Illumina .txt file.
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
    \item{verbose}{
    : logical. When \code{TRUE} (default), progress is printed.
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
\code{\link{readAgilent}}, \code{\link{readAffySNP6}},
\code{\link{readAffyCytoScan}}, \code{\link{readGeneric}},
\code{\link{setInfo}}, \code{\link{getInfo}}
}

\examples{

filePath <- system.file("extdata", "Illumina.txt.bz2", package = "rCGH")
cgh <- readIllumina(filePath, sampleName = "demo")
cgh

}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }

