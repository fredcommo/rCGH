\name{readAffyCytoScan}
\alias{readAffyCytoScan}

\title{Affymetrix CytoScanHD \code{"\linkS4class{rCGH-cytoScan}"} Constructor}

\description{
A constructor function which takes an Affymetrix cytoScanHD cychp.txt
(or cnchp.txt) file as input, possibly in a compressed format (.bz2 or .gz).\cr
These files are exported from Chromosome Analysis Suite (ChAS) or Affymetrix
Power Tools (see \code{References} section).
}

\usage{
readAffyCytoScan(filePath, sampleName=NA, labName=NA,
    useProbes=c("snp", "cn", "all"), genome = c("hg19", "hg18", "hg38"),
    verbose=TRUE)
}

\arguments{
    \item{filePath}{
    : string. A path to an Affymetrix cytoScanHD cychp.txt (or cnchp.txt) file.
    }
    \item{sampleName}{
    : string. A sample Id. Optional.
    }
    \item{labName}{
    : string. A lab Id. Optional.
    }
    \item{useProbes}{
    : character. What probes to consider. Possible choices are SNP probes 
    only ("snp", default), CN probes only ("cn"), or all the probes ("all").
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

\references{
\href{http://www.affymetrix.com/estore/partners_programs/programs/developer/tools/powertools.affx}{Affymetrix Power Tools}
}

\author{
Frederic Commo
}

\seealso{
\code{\link{readAgilent}}, \code{\link{readAffySNP6}},
\code{\link{readGeneric}},
\code{\link{setInfo}}, \code{\link{getInfo}}
}

\examples{

filePath <- system.file("extdata", "Affy_cytoScan.cyhd.CN5.CNCHP.txt.bz2",
    package = "rCGH")
cgh <- readAffyCytoScan(filePath, sampleName = "AffyScHD")
cgh

}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }

