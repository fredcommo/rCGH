\name{byGeneTable}
\alias{byGeneTable}

\title{Converting a Segmentation Table Into a By-Gene Table}

\description{
This function creates a by-gene table by listing all the genes contained in each
of the segments in the segmentation table.\cr
Gene annotations (symbol, location,...), segmented Log2Ratios, and segment
legnth are reported in the final table.\cr
A supplementary score is the \code{relativeLog}: the magnitude, in Log2, from
the closest centromere.
}

\usage{
    byGeneTable(segTable, symbol = NULL,
        genome = c("hg19", "hg18", "hg38"), columns = NA, verbose = TRUE)
}

\arguments{
    \item{segTable}{
    : data frame. A segmentation table exported from an object of class
    \code{"\linkS4class{rCGH}"}
    }
    \item{symbol}{
    : character. A valid HUGO symbol. When \code{NULL} the full gene table is
    returned, the corresponding gene information only o/w.
    }
    \item{genome}{
    : string. The genome build to use. Supported genomes are \code{hg18},
    \code{hg19} (default), and \code{hg38}.
    }
    \item{columns}{
    : string. what supplementary genes annotations to export.
    Allowed annotations are those supported by the
    \code{\link[AnnotationDbi]{select}} method from the \code{AnnotationDbi}
    bioconductor package.\cr
    When \code{NA} (default), 'SYMBOL','ENTREZID', 'GENENAME' and 'MAP' are
    exported.
    }
    \item{verbose}{
    : logical. When \code{TRUE} progress is printed.
    }
}

\details{
For gene annotations, Hg19/GRCh37 annotations downloaded from \emph{NCBI} are
considered.
}

\value{
An object of class \code{"\linkS4class{rCGH}"}
}

\author{
Frederic Commo
}

\seealso{
\code{\link{getSegTable}}
}

\examples{

filePath <- system.file("extdata", "Affy_cytoScan.cyhd.CN5.CNCHP.txt.bz2",
    package = "rCGH")
cgh <- readAffyCytoScan(filePath, sampleName = "AffyScHD")
cgh <- adjustSignal(cgh, nCores=1)
cgh <- segmentCGH(cgh, nCores=1)
cgh <- EMnormalize(cgh)
st <- getSegTable(cgh)
bygene <- byGeneTable(st)
head(bygene)

}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
