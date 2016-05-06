\name{geneDB}
\alias{geneDB}
\docType{data}

\title{Gene Annotations, in Hg19}

\description{
GRanges object with 23056 ranges and 1 metadata column, according to
UCSC Hg19.
}

\usage{geneDB}

\value{A GRanges object}

\format{
A GRanges object including:
    \itemize{
        \item seqnames, class 'Rle': chromosome number.
        \item ranges:, class 'IRanges': gene [start, end] position.
        \item strand, class 'Rle': strand.
        \item elementMetadata, gene_id, class 'DataFrame': entrezIds.
    }
}

\note{
    This object is automaticaly generated from
    \code{\link{TxDb.Hsapiens.UCSC.hg19.knownGene}} when
    \code{\link{rCGH}} is loaded.
}

\author{Frederic Commo}

\seealso{
\code{\link{TxDb.Hsapiens.UCSC.hg19.knownGene}}
\code{\link[GenomicFeatures]{genes}}
}

\examples{
geneDB

# Below, how rCGH builds geneDB:
\dontrun{
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library( org.Hs.eg.db)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
geneDB  <- genes(txdb, columns=c("gene_id"))
geneDB
}
}

\keyword{datasets}
