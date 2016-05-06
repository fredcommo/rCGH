\name{hg18}
\alias{hg18}
\docType{data}
\title{Hg18 Chromosome Lengths and Centromere Locations}
\description{
A data set containing lengths and centromere locations for each of the 24
chromosomes, according to Hg18.
}
\usage{hg18}
\value{a data set.}
\format{
A data set with 24 rows and 5 columns:
    \itemize{
        \item chrom: chromosome number.
        \item length: chromosome length.
        \item centromerStart: centromere start position.
        \item centromerEnd: centromere start position.
        \item cumlen: cumulative length (where the previous chromosome ends).
    }
}
\source{
These data derived from the Hg18 gap UCSC table, freely available at:
\href{http://genome.ucsc.edu/cgi-bin/hgTables}{UCSC}\cr
Access date: 10-10-2015\cr
Within the browser, select:\cr
group: All Tables\cr
database: hg18\cr
table: gap
}
\examples{
# For users convenience, we provide a prebuilt dataset
# containing the Hg18 chr lengths, and centromeres location.
hg18

# The same dataset can be obtained as follow:
\dontrun{
library(BSgenome)
library(rtracklayer)

getChrLength <- function(genome){
    genome <- sprintf("BSgenome.Hsapiens.UCSC.%s", genome)
    g <- getBSgenome(genome, masked=FALSE)
    data.frame(chrom=1:24, length=seqlengths(g)[1:24])
}
.chrAsNum <- function(tbl){
    tbl$chrom <- gsub("chr", "", tbl$chrom)
    tbl$chrom[tbl$chrom=="X"] <- 23
    tbl$chrom[tbl$chrom=="Y"] <- 24
    tbl$chrom <- as.numeric(tbl$chrom)
    tbl[order(tbl$chrom),]
}
getCentromeres <- function(genome){
    mySession <- try(browserSession("UCSC"), silent=TRUE)
    # In case it fails, use another mirror
    if(inherits(mySession, "try-error"))
        mySession <- browserSession("UCSC",
                                    url="http://genome-euro.ucsc.edu/cgi-bin/")
    genome(mySession) <- genome
    obj <- ucscTableQuery(mySession, table="gap")
    tbl <- getTable(obj)
    tbl <- tbl[tbl$type=="centromere", c("chrom", "chromStart", "chromEnd")]
    colnames(tbl)[2:3] <- c("centromerStart", "centromerEnd")
    .chrAsNum(tbl)
}
makeHg <- function(genome){
    chrL <- getChrLength(genome)
    ctm <- getCentromeres(genome)
    tbl <- merge(chrL, ctm, by="chrom")
    cumlen <- c(0, cumsum(as.numeric(tbl$length))[-nrow(tbl)])
    cbind.data.frame(tbl, cumlen=cumlen)    
}
hg18 <- makeHg("hg18")
hg18
}
}
\author{Frederic Commo}
\keyword{datasets}
