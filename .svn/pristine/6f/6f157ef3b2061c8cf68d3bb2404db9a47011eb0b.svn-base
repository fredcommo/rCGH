\name{rCGH-package}
\alias{rCGH-package}
\alias{rCGH}
\docType{package}

\title{
Comprehensive Pipeline for Analyzing and Visualizing Array-Based CGH Data
}

\description{
    A comprehensive pipeline for analyzing and interactively visualizing genomic
    profiles generated through commercial or custom aCGH arrays. As inputs, rCGH
    supports Agilent dual-color Feature Extraction files (.txt), from 44 to
    400K, Affymetrix SNP6.0 and cytoScanHD probeset.txt, cychp.txt, and
    cnchp.txt files exported from ChAS or Affymetrix Power Tools.\cr
    rCGH also supports custom arrays, provided data is
    in a suitable format. This package takes over all the steps required for
    individual genomic profiles analysis, from reading files to segmenting
    and annotating genes. This package provides several visualization functions
    (static or interactive) which facilitate individual profiles interpretation.
    Input files can be in compressed format, e.g. .bz2 or .gz.
}

\author{
Frederic Commo <frederic.commo@gustaveroussy.fr>
}

\keyword{package}
\keyword{aCGH}
\keyword{CopyNumberVariation}
\keyword{Preprocessing}
\keyword{FeatureExtraction}

\seealso{
\code{\link{readAgilent}}, \code{\link{readAffySNP6}},
\code{\link{readAffyCytoScan}}, \code{\link{readGeneric}}
}

\examples{

filePath <- system.file("extdata", "Affy_cytoScan.cyhd.CN5.CNCHP.txt.bz2",
    package = "rCGH")
cgh <- readAffyCytoScan(filePath, sampleName = "AffyScHD")

cgh <- adjustSignal(cgh, nCores=1)
cgh <- EMnormalize(cgh)
cgh <- segmentCGH(cgh, nCores=1)

# Static visalizations
plotDensity(cgh)
multiplot(cgh)

\dontrun{
# Interactive visalizations
view(cgh)
}

}
