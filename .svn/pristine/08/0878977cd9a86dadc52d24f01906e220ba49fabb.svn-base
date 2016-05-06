\name{rCGH-package}
\alias{rCGH-package}
\alias{rCGH}
\docType{package}

\title{
Comprehensive Pipeline for Analyzing and Visualizing Agilent and
    Affymetrix Array-Based CGH Data
}

\description{
A comprehensive pipeline for analyzing and interactively
visualizing genomic profiles generated through Agilent and Affymetrix
microarrays.\cr
As inputs, rCGH supports Agilent dual-color Feature Extraction files (.txt),
from 44 to 400K, and Affymetrix SNP6.0 and cytoScan probeset.txt, cychp.txt,
and cnchp.txt files, exported from ChAS or Affymetrix Power Tools.\cr
This package takes over all the steps required for a genomic profile analysis,
from reading the files to the segmentation and genes annotations, and provides
several visualization functions (static or interactive) which facilitate
profiles interpretation.\cr
Input files can be in compressed format, e.g. \code{.bz2} or \code{.gz}.
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
\code{\link{readAffyCytoScan}}
}

\examples{

filePath <- system.file("extdata", "Affy_cytoScan.cyhd.CN5.CNCHP.txt.bz2",
    package = "rCGH")
object1 <- readAffyCytoScan(filePath, sampleName = "AffyScHD")

object2 <- adjustSignal(object1, nCores=1)
object3 <- EMnormalize(object2)
object4 <- segmentCGH(object3, nCores=1)

# Static visalizations
plotDensity(object4)
multiplot(object4)

\dontrun{
# Interactive visalizations
view(object4)
}

}
