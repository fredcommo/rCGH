\name{view}
\alias{view}
\alias{view-methods}
\alias{view,rCGH-method}

\title{Interactive Genomic Profile Visualization}

\description{
    This function is build on top of \code{\link[shiny]{shiny}}, and provides an
    iteractive way for visualizing a genomic profile, and exploring the list of
    genes.\cr
    From a command panel, the user can interact with the graph in different 
    ways. See \code{details}.
}

\usage{
\S4method{view}{rCGH}(object, browser = TRUE, ...)
}

\arguments{
    \item{object}{
    : An object of class \code{"\linkS4class{rCGH}"}
    }
    \item{browser}{
    : logical. When \code{TRUE} (default), the system's default web browser will
    be launched automatically.
    }
    \item{...}{
    : Optional parameters used by \code{\link[shiny]{runApp}}.
    }
}

\details{
    The left command panel allows the user several actions:
    \itemize{
        \item displaying a specific gene by calling its HUGO symbol.
        \item showing all or one unique chromosome.
        \item merging segments shorter than a specified value, in Kb.
        \item recentering the entire profile.
        \item rescaling the y-axis.
        \item specifying the Log2Ratio cut offs for defining gains and losses.
        \item specifying a segment lenght cut off, in Mb.
        \item exporting the genomic plot.
        \item exporting the genes list.
    }
    Some actions, such as showing one unique chromosome or specifying cut offs
    (gain, loss, segment length), automatically update the gene table available
    in the \emph{"Genes table"} tab.
}

\value{
    None.
}

\author{
Frederic Commo
}

\seealso{
\code{\link{plotProfile}}, \code{\link{plotLOH}}, \code{\link{multiplot}},
\code{\link[shiny]{runApp}}
}

\examples{

filePath <- system.file("extdata", "Affy_cytoScan.cyhd.CN5.CNCHP.txt.bz2",
    package = "rCGH")
cgh <- readAffyCytoScan(filePath, sampleName = "AffyScHD")
cgh <- adjustSignal(cgh, nCores=1)
cgh <- EMnormalize(cgh)
cgh <- segmentCGH(cgh, nCores=1)

\dontrun{
# Interactive visalizations
view(cgh)
}

}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
