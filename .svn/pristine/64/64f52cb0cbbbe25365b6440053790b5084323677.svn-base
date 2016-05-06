\name{rCGH-class}
\Rdversion{1.1}
\docType{class}
\alias{rCGH-class}

\title{Class \code{"rCGH"}}

\description{
Class \code{"\linkS4class{rCGH}"} is a superclass living on top of
\code{"\linkS4class{rCGH-Agilent}"}, \code{"\linkS4class{rCGH-SNP6}"}, and 
\code{"\linkS4class{rCGH-cytoScan}"}. These objects inherit most of the
properties of the superclass, but specific parametrizations used during the
analysis process.\cr
Objects are created by platform-specific \code{read} functions:
\code{"\link{readAgilent}"}, \code{"\link{readAffySNP6}"}, and
\code{"\link{readAffyCytoScan}"}, each corresponding to their matched file
format.\cr
Slots described below are used to store sample information and analysis
parameters, as well as segmentation results. All are accessible through specific
\code{"\link{Accessors}"} functions.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("rCGH", ...)}.\cr
Slots content are updated at each different analysis step, and are accessible
through specific \code{get} functions.
}

\section{Slots}{
    \describe{
        \item{\code{info}:}{
        Object of class \code{"character"}: where sample information can be 
        stored. See \code{"\link{getInfo}"} and \code{"\link{setInfo}"}.
        }
        \item{\code{cnSet}:}{
        Object of class \code{"data.frame"}: the full data set.
        See \code{"\link{getCNset}"}.
        }
        \item{\code{param}:}{
        Object of class \code{"list"}: the analysis parameters stored for
        traceability. See \code{"\link{getParam}"}.
        }
        \item{\code{segTable}:}{
        Object of class \code{"data.frame"}: the segmentation table.
        See \code{"\link{getSegTable}"}.
        }
    }
}

\section{Methods}{
    \describe{
        \item{show}{\code{signature(object = "rCGH")}: ...}
    }
}

\author{
Frederic Commo
}

\seealso{
\code{"\linkS4class{rCGH-Agilent}"}, \code{"\linkS4class{rCGH-SNP6}"},
\code{"\linkS4class{rCGH-cytoScan}"}
}

\examples{
showClass("rCGH")
}

\keyword{classes}
