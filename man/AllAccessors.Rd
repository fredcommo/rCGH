\name{AllAccessors}
\alias{AllAccessors}
\alias{Accessors}
\alias{getInfo}
\alias{getInfo-methods}
\alias{getInfo,rCGH-method}
\alias{getCNset}
\alias{getCNset-methods}
\alias{getCNset,rCGH-method}
\alias{getParam}
\alias{getParam-methods}
\alias{getParam,rCGH-method}
\alias{getSegTable}
\alias{getSegTable-methods}
\alias{getSegTable,rCGH-method}

\title{ \code{"\linkS4class{rCGH}"} Accessor Functions }

\description{
Methods for extracting information from an object of
class \code{"\linkS4class{rCGH}"}.\cr
Each of the below methods are simply convenience functions which extract the
corresponding slots (as the name of each method suggests) from an object of
class \code{"\linkS4class{rCGH}"}.
}

\usage{
\S4method{getInfo}{rCGH}(object, item = NULL)
\S4method{getCNset}{rCGH}(object)
\S4method{getParam}{rCGH}(object)
\S4method{getSegTable}{rCGH}(object, minLen = NULL)
}

\arguments{
    \item{object}{
    : An object of class \code{"\linkS4class{rCGH}"}
    }
    \item{item}{
    : character. Can be one, or a vector of items.
    When \code{NULL}, the full available information is returned.
    If \code{item} is specified, and exists, the corresponding value(s) only
    is(are) returned.
    }
    \item{minLen}{
    : numeric. The mininal length for a segment, in Kb.
    When \code{NULL} (default), the segmentation table is exported, as it has
    been computed with \code{\link{segmentCGH}}, segments shorter than the
    specified value are re-merged otherwise.
    }
}

\section{Methods}{
    \describe{
        \item{\code{"\linkS4class{rCGH}"}}{
        \itemize{
            \item getInfo(object, item = NULL): returns the values of the 
            specified items, all the information otherwise.
            \item getCNset(object): returns the full by-probe dataset.
            \item getParam(object): returns the analysis parameters.
            \item getSegTable(object, minLen = NULL): returns the segmentation
            table - one row per segment.
            }
        }
    }
}

\value{
    \itemize{
        \item getInfo(object, item = NULL): character.
        \item getCNset(object): a data frame.
        \item getParam(object): a list of parameters.
        \item getSegTable(object, minLen = NULL): a data frame.
    }
}

\author{
Frederic Commo
}

\seealso{
\code{\link{setInfo}}, \code{\link{segmentCGH}}
}

\examples{

filePath <- system.file("extdata", "Agilent4x180K.txt.bz2", package = "rCGH")
cgh <- readAgilent(filePath, sampleName = "Agilent4x180K", labName = "myLab")

# Getting all the information
getInfo(cgh)

# Getting specific items
getInfo(cgh, c("sampleName", "labName"))

}

