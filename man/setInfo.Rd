\name{setInfo<--methods}
\docType{methods}
\alias{setInfo}
\alias{setInfo<-}
\alias{setInfo<--methods}
\alias{setInfo<-,rCGH-method}

\title{
Adding Information In An Object Of Class \code{"\linkS4class{rCGH}"}
}

\description{
This function allows the user to store any type of supplementaty information in
an object of class \code{"\linkS4class{rCGH}"}.
}

\usage{
\S4method{setInfo}{rCGH}(object, item = NULL) <- value
}

\arguments{
    \item{object}{
    : An object of class \code{"\linkS4class{rCGH}"}
    }
    \item{item}{
    : character. An item name to store. Default is \code{NULL}.
    }
    \item{value}{
    : any. A value to store.
    }
}

\value{
An object of class \code{"\linkS4class{rCGH}"}
}

\author{
Frederic Commo
}

\section{Warning}{
When either \code{item} or \code{value} is \code{NULL}, an error is returned.
}

\seealso{
\code{\link{getInfo}}
}

\examples{

filePath <- system.file("extdata", "Affy_cytoScan.cyhd.CN5.CNCHP.txt.bz2",
    package = "rCGH")
cgh <- readAffyCytoScan(filePath, sampleName = "AffyScHD")

# When supplementary information is added,
# numerical, logical, or strings are supported
setInfo(cgh, "someItem1") <- 35
setInfo(cgh, "someItem2") <- TRUE
setInfo(cgh, "someItem3") <- "someComment"
getInfo(cgh)

# or to get back specific items
getInfo(cgh, c("someItem1", "someItem3"))

}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
