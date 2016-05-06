\name{adjustSignal}
\alias{adjustSignal}
\alias{adjustSignal-methods}
\alias{adjustSignal,rCGH-method}

\title{Array-based CGH Preprocessing}

\description{
This function performs several preprocessing steps: local regressions
(loessFit) are used to correct cy3/cy5 and GC\% bias, when Agilent dual-color
hybridization are used only.\cr
In case of Affymetrix cychp.txt (or cnchp.txt) data are used, these steps have
been already done in ChAS or Affymetrix Power Tools.
}

\usage{
\S4method{adjustSignal}{rCGH}(object, Scale=TRUE, Cy=TRUE, GC=TRUE, Ref="cy3",
    suppOutliers=TRUE, nCores=NULL, verbose=TRUE)
}

\arguments{
    \item{object}{
    : An object of class \code{"\linkS4class{rCGH}"}
    }
    \item{Scale}{
    : logical. If \code{TRUE} (default), Log2Ratios are standardized.
    }
    \item{Cy}{
    : logical. If \code{TRUE} (default), cy3/cy5 bias is corrected using a local
    regression (loessFit). For Agilent dual-color hybridization only.
    Notice that, in case of Affymetrix files (cychp.txt or cnchp.txt), this
    argument is automatically set to \code{FALSE}, since this step is 
    managed when files are exported from ChAS or APT.
    }
    \item{GC}{
    : logical. If \code{TRUE} (default), the GC\% bias is corrected using a
    local regression (loessFit). For Agilent dual-color hybridization only.
    Notice that, in case of Affymetrix files (cychp.txt or cnchp.txt), this
    argument is automatically set to \code{FALSE}, since this step is 
    managed when files are exported from ChAS or APT.
    }
    \item{Ref}{
    : string. Which cyanine was used as the reference. Possible values are "cy3"
    (default) and "cy5". For Agilent dual-color hybridization only.
    }
    \item{suppOutliers}{
    : logical. If \code{TRUE} (default), outliers are removed using an
    Expectation-Maximization algorithm (EM). See \code{details} below.
    }
    \item{nCores}{
    : numeric. The number of cores to use in order to speed up the computation.
    When \code{NULL} (default), half of the available cores are used. If a
    value greater than the number of available cores is passed, this latter will
    be used. See \code{\link[parallel]{detectCores}}.
    }
    \item{verbose}{
    : logical. When \code{TRUE} (default), progress is printed.
    }
}

\details{
When \code{suppOutliers} is \code{TRUE} (default), the Log2Ratios are splitted
with respect to chromosomes. The main regions within each chromosome are 
identified using EM. Within each region \emph{r_i}, x[r_i] values are redifined 
according to the corresponding model parameters.
as: \deqn{x[r_i] ~ N(mu_i, sigma_i)}
Notice that this step may substantially increase the computation time.
}

\value{
An object of class \code{"\linkS4class{rCGH}"}
}

\author{
Frederic Commo
}

\seealso{
\code{\link[parallel]{detectCores}}, \code{\link[parallel]{mclapply}}
}

\examples{

filePath <- system.file("extdata", "Affy_cytoScan.cyhd.CN5.CNCHP.txt.bz2",
    package = "rCGH")
cgh <- readAffyCytoScan(filePath, sampleName = "AffyScHD")
cgh <- adjustSignal(cgh, nCores=1)
getParam(cgh)

}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
