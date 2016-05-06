\name{EMnormalize}
\alias{EMnormalize}
\alias{EMnormalize-methods}
\alias{EMnormalize,rCGH-method}

\title{Genomic Profile Centralization}

\description{
This function analyses the Log2Ratios as a mixture of several gaussian
populations, using an Expectation-Maximization algorithm (EM).\cr
The \code{peakThresh} argument specifies what proportion of the main density
peak is allowed for choosing a neutral 2-copies population. The mean of the 
chosen population is used for centralizing the profile.\cr
See \code{\link[mclust]{Mclust}}.
}

\usage{
\S4method{EMnormalize}{rCGH}(object, cut = c(0.01, 0.99), G = 2:6, useN = 25e3,
peakThresh = 0.5, ksmooth = NA, mergeVal = 0.1, Title = NA, verbose=TRUE)
}

\arguments{
    \item{object}{
    : An object of class \code{"\linkS4class{rCGH}"}
    }
    \item{cut}{
    : numeric. A vector of 2 alpha values (between 0 and 1). Log2Ratios outside
    the corresponding quantiles will be excluded for the gaussian mixture
    estimation. Default quantiles are \eqn{q_{0.01}} and \eqn{q_{0.99}}.
    }
    \item{G}{
    : numeric. The number of groups to test during the gaussian mixture
    estimation. Default is from 2 to 6.
    }
    \item{useN}{
    : numeric. The number of probes to use for estimating the mixture 
    parameters. Default is 25e3.
    }
    \item{peakThresh}{
    : numeric. The proportion of the highest peak to consider as a peak
    selection threshold. Default is 0.5.
    }
    \item{ksmooth}{
    : numeric. A smoothing value applied to Log2Ratios before modeling the
    gaussian mixture. When \code{NA} (default) \code{ksmooth} is estimated from
    the median absolute deviation of the Log2Ratios.
    }
    \item{mergeVal}{
    : numeric. Populations with means closer than \code{mergeVal} will be
    pooled together, default is 0.1. 
    Set \code{mergeVal} to zero to not pool closed sub-populations.
    }
    \item{Title}{
    : character string. A title for the density plot. If \code{NA} (default),
    the sample name (when exists in object info) will be used.
    }
    \item{verbose}{
    : logical. When \code{TRUE} (default), progress is printed.
    }
}
\details{
Depending on \code{peakThresh}, the mean of the highest density, or a lower
value, will be chosen for centering the Log2Ratios before the segmentation.\cr
When a \code{peakThresh} value is specified, heights of density peaks are
compared: the lowest peak mean among the peaks respecting the criteria: 
\code{peakHeight > max(peaks)*peakThresh}, is chosen for centralizing the data.
See \code{References}
}

\value{
An object of same class as the input.
}

\references{
\href{http://www.ncbi.nlm.nih.gov/pubmed/25538175}{Commo et al. Impact of
centralization on aCGH-based genomic profiles for precision medicine in
oncology. Ann Oncol. 2014}
}

\author{
Frederic Commo
}

\seealso{
\code{\link{plotDensity}}, \code{\link[mclust:mclust-package]{mclust}}
}

\examples{

filePath <- system.file("extdata", "Affy_cytoScan.cyhd.CN5.CNCHP.txt.bz2",
    package = "rCGH")
cgh <- readAffyCytoScan(filePath, sampleName = "AffyScHD")
cgh <- adjustSignal(cgh, nCores=1)
cgh <- EMnormalize(cgh)
getParam(cgh)

}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }

