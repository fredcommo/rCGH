\name{agilentDB}
\alias{agilentDB}
\docType{data}
\title{aCGH Agilent Probes GC Fraction}
\description{
A dataset containing the Agilent probe Ids, and their GC content.\cr
These data allow \code{\link{rCGH}} to support Agilent aCGH arrays, from 44K to
400K. See the \code{source} section.\cr
This information is used to correct the GC\% bias. For Agilent data only.
}
\usage{agilentDB}
\value{A dataset}
\format{
A data frame with 411056 rows and 2 columns:
    \itemize{
        \item ProbeID: Official Agilent probe ids.
        \item GC: The GC content in each probe sequence, expressed as the
        GC fraction.
    }
}
\author{Frederic Commo}
\source{
These data derive from the official Agilent Sequence List,
SurePrint_G3_Human_CGH_Microarray 2x400K_021850_D_SequenceList_20111015.txt,
freely available at: 
\href{https://earray.chem.agilent.com/suredesign/index.htm}{Agilent SureDesign}
\cr
Access date: 3-2-2015\cr
Notice that a User ID and Password are required to sign in.
}
\note{
    Probe sequences are not used and have been removed after computing the GC
    fractions as GC = sum(G or C)/length(sequence), for each sequence in file.
}
\keyword{datasets}
