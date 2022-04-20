\name{ALShctObjects}
\alias{ALShctObjects}
\alias{hct220ALSFRS}
\alias{hct220Delta}
\docType{data}
\title{Objects Created by hct for the ALS of Clinical Trial Results}

\description{These are objects of class \code{hct} for the example an ALS clinical trials}
\format{
   Two \code{hct} objects  \code{hct220ALSFRS,hct220Delta.} These are fit
   from the collection of ALS clinical trials without and with covariate
   predictions, respectively. They allow the use of the
   functions \code{criteria} and \code{power} functions to
   calculate statistics of interest, such as p-values for an
   historically controlled ALS trial}

\examples{
data(ALShctObjects)
#calculate the p-value of an ALS clinical trial where he
#slope of the ALSRS was -2 with standard Error 0.2
hct220ALSFRS$power(-.2,0,.2)
#Calculate the criteria for a p-value of 0.10 one sided to us in #a pilot study
hct220ALSFRS$criteria(0.1,.2)
}
\keyword{datasets}
