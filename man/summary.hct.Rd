\name{summary.hct}
\alias{summary.hct}
\alias{summary}
\title{Summary function for an hct object.
}
\description{
Summary of a hct object.
}
\usage{
summary(object,...)}

\arguments{
  \item{object}{hct object.}
  \item{...}{arguments passed to \code{summary.data.frame}}}
\value{A list with entries \code{effective.SD} and \code{prior.distribution}, the
first is the effective.SD and the second is the \code{summary.data.frame} applied to the
\code{prior.distribution} element of the hct object. Note that the \code{prior.distribution} can be considered a posterior.distribution as it is the posterior distribution of the parameters with the input data base, however it is used as a prior distribution for the purpose of analyzing the historical trial.}
\author{David A. Schoenfeld}
\seealso{\code{\link{hct}}}
\keyword{Summary }
