\name{hct}
\alias{hct}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Comparison of a single armed study to a collection of study results from multiple previous clinical trials.
%%  ~~function to do ... ~~
}
\description{ When comparing a single armed study to
historical controls it is necessary to take
into account that there may be variation in the underlying
treatment/placebo effect from study to study.  If this among-study
variability is not accounted for the type one and two errors may
be inaccurate. Given a historical database of study data,
such as one might have in a meta-analysis the program
calculates the criteria for
significance for a new study that uses the database as an
historical control and calculates the power of such as study
as a function of sample size and difference to be detected.}
\usage{
hct(data, estimate, standardError, N, iter = 2000, rseed = NA,
silent=TRUE,constantStderr=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A data frame of historical data one study per row.}
  \item{estimate}{The name or column number of the variable in
  \code{data} containing the estimated outcome.}
  \item{standardError}{The name or column number of the variable in \code{data}
  containing the standard error of the estimated outcome.}
  \item{N}{ The name or column number of the variable in \code{data}
  containing the sample size of the study.}
  \item{iter}{
The number of interations to use in the MCMC to
calculate the posterior distribution of the among-study
variation and mean outcome measure.}
  \item{rseed}{Seed for random number generator}
  \item{silent}{Suppresses STAN output to the console}
  \item{constantStderr}{If \code{TRUE} it assumes that that standard
  deviation for each study is known exactly this is
  appropriate for larger studies.
  If \code{FALSE} it assumes that they are proportional to a
  chi-square distribution with \code{N-1} degrees of freedom and
  uses a hierarchical model for the patient-level variance.
}
}
\value{
A \code{hct} object which is a list of four elements.

\item{criteria}{A function with signature \code{(p,se,df==NULL)} to calculate the cut-off value
for statistical significance at a one sided p-value of \code{p},  with standard
error of the estimate equal to \code{se}. When \code{df=NULL} a normal test is used otherwise a t-test with \code{df} degrees of freedom.}

\item{power}{A function to calculate the power of a study with signature \code{(t,delta,se,df=NULL)},
where \code{t} is the cuttoff value, \code{se} is the standard error of the estimate
and \code{delta}
is the treatment effect}
\item{effective.SD}{Which is the effective standard deviation of the outcome measure.
It is calculated as
\code{sqrt(sum(data[,standardError]^2*data[,N]*(data[,N]-1))/(sum(data[,N])-1)))},
which is what it would be if the parameter estimates were sample means. When \code{df=NULL} a normal test is used otherwise a t-test with \code{df} degrees of freedom.}
\item{fit}{An object of class \code{stanfit} with the fit of the data}

A generic summary function prints out the value of \code{effective.SD} and uses the data
frame \code{summary} function for
\code{data.frame(extract(fit,c("mu"","sig"","y_pred")))}. The generic \code{print}
function prints this summary.
}
\references{
\href{http://hedwig.mgh.harvard.edu/biostatistics/node/98}{Design and
analysis of a clinical trial using previous trials as historical control}
}
\author{
David A. Schoenfeld
}
\examples{
\donttest{
als=data.frame(estimate=c(3.5,2.6,2.3),SE=c(.4,.3,.6),N=c(100,150,76))
ts=hct(als,'estimate','SE','N')
print(ts)
us=ts$criteria(0.025,.3)
ts$power(us,5,.4)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models }% use one of  RShowDoc("KEYWORDS")
