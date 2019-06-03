\name{gen2single}
\title{Generalized exact single stage phase II design}
\alias{gen2single}
\alias{print.gen2single}
\keyword{design}
\description{
  Calculates the generalized exact one stage phase II design based
  on the R function ph2single.
}
\usage{
gen2single(pu,pa,ep1,ep2,nsoln=5)
\method{print}{gen2single}(x, \dots)
}
\arguments{
  \item{pu}{unacceptable response/toxicity rate}
  \item{pa}{response/toxicity rate that is desirable}
  \item{ep1}{threshold for the probability of declaring drug/treatment desirable
    under p0}
  \item{ep2}{threshold for the probability of rejecting the drug/treatment under
    p1}
  \item{nsoln}{number of designs with given alpha and beta}
  \item{x}{object returned by gen2single}
  \item{...}{arguments to be passed onto print command called
    within}}
\value{
  gen2single returns the optimal design with pu, pa, alpha, and beta as above
  and:
  \item{out}{matrix of the single-stage designs up to nsoln.  The 4 columns are: r, n, alpha (type I error),
    beta (type II erro)}

  The "print" method formats and returns the optimal
  design.
}
\seealso{
  \code{\link{gen2simon}}
}
\examples{
  # Single-stage safety design with pu (p0) = 0.33 vs. pa (p1) = 0.20
  gen2single(0.33, 0.20, 0.05, 0.20)

  # Single-stage efficacy design with pu (p0) = 0.67 vs. pa (p1) = 0.80
  gen2single(0.67, 0.80, 0.05, 0.20)

  # save and print
  safety1 = gen2single(0.33, 0.20, 0.05, 0.20)
  print(safety1)
}

\references{
  Kim S and Wong WK. Phase II Two-Stage Single-Arm Clinical Trials for Testing Toxicity Levels. \emph{Commun Stat Appl Methods. 2019 Mar;26(2):163-173.} \url{https://www.ncbi.nlm.nih.gov/pubmed/31106162}.
}

\keyword{design}
