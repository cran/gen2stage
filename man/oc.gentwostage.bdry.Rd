\name{oc.gentwostage.bdry}
\title{Two-stage boundary operating characteristics}
\alias{oc.gentwostage.bdry}
\keyword{design}
\description{
  Calculates the operating characteristics of a two-stage boundary based
  on the R function oc.twostage.bdry.
}
\usage{
oc.gentwostage.bdry(pu, pa, r1, n1, r, n)
}
\arguments{
  \item{pu}{unacceptable response rate}
  \item{pa}{response rate that is desirable}
  \item{r1}{first stage threshold to declare treatment undesirable}
  \item{n1}{first stage sample size}
  \item{r}{overall threshold to declare treatment undesirable}
  \item{n}{total sample size}
}
\value{
  oc.gentwostage.bdry returns the type I and II error rates as well as the
  probability of early temination and expected sample size under pu for
  a specific boundary.
}
\seealso{
  \code{\link{gen2simon}}
}
\examples{
  # Optimal two-stage safety design with pu (p0) = 0.33 vs. pa (p1) = 0.20
  oc.gentwostage.bdry(0.33, 0.20, 8, 26, 22, 85)

  # Optimal two-stage efficacy design with pu (p0) = 0.67 vs. pa (p1) = 0.80
  oc.gentwostage.bdry(0.67, 0.80, 18, 26, 63, 85)
}

\references{
  Kim S and Wong WK. Phase II Two-Stage Single-Arm Clinical Trials for Testing Toxicity Levels. \emph{Commun Stat Appl Methods. 2019 Mar;26(2):163-173.} \url{https://www.ncbi.nlm.nih.gov/pubmed/31106162}.
}

\keyword{design}
