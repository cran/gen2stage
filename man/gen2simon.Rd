\name{gen2simon}
\title{Generalized Simon's 2-stage phase II design}
\alias{gen2simon}
\alias{print.gen2simon}
\alias{plot.gen2simon}
\keyword{design}
\description{
  Calculates generalized optimal and minimax 2-stage phase II designs based
  on the R function ph2simon.
}
\usage{
gen2simon(pu, pa, ep1, ep2, nmax=100)
\method{print}{gen2simon}(x, \dots)
\method{plot}{gen2simon}(x, \dots)
}
\arguments{
  \item{pu}{unacceptable response/toxicity rate}
  \item{pa}{response/toxicity rate that is desirable}
  \item{ep1}{threshold for the probability of declaring drug/treatment desirable
    under p0}
  \item{ep2}{threshold for the probability of rejecting the drug/treatment under
    p1}
  \item{nmax}{maximum total sample size (default 100; can be at most
    500)}
  \item{x}{object returned by gen2simon}
  \item{...}{arguments to be passed onto plot and print commands called
    within}
}
\value{
  gen2simon returns a list with pu, pa, alpha, beta and nmax as above
  and:
  \item{out}{matrix of best 2 stage designs for each value of total
    sample size n.  The 6 columns are: r1, n1, r, n, EN(p0), PET(p0), alpha, beta}

  The "print" method formats and returns the minimax and optimal
  designs.  The "plot" plots the expected sample size agains the maximum
  sample size as in Jung et al., 2001
}

\seealso{
  \code{\link{oc.gentwostage.bdry}}
}

\examples{
  # Two-stage safety design with pu (p0) = 0.33 vs. pa (p1) = 0.20
  gen2simon(0.33, 0.20, 0.05, 0.20)
  gen2simon(0.33, 0.20, 0.05, 0.10, nmax=150)

  # Two-stage efficacy design with pu (p0) = 0.67 vs. pa (p1) = 0.80
  gen2simon(0.67, 0.80, 0.05, 0.20)
  gen2simon(0.67, 0.80, 0.05, 0.10, nmax=150)

  # save, print and plot
  safety2 = gen2simon(0.33, 0.20, 0.05, 0.20)
  print(safety2)
  plot(safety2)
}

\references{
  Kim S and Wong WK. Phase II Two-Stage Single-Arm Clinical Trials for Testing Toxicity Levels. \emph{Commun Stat Appl Methods. 2019 Mar;26(2):163-173.} \url{https://www.ncbi.nlm.nih.gov/pubmed/31106162}.

  Jung SH, Carey M and Kim KM. (2001). Graphical Search for Two-Stage
  Designs for Phase II Clinical Trials. \emph{Controlled Clinical
    Trials}  22, 367-372.

  Simon R. (1989).  Optimal Two-Stage Designs for Phase II Clinical
  Trials. \emph{Controlled Clinical Trials} 10, 1-10.
}

\keyword{design}
