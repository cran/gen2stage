\name{gen2stage-package}
\alias{gen2stage-package}
\alias{gen2stage}
\docType{package}
\title{
Generalized Two-Stage Designs for Phase II Single-Arm Studies
}
\description{
  The R package \code{gen2stage} can generate single-stage and two-stage designs
  for phase II single-arm efficacy or safety studies.
}
\details{
\tabular{ll}{
Package: \tab gen2stage\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2017-10-05\cr
License: \tab GPL-2\cr
}
}
\author{
Seongho Kim <biostatistician.kim@gmail.com>
}
\references{
  Kim S and Wong WK. Phase II Two-Stage Single-Arm Clinical Trials for Testing Toxicity Levels. \emph{Commun Stat Appl Methods. 2019 Mar;26(2):163-173.} \url{https://www.ncbi.nlm.nih.gov/pubmed/31106162}.
}
\keyword{design}
\examples{
  # Single-stage safety design with pu (p0) = 0.33 vs. pa (p1) = 0.20
  gen2single(0.33, 0.20, 0.05, 0.20)

  # Single-stage efficacy design with pu (p0) = 0.67 vs. pa (p1) = 0.80
  gen2single(0.67, 0.80, 0.05, 0.20)

  # save and print
  safety1 = gen2single(0.33, 0.20, 0.05, 0.20)
  print(safety1)

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
