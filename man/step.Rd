\name{step}
\alias{step}
\alias{stepAIC}
\title{
  Stubs for non-working step, stepAIC functions
}
\description{
  \code{\link{step}} (from the \code{stats} package in base R)
  and \code{\link{stepAIC}} (from the \code{MASS} package) do
  not currently work properly with \code{glmmadmb} fits.  These
  functions are an attempt to intercept such use.
}
\usage{
step(...)
stepAIC(...)
}
\arguments{
  \item{\dots}{
    stub
}
}
\value{
  None (error).
}
\note{
  \itemize{
    \item \code{\link{drop1}} \emph{does} work for single-step deletions
    (presumably \code{\link{add1}} works as well, although it has not
    been tested).
    \item One can still shoot oneself in the foot by loading the
    \code{MASS} package after \code{glmmADMB} and ignoring or suppressing
    the warnings about masking of the \code{stepAIC} function.
  }
}
\keyword{misc}
