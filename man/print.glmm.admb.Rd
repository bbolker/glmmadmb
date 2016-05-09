\name{print.glmmadmb}
\alias{print.glmmadmb}
\alias{coef.glmmadmb}
\alias{fixef}
\alias{fixef.glmmadmb}
\alias{residuals.glmmadmb}
\alias{fitted.glmmadmb}
\alias{vcov.glmmadmb}
\alias{VarCorr}
\alias{VarCorr.glmmadmb}
\alias{VarCorr.summary.glmmadmb}
\alias{VarCorr,glmmadmb-method}
\alias{VarCorr,summary.glmmadmb-method}
\alias{stdEr}
\alias{stdEr.glmmadmb}
\title{Basic accessor methods for \code{glmmadmb} objects}
\description{Basic accessor methods for \code{glmmadmb} objects:
  print, retrieve coefficients, etc.}
\usage{
stdEr(object,\ldots)
fixef(object,\ldots)
\method{print}{glmmadmb}(x,\dots)
\method{coef}{glmmadmb}(object,\dots)
\method{fixef}{glmmadmb}(object,\dots)
\method{residuals}{glmmadmb}(object,type=c("pearson", "response"),\dots)
\method{fitted}{glmmadmb}(object,\dots)
\method{vcov}{glmmadmb}(object,\dots)
\method{stdEr}{glmmadmb}(object,\dots)
\method{VarCorr}{glmmadmb}(x,sigma=1,rdig=3)
\method{VarCorr}{summary.glmmadmb}(x,sigma=1,rdig=3)
}
\arguments{
  \item{x}{a fitted \code{glmmadmb} model}
  \item{object}{a fitted \code{glmmadmb} model}
  \item{type}{type of residuals to retrieve}
  \item{sigma}{for compatibility with \code{\link[nlme]{VarCorr}}:
    ignored}
  \item{rdig}{for compatibility with \code{\link[nlme]{VarCorr}}: ignored}
  \item{\dots}{further arguments. If \code{sd_S_print=TRUE}
    for \code{print.glmmadmb},
    print the covariance matrix of the random effects vector and corresponding
    standard deviations}
}
\details{
\itemize{
\item \code{coef} extracts (fixed-effect) coefficients from a
\code{glmmadmb} object, or the coefficient table (estimate, standard
error, Z-score, p-value) from a \code{summary.glmmadmb} object;
\code{fixef} is a synonym, for \code{nlme}/\code{lme4} compatibility
\item \code{vcov} extracts the approximate variance-covariance matrix
 of the fixed-effect parameters, based on the information matrix
\item \code{stdEr} extracts the approximate standard errors of the
fixed-effect parameters, as above
\item \code{VarCorr} extracts the random effects variance-covariance matrices
}
}

\seealso{
  \code{\link{glmmadmb}}, \code{\link{anova.glmmadmb}},
  \code{\link{logLik.glmmadmb}}
}
\keyword{models}
