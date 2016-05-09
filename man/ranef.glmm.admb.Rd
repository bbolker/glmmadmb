\name{ranef.glmmadmb}
\alias{ranef}
\alias{ranef.glmmadmb}
\title{Extract Random Effects}
\description{
  Extract random effects from fitted \code{glmmadmb} models, and
  (optionally) their standard deviations.
}
\usage{
\method{ranef}{glmmadmb}(object, sd=FALSE, scale=TRUE, \dots)
}
\arguments{
  \item{object}{a fitted \code{glmmadmb} model.}
  \item{sd}{whether the standard deviations of the random effects should
    be returned, rather than point estimates.}
  \item{scale}{(logical) scale random effects by the relevant
    random-effects standard deviations (i.e., the resulting effects
    should have variance approximately equal to the estimated RE
    variance if \code{scale=TRUE}, or approximately equal to 1 if
    \code{scale=FALSE})}
  \item{\dots}{further arguments (currently ignored).}
}
\details{
  The point estimates of the random effects are \sQuote{empirical Bayes}
  estimates in the sense that the hyperparameters (fixed effects,
  variance components, etc.) are fixed at their MLE's when the joint
  likelihood is maximized with respect to the random effects. However,
  uncertainty about the hyperparameters is taken into account when the
  standard deviations of the random effects are calculated.  The
  standard deviations of the random effects are quadratic
  approximations, so should be taken with a grain of salt for small
  sample sizes.
}
\value{
  A matrix of estimated random effects (or standard deviations).
}
\seealso{
  \code{\link{glmmadmb}}.

  \code{\link{glmmADMB-package}}.
}
\keyword{models}
