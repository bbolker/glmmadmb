\name{anova.glmmadmb}
\alias{anova.glmmadmb}
\title{Model Comparison Table}
\description{Compute a model comparison (likelihood ratio test) table for a series of fitted models.}
\usage{
\method{anova}{glmmadmb}(object, \dots)
}
\arguments{
  \item{object}{a fitted \code{glmmadmb} model.}
  \item{\dots}{More complex fitted models. These may differ in fixed
    or random effects, but note that the likelihood ratio test is
    generally conservative for testing among models with differing
    random terms, because of boundary effects.}
}
\details{
  Basic table for testing the fixed effects, using
  likelihood ratio tests (may be extended
  later).
}
\note{The models must be nested and arranged in order
  of increasing complexity.}
\value{
  An object of class \code{"anova"}, inheriting from class
  \code{"data.frame"}.
}
\seealso{
  \code{\link{glmmadmb}},
  \code{\link{glmmADMB-package}}.
}
\examples{
  data(epil2)
  epil2$subject <- factor(epil2$subject)
\dontrun{
  ## currently (June 2013) having trouble with Windows checks on r-forge, CRAN
  f1 <- glmmadmb(y~Base*trt+Age+Visit+(Visit|subject),
                data=epil2, family="nbinom")
  f2 <- glmmadmb(y~Base*trt+Age+(Visit|subject),
                data=epil2, family="nbinom")
  anova(f2,f1)
}
}
\keyword{models}
