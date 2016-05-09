\name{glmmADMB-package}
\alias{glmmADMB-package}
\alias{glmmADMB}
\docType{package}
\title{Generalized Linear Mixed Models Using AD Model Builder}
\description{
  Fits (extended)
  generalized linear mixed-effects models to data using a variety of
  distributions and link functions, including zero-inflated models.
}
\details{
  \emph{Main function:}
  \tabular{ll}{
    \code{\link{glmmadmb}} \tab fit a model
  }
  \emph{Methods:}
  \tabular{ll}{
    \code{\link{coef.glmmadmb}}   \tab extract fixed effect coefficients \cr
    \code{\link{anova.glmmadmb}}   \tab likelihood ratio tests between
    models \cr
    \code{\link{logLik.glmmadmb}}  \tab log likelihood\cr
    \code{\link{plot.glmmadmb}}    \tab plot\cr
    \code{\link{print.glmmadmb}}   \tab print\cr
    \code{\link{ranef.glmmadmb}}   \tab random effect coefficients \cr
    \code{\link{summary.glmmadmb}} \tab summary \cr
    \code{\link{fitted.glmmadmb}} \tab fitted values \cr
    \code{\link{predict.glmmadmb}} \tab predicted values \cr
    \code{\link{residuals.glmmadmb}} \tab Pearson residuals  \cr
    \code{\link{stdEr.glmmadmb}} \tab standard errors of coefficients \cr
    \code{\link{vcov.glmmadmb}} \tab variance-covariance matrix of
    fixed effects \cr
    \code{\link{VarCorr.glmmadmb}} \tab random effects variances and covariances
  }
  \emph{Example data:}
  \tabular{ll}{
    \code{\link{epil2}} \tab seizure counts
  }
}
\note{
  If you find that the \pkg{glmmADMB} package does not do exactly what
  you want, you may want to modify the underlying ADMB \file{.tpl} files
  being called by the \code{glmmadmb} function, which can be found via
  \code{system.file("tpl","glmmadmb.tpl",package="glmmADMB")}.

}
\author{
  Hans Skaug, Dave Fournier, Anders Nielsen, Arni Magnusson, and Ben
  Bolker.
}
\seealso{\url{http://admb-project.org}}
\keyword{models}
