\name{predict.glmmadmb}
\alias{predict.glmmadmb}
\title{
  Prediction for glmmADMB models
}
\description{
  Generates predicted values, standard errors and confidence intervals
  for glmmADMB fits
}
\usage{
\S3method{predict}{glmmadmb}(object,
 newdata = NULL, type = c("link", "response"),
  se.fit = FALSE, interval = c("none", "confidence"),
  random = ~0, level = 0.95, \dots)
}
\arguments{
  \item{object}{
    A fitted model from \code{\link{glmmadmb}}
  }
  \item{newdata}{
    An optional data frame in which to look for variables with
    which to predict.  If omitted, the fitted values are used.
  }
  \item{type}{
    Whether to return predictions on the
    scale of the linear predictor (\code{"link"}) or
    the scale of the data (\code{"response"}).
  }
  \item{se.fit}{
    A switch indicating if standard errors are required.
  }
  \item{interval}{
    What type of confidence intervals to return
  }
  \item{random}{
    Specification of random effects to use (STUB)
  }
  \item{level}{
    Tolerance/confidence level
  }
  \item{\dots}{
    further arguments passed to or from other methods.
  }
}
\details{
  Produces predicted values from the model, based on evaluating
  the model with 'newdata' (if specified).  If the logical
  'se.fit' is 'TRUE', standard errors are calculated (this is
  only available on the link scale).  It 
}
\value{
  \code{predict.lm} produces a vector of predictions or a data frame of
  predictions and bounds with column names \code{fit}, \code{lwr}, and
  \code{upr} if \code{interval} is set.  If \code{se.fit} is
  \code{TRUE}, a list with the components \code{fit}
  (vector or matrix as above) and \code{se.fit}
  (standard error of predicted means) is returned.
  }
\author{
  Ben Bolker (partially copied from \code{\link{predict.lm}})
}
\note{
  \itemize{
    \item Prediction is currently available only at the population level (i.e.
    with all random effects terms set to zero).
    \item Confidence intervals and standard errors neglect uncertainty
    due to random effects, as well as the uncertainty of random
    effects. For more accurate confidence intervals, use MCMC-sampling
    methods.
  }
}
\examples{
getdata("OwlModel")
pp <- predict(OwlModel,interval="confidence",type="response")
pframe <- data.frame(obs=Owls$NegPerChick,pp)
with(pframe,plot(obs,fit,ylim=c(0,8)))
with(pframe,arrows(obs,lwr,obs,upr,length=0))
abline(a=0,b=1,col=2,lty=2) ## 1-to-1 line
}
\keyword{misc}

