logLik.glmmadmb <- function(object, ...)
{
  ret <- object$loglik
  class(ret) <- "logLik"
  attr(ret,"df") <- object$npar
  ## length(object$b)+length(object$S)+length(object$pz)
  ## or object$npar?
  attr(ret,"nobs") <- length(object$fitted)
  return(ret)
}
