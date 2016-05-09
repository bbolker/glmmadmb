## modeled after summary.glm, print.summary.glm 
summary.glmmadmb <- function(object, ...) 
{
  ## print.glmmadmb(object, ...)
  ## calculate coef table

  ## for now, assume dispersion KNOWN
  ##  glm.nb inherits from glm, so summary.glm is used
  ##    don't allow for uncertainty from estimating theta??
  ## est.disp <- object$family %in% c("binom","poisson")
  est.disp <- FALSE
  
  coef.p <- object$b
  s.err <- object$stdbeta
  tvalue <- coef.p/s.err

  dn <- c("Estimate", "Std. Error")
  if(!est.disp) { # known dispersion
    pvalue <- 2*pnorm(-abs(tvalue))
    coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
    dimnames(coef.table) <- list(names(coef.p),
                                 c(dn, "z value","Pr(>|z|)"))
  }
  ans <- c(object,
           list(aic=AIC(object),coefficients=coef.table))
  class(ans) <- "summary.glmmadmb"
  ## modeled after summary.glm
  ans
}

print.summary.glmmadmb <- function(x, digits=max(3, getOption("digits") - 4),
                                   symbolic.cor=x$symbolic.cor,
                                   signif.stars=getOption("show.signif.stars"),
                                   ...)
{
    cat("\nCall:\n",
        paste(deparse(x$call), sep="\n", collapse = "\n"), "\n\n", sep="")
    ## print deviance residuals?
    ## cat("Deviance Residuals: \n")
    ## if(x$df.residual > 5) {
    ##     x$deviance.resid <- quantile(x$deviance.resid,na.rm=TRUE)
    ##     names(x$deviance.resid) <- c("Min", "1Q", "Median", "3Q", "Max")
    ## }
    ## xx <- zapsmall(x$deviance.resid, digits + 1)
    ## print.default(xx, digits=digits, na.print = "", print.gap = 2)
    cat("AIC:",round(c(x$aic),1),"\n")
    cat("\nCoefficients:\n")
    coefs <- x$coefficients
    printCoefmat(coefs, digits=digits, signif.stars=signif.stars,
                 na.print="NA", ...)
    cat("\n","Number of observations: total=",x$n,sep="")
    if (!is.null(x$S)) {
        ## has random effects
        cat(", ")
        cat(paste(names(x$q),x$q,sep="=",collapse=", "),"\n")
        cat("Random effect variance(s):\n")
        print(VarCorr(x))
        ## FIXME: prettier? standard errors?
    }
    cat("\n")
    if (!is.null(x$alpha)) {
      label <- switch(x$family,truncnbinom=,nbinom="Negative binomial dispersion parameter",
                      gamma="Gamma shape parameter",
                      beta="Beta dispersion parameter",
                      betabinom="Beta-binomial dispersion parameter",
                      gaussian="Residual variance",
                      logistic="Scale parameter",
                      "Unknown dispersion parameter")
      cat(label,": ",x$alpha," (std. err.: ",x$sd_alpha,")\n",
          sep="")
    }

    if (!is.null(x$pz)) {
      cat("Zero-inflation:",x$pz," (std. err.: ",x$sd_pz,")\n")
    }

    cat("\nLog-likelihood:",x$loglik,"\n")
    ## offset
    ## cat("\nOffset:\
}
