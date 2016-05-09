## other accessor methods (some trivial)
coef.glmmadmb <- function(object, ...) {
  object$b
}

## for lme4/nlme compatibility
fixef.glmmadmb <- function(object, ...) {
  object$b
}

## need to make sure that this plays nicely with
##  lme4 (S4 methods).  Not sure how.

## ranef <- function(object, ...) {
##   UseMethod("ranef")
##}

## setGeneric("ranef", function(object, ...) {
##     standardGeneric("ranef")
## })

## setMethod("ranef","glmmadmb",
##           function(object, sd=FALSE, ...) {
##    if(sd) return(object$sd_U)
##    mapply(sweep,object$U,lapply(object$S,function(z)sqrt(diag)),
##           MoreArgs=list(MARGIN=2,FUN="*"),SIMPLIFY=FALSE)
##  })

ranef.glmmadmb <- function(object, sd=FALSE, scale=TRUE, ...) {
  sdvals <- lapply(object$S,function(z)sqrt(diag(z)))
  X <- if (sd) object$sd_U else object$U
  if (scale) X <- mapply(sweep,X,sdvals,
                         MoreArgs=list(MARGIN=2,FUN="*"),SIMPLIFY=FALSE)
  return(X)
}

residuals.glmmadmb <- function(object, type=c("pearson", "response"), ...) {
  type <- match.arg(type)
  if (type=="response") {
    object$residuals
  } else {
    object$residuals/object$sd.est
  }
}

fitted.glmmadmb <- function(object, ...) {
  object$fitted
}

## generic is defined in R2admb
## stdEr <- function(x, ...) {
##   UseMethod("stdEr")
## }

stdEr.glmmadmb <- function(object, ...) {
  object$stdbeta
}

vcov.glmmadmb <- function(object, ...) {
  outer(object$stdbeta,object$stdbeta)*object$corMat
}

nobs.glmmadmb <- function(object,...) {
  length(object$fitted)
}

## VarCorr <- function(x,...) {
##   UseMethod("VarCorr")
##}

## big difficulty here with nlme (S3 method, arguments x, sigma=1, rdig=3)
##  and lme4 (S4 methods, arguments x, ...)
VarCorr.glmmadmb <- function(x,sigma=1,rdig=3) {
  if (!missing(sigma) || !missing(rdig)) warning("'sigma' and 'rdig' arguments are present for compatibility only: ignored")
  vc <- x$S
  class(vc) <- "VarCorr"
  vc
}

print.VarCorr <- function(x, digits=4, ...) {
  for (i in seq_along(x)) {
    cat("Group=",names(x)[i],"\n",sep="")
    vc <- x[[i]]
    v <- diag(vc)
    vmat <- cbind(Variance=v,StdDev=sqrt(v))
    if (nrow(vc)==1 || all(vc[lower.tri(vc)]==0)) {
      print(vmat,digits=digits)
    } else {
      cmat <- matrix("",nrow=nrow(vmat),ncol=nrow(vmat)-1)
      cc <- cov2cor(vc)
      cmat[lower.tri(cmat)] <- format(cc[lower.tri(cc)],digits=digits)
      colnames(cmat) <- c("Corr",rep("",ncol(cmat)-1))
      cmat[1,] <- abbreviate(rownames(vc)[-nrow(vc)],digits+2)
      vmat <- format(vmat,digits=digits)
      print(cbind(vmat,cmat),quote=FALSE)
    }
  }
}

VarCorr.summary.glmmadmb <- VarCorr.glmmadmb

## want to make this work when lme4 is loaded, too ... needs S4 method
setOldClass("glmmadmb")
setOldClass("summary.glmmadmb")
setMethod("VarCorr", signature(x="glmmadmb"), VarCorr.glmmadmb)
setMethod("VarCorr", signature(x="summary.glmmadmb"), VarCorr.glmmadmb)
## FIXME:
##   needed:
##    update (for general convenience & to make drop1 work)
##    terms, extractAIC  (to make drop1 work)
##      for terms, do we want to save model frame? save_frame
##          (or save.frame or saveFrame) ?


model.frame.glmmadmb <- function(formula,...) {
    formula$frame
}

df.residual.glmmadmb <- function(object,...) {
    nparams <- nrow(object$frame)-object$npar
    ## FIXME: is npar correct ???
}

## for drop1 etc.
extractAIC.glmmadmb <- function(fit,scale,k=2,...) {
    if (!missing(scale) && scale!=1) warning("ignored explicit specification of scale")
    L <- logLik(fit)
    edf <- attr(L,"df")
    c(edf=edf,AIC=-2*L+k*edf)
}

step <- stepAIC <- function(...) {
    stop("functions step and MASS::stepAIC are **not** currently compatible with glmmADMB.  Sorry.")
}
