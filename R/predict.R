nobs.glmmadmb <- function(object,...) {
  nrow(object$frame)
}

predict.glmmadmb <- function(object, newdata=NULL,
                             type=c("link","response"),
                             se.fit=FALSE,
                             interval=c("none","confidence"),
                             random=~0,
                             level=0.95, ...) {
  if (se.fit && type=="response") {
    warning("se.fit && type='response': setting se.fit to NA")
  }
  type <- match.arg(type)
  interval <- match.arg(interval)
  ## Construct model matrix, nobs x np
  if (missing(newdata) || is.null(newdata)) {
    newdata <- object$frame
    X <- model.matrix(object, data=newdata)
    offset <- object$offset
  } else {
    form <- as.formula(as.character(object$fixed)[-2])
    X <- model.matrix(form, data=newdata)
    tt <- object$terms
    ## handle offset
    offset <- rep(0, nrow(X))
    off.num <- attr(tt, "offset")
    if (!is.null(off.num))  {
      for (i in off.num) {
        cur.offset <- eval(attr(tt, "variables")[[i+1]], newdata)
        offset <- offset + cur.offset
      }
    }
  }
  beta <- as.vector(object$b)
  phat <- X %*% beta
  ## interpret random-effects formula
  if (!identical(random,~0)) {
    stop("random-effects prediction not yet implemented")
  }  
  if (!is.null(offset)) phat <- c(phat + offset)
  if (se.fit || interval!="none") {
    stderr <- c(sqrt(diag(X %*% vcov(object) %*% t(X))))
  }
  if (interval=="confidence") {
    qq <- qnorm((1+level)/2)
    phat <- cbind(phat,phat-qq*stderr,phat+qq*stderr)
    colnames(phat) <- c("fit","lwr","upr")
  }
  if (type=="response") {
    phat <- object$ilinkfun(phat)
    stderr <- NA
  }
  if (interval=="confidence") phat <- as.data.frame(phat)
  if (se.fit) list(fit=phat,se.fit=stderr) else phat
}
