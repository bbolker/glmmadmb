anova.glmmadmb <- function(object, ...)
{
  objects <- list(object, ...)

  if(length(objects) < 2)
    stop("Two or more model fits required.")
  ## FIXME: no longer effective (random is now NULL).
  ##  check equality of random terms??
  ## if(length(unique(paste(lapply(objects,function(x) x$random)))) > 1)
  ## stop("Random effects are not identical")

  npar <- as.numeric(lapply(objects, function(x) x$npar))
  if (any(diff(npar)<0)) {
    warning("rearranging models in order of increasing complexity")
    npar_order <- order(npar)
    objects <- objects[npar_order]
    npar <- npar[npar_order]
  }
  logLik <- as.numeric(lapply(objects, function(x) x$loglik))

  if (any(diff(logLik)<0)) warning("something's wrong: models should be nested, increasing complexity should imply increasing log-likelihood")

  df <- c(NA, diff(npar))

  n2logQ <- 2 * c(NA,diff(logLik))
  P.value <- c(NA, 1-pchisq(n2logQ[-1],df[-1]))
  table <- data.frame(npar, logLik, df, n2logQ, P.value)
  variables <- lapply(objects, function(x) x$fixed)

  dimnames(table) <- list(1:length(objects), c("NoPar","LogLik","Df","Deviance","Pr(>Chi)"))
  title <- "Analysis of Deviance Table\n"
  topnote <- paste("Model ", format(1:length(objects)), ": ", variables, sep="", collapse="\n")

  output <- structure(table, heading=c(title,topnote), class=c("anova","data.frame"))

  return(output)
}
