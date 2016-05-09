library(glmmADMB)

## copied from emdbook package:
rzinbinom <- function (n, mu, size, zprob)  {
    ifelse(runif(n) < zprob, 0, rnbinom(n, mu = mu, size = size))
}

dzinbinom <- function (x, mu, size, zprob, log = FALSE)  {
    logv <- log(1 - zprob) + dnbinom(x, mu = mu, size = size, 
        log = TRUE)
    logv <- ifelse(x == 0, log(zprob + exp(logv)), logv)
    if (log) logv else exp(logv)
  }

set.seed(1001)
y0 <- rnbinom(500,mu=2,size=0.5)
y <- rzinbinom(500,mu=2,size=0.5,zprob=0.5)

if (!check_rforge()) {
g0 <- glmmadmb(y0~1,family="nbinom")
logLik(g0)
sum(dnbinom(y0,mu=exp(coef(g0)),size=g0$alpha,log=TRUE))

m1 <- MASS:::fitdistr(y0,"negative binomial")
coef(m1)
logLik(m1)

g1 <- glmmadmb(y0~1,family="nbinom",zeroInflation=TRUE)

logLik(g0)
logLik(g1)

g3 <- glmmadmb(y~1,family="nbinom")
g4 <- glmmadmb(y~1,family="nbinom",zeroInflation=TRUE)

##  library(bbmle)
##  m2 <- mle2(y~dzinbinom(mu=mu,size=alpha,zprob=zprob),
##             method="L-BFGS-B",
##             start=list(mu=2,alpha=0.5,zprob=0.2),
##             lower=rep(0.002,3),
##             upper=c(mu=Inf,alpha=Inf,zprob=0.998),
##             data=data.frame(y))

##  m2P <- profile(m2,which="zprob",std.err=0.025)
##  plot(m2P,show.points=TRUE)

m4 <- fitdistr(y,dzinbinom,start=list(mu=2,size=0.5,zprob=0.5))

ae <- function(x,y,tolerance=1e-2) {
  all.equal(x,y,check.attr=FALSE,tolerance=tolerance)
}

if (.Platform$OS.type=="unix") {
    stopifnot(ae(unname(exp(coef(g4))),unname(coef(m4)["mu"])),
              ae(g4$alpha,unname(coef(m4)["size"])),
              ae(g4$pz,unname(coef(m4)["zprob"]),tol=3e-2))
}
}
