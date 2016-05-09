library(glmmADMB)
## source("glmmadmb.R")

options(digits=3)
set.seed(1002)
nblock <- 10
nperblock <- 20
sd.u <- 1
ntot <- nblock*nperblock
d <- data.frame(x=runif(ntot),f=factor(rep(LETTERS[1:nblock],each=nperblock)))
r <- rnorm(nblock,mean=0,sd=sd.u)
d$eta <- with(d,0.2+0.5*x+r[f])
d$mu <- exp(d$eta)
d$y <- rpois(ntot,lambda=d$mu)

if (!check_rforge()) {
g1 <- glmmadmb(y~x+(1|f),family="poisson",data=d)
coef(g1)
VarCorr(g1)

if (.Platform$OS.type=="unix") {
g1M <- glmmadmb(y~x+(1|f),family="poisson",data=d,mcmc=TRUE,
                 mcmc.opts=mcmcControl(mcmc=100))

library(coda)
HPDinterval(g1M$mcmc)
summary(g1M$mcmc)
fixef(g1M)
## xyplot(as.mcmc(g1M$mcmc),layout=c(4,4),aspect="fill")
## densityplot(as.mcmc(g1M$mcmc),layout=c(4,4),aspect="fill")

### try example from simon.chamaille@cefe.cnrs.fr

nblock <- 8
indperblock <- 10
nperind <- 5
sd.u <- 1
ntot <- nblock*indperblock*nperind
d <- expand.grid(f=factor(rep(LETTERS[1:nblock],each=nperblock)),
                 ind=factor(rep(1:indperblock)),
                 rep=rep(1:nperind))
d$x <- runif(ntot)
u1 <- rnorm(nblock,mean=0,sd=sd.u)
u2 <- rnorm(nblock*indperblock,mean=0,sd=sd.u)
d$eta <- with(d,0.2+0.5*x+u1[f]+u2[interaction(f,ind)])
d$mu <- exp(d$eta)
d$y <- rpois(ntot,lambda=d$mu)

if (FALSE) {
    ## slow!
    mod.admb <- glmmadmb(formula=y~1+(1|f/ind),
                         data=d,
                         family="nbinom1",link="log",
                         admb.opts=admbControl(shess=FALSE,noinit=FALSE),
                         mcmc=TRUE,mcmc.opts=mcmcControl(mcmc=100),
                         extra.args="-ndi 60000",
                         verbose=TRUE)
}
}
}
