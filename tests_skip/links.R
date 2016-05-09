## testing identity, cloglog links, gaussian family
library(glmmADMB)

set.seed(1002)
nblock <- 10
nperblock <- 50
sd.u <- 1
beta <- c(0.2,0.5)
ntot <- nblock*nperblock
d <- data.frame(x=runif(ntot),f=factor(rep(LETTERS[1:nblock],each=nperblock)))
r <- rnorm(nblock,mean=0,sd=sd.u)
gshape <- 1.5
d$offset <- rgamma(ntot,1,1)
d <- within(d,
            {
              eta0 <- beta[1]+beta[2]*x+offset
              eta <- eta0+r[f]
            })


## cloglog:
cc <- binomial(link="cloglog")
d$mu0 <- cc$linkinv(d$eta0)
d$mu <- cc$linkinv(d$eta)

d$y0 <- rbinom(ntot,prob=d$mu0,size=1)
d$y <- rbinom(ntot,prob=d$mu,size=1)

if (!check_rforge()) {
## A. no random effects (vs glm)
g0 <- glmmadmb(y0~x+offset(offset),data=d,
               family="binomial",link="cloglog")

g0A <- glm(y0~x+offset(offset),data=d,
           family=binomial(link="cloglog"))

mlist <- list(glmmadmb0=g0,glm0=g0A)

stopifnot(all.equal(coef(g0),coef(g0A),tol=1e-4))

t(sapply(mlist,coef))
sapply(mlist,logLik)

## B. random effects (vs glmer)
g1 <- glmmadmb(y~x+(1|f)+offset(offset),data=d,
               family="binomial",link="cloglog")


p0 <- predict(g1)
pb <- predict(g1,type="response")


### GAMMA/LOG LINK
gshape <- 1.5

cc <- Gamma(link="log")
d$mu0 <- cc$linkinv(d$eta0)
d$mu <- cc$linkinv(d$eta)

d$y0 <- rgamma(ntot,shape=gshape,scale=d$mu0/gshape)
d$y <- rgamma(ntot,shape=gshape,scale=d$mu/gshape)

## glmmadmb vs glm, no random effect
g2 <- glmmadmb(y0~x,data=d,family="gamma",link="log")
g2L <- glm(y0~x,data=d,family=Gamma(link="log"))

stopifnot(all.equal(coef(g2),coef(g2L),tol=1e-5))

g3 <- glmmadmb(y~x+(1|f),data=d,family="gamma",link="log")
## "matrix not pos definite in sparse choleski" warning

## POISSON/identity link
## FIXME: allow on Windows
## FIXME: fails on Fedora 18 (A. Magnusson) ?
if (.Platform$OS.type=="unix") {
  dd <- data.frame(y=rpois(20,lambda=10),f=factor(rep(1:5,each=4)))
   g5 <- glmmadmb(y~1,data=dd,
         start=list(fixed=10),
         family="poisson",link="identity")
   g5R <- glmmadmb(y~1+(1|f),data=dd,
         start=list(fixed=10),
         family="poisson",link="identity")

   coef(g5)
coef(g5R)

### GAUSSIAN/IDENTITY LINK
d$mu0 <- d$eta0
d$mu <- d$eta

d$y0 <- rnorm(ntot,d$mu0,sd=1)
d$y <- rnorm(ntot,d$mu,sd=1)

g4 <- lm(y0~x,data=d)
g4B <- glm(y0~x,data=d)
g4C <- glmmadmb(y0~x,data=d,family="gaussian")
stopifnot(all.equal(coef(g4),coef(g4B),coef(g4C)))
stopifnot(all.equal(c(logLik(g4)),c(logLik(g4B)),c(logLik(g4C)),tol=1e-4))

}
}
