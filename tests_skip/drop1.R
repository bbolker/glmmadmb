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
drop1(g1,test="Chisq")
}
