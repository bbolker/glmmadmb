library("glmmADMB")

set.seed(1002)
nblock <- 10
nperblock <- 50
sd.u <- 0.2
ntot <- nblock*nperblock
d <- data.frame(x=runif(ntot),f=factor(rep(LETTERS[1:nblock],each=nperblock)))
r <- rnorm(nblock,mean=0,sd=sd.u)
phi <- 2
d$eta0 <- with(d,-3+5*x)
d$eta <- with(d,eta0+r[f])
d$mu0 <- plogis(d$eta0)
d$mu <- plogis(d$eta)
theta <- 2
shape1 <- theta * d$mu
shape2 <- theta * (1 - d$mu)
shape1_0 <- theta * d$mu0
shape2_0 <- theta * (1 - d$mu0)
d$y0 <- rbinom(ntot, size = 20, prob = rbeta(ntot, shape1_0, shape2_0))
d$y <- rbinom(ntot, size = 20, prob = rbeta(ntot, shape1, shape2))
if (FALSE) {
    with(d,plot(x,y0))
    with(d,curve(20*plogis(-3+5*x),col=2,add=TRUE))
}
if (FALSE) {
  ## comment out until binaries are updated on all platforms
g0 <- glmmadmb(cbind(y0,20-y0)~x,data=d,family="betabinomial")

coef(g0)
summary(g0)
if (FALSE) {
    g1 <- glmmadmb(cbind(y,20-y)~x+(1|f),data=d,family="betabinomial",
                   save.dir="bbtmp")
    coef(g1)
    VarCorr(g1)
    summary(g1)
}
}
