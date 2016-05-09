library(glmmADMB)

set.seed(1002)
nblock <- 10
nperblock <- 50
sd.u <- 1
ntot <- nblock*nperblock
d <- data.frame(x=runif(ntot),f=factor(rep(LETTERS[1:nblock],each=nperblock)))
r <- rnorm(nblock,mean=0,sd=sd.u)
phi <- 2
d$eta <- with(d,0.2+0.5*x+r[f])
d$mu <- plogis(d$eta)
d$y <- rbeta(ntot,shape1=d$mu*phi,shape2=(1-d$mu)*phi)

if (!check_rforge()) {
    g1 <- glmmadmb(y~x+(1|f),data=d,family="beta")

coef(g1)
VarCorr(g1)
summary(g1)
}

