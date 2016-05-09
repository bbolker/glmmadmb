library(glmmADMB)

set.seed(1002)
nblock <- 10
nperblock <- 50
sd.u <- 1
ntot <- nblock*nperblock
d <- data.frame(x=runif(ntot),f=factor(rep(LETTERS[1:nblock],each=nperblock)))
r <- rnorm(nblock,mean=0,sd=sd.u)
gshape <- 1.5
d$eta <- with(d,0.2+0.5*x+r[f])
d$mu <- exp(d$eta)
d$y <- rgamma(ntot,shape=gshape,scale=d$mu/gshape)

if (!check_rforge()) {

g1 <- glmmadmb(y~x+(1|f),data=d,family="gamma")
g1L <- glmmPQL(y~x,random=~1|f,data=d,family=Gamma(link="log"))

coef(g1)
fixef(g1L)
VarCorr(g1)
VarCorr(g1L)

d2 <- d
d2$eta <- with(d2,0.2+0.5*x) ## 0.5 ## 
d2$mu <- exp(d2$eta)
d2$y <- rgamma(ntot,shape=gshape,scale=d2$mu/gshape)

g0 <- glmmadmb(y~x,data=d2,family="gamma")
g0L <- glm(y~x,data=d2,family=Gamma(link="log"))
coef(g0)
coef(g0L)
}
