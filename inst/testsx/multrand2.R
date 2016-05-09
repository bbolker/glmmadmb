set.seed(1001)
nb <- c(5,8,10)
d <- data.frame(x=runif(prod(nb)),expand.grid(f1=LETTERS[1:nb[1]],f2=letters[1:nb[2]]))
d$inter <- interaction(d$f1,d$f2)
sdvec <- c(2,1)
u1 <- rnorm(nb[1],sd=sdvec[1])
u2 <- rnorm(nb[1]*nb[2],sd=sdvec[2])
d$eta <- with(d,1+0.5*x+u1[f1]+u2[inter])
d$y <- rpois(prod(nb),exp(d$eta))

d <- subset(d,select=c(x,f1,f2,inter,y))

library(lme4)
t1 <- system.time(g1 <- glmer(y~x+(1|f1/f2),data=d,family=poisson))
t1B <- system.time(g1B <- glmer(y~x+(1|f1)+(1|inter),data=d,family=poisson))
fix1 <- fixef(g1)
fix1B <- fixef(g1B)
all.equal(fix1,fix1B)
ran1 <- ranef(g1)
ran1B <- ranef(g1B)
all(unlist(ran1)-unlist(ran1B)==0)

## 1 second
## NB 'fixef' in glmmADMB screws up/masks the one in lme4 ...

library(glmmADMB)
t2 <- system.time(g2 <- glmmadmb(y~x+(1|f1/f2),data=d,family="poisson"))
## 25-30 seconds: warning about non-pos-dev cov matrix
t2B <- system.time(g2B <- glmmadmb(y~x+(1|f1)+(1|inter),data=d,family="poisson"))
## about the same time
t2C <- system.time(g2C <- glmmadmb(y~x+(1|f1/inter),data=d,family="poisson"))

fix2 <- fixef(g2)
fix2B <- fixef(g2B)
fix2C <- fixef(g2C)

all.equal(fix2,fix2B)
ran2 <- ranef(g2)
ran2B <- ranef(g2B)
ran2C <- ranef(g2C)
## these end up in a different order but are otherwise identical ...
all(sort(unlist(ran2))-sort(unlist(ran2B))==0)

max(abs(fix1-fix2))
max(abs(ran1[[2]]-ran2[[1]]))  ## different order
max(abs(sort(unlist(ran1[[1]]))-sort(ran2[[2]])))  ## different order

## NB names of 'ran2' need work
