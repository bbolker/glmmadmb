library(glmmADMB)  ## testing version

## make sure file doesn't automatically get run during testing ...
if (file.exists("multrand_batch.RData")) {
  load("multrand_batch.RData")
} else {
  set.seed(101)
  ## nblock <- 10
  ## nrep <- 10

  ## smaller (5 x 50) for quicker testing
  simfun <- function(seed=101,
                   nblock=5,
                   nrep=50,
                   rsd=c(f=1,g=2),
                   beta=c(1,2)) {
  if (!is.null(seed)) set.seed(seed)
  d <- expand.grid(f=factor(LETTERS[1:nblock]),
                   g=factor(letters[1:nblock]),
                   rep=1:nrep)
  N <- nrow(d)
  d$x <- runif(N)
  u_f <- rnorm(nblock,sd=rsd["f"])
  u_g <- rnorm(nblock,sd=rsd["g"])
  eta <- model.matrix(~x,data=d) %*% beta +
    u_f[as.numeric(d$f)]+u_g[as.numeric(d$g)]
  d$y <- rpois(N,exp(eta))
  attr(d,"reff") <- rbind(data.frame(eff="f",block=levels(d$f),u=u_f),
                          data.frame(eff="g",block=levels(d$g),u=u_g))
  d
}
                          
                               
library(lme4)
## kluge, for passing tests until I can get this sorted out
setMethod("VarCorr", signature(x="glmmadmb"), glmmADMB:::VarCorr.glmmadmb)
setMethod("VarCorr", signature(x="summary.glmmadmb"), glmmADMB:::VarCorr.glmmadmb)

d1 <- simfun()
t1_lme4 <- system.time(g1_lme4 <- glmer(y~x+(1|f)+(1|g),family="poisson",data=d1))
t1_GA <- system.time(g1_GA <- glmmadmb(y~x+(1|f)+(1|g),family="poisson",data=d1))

save.image("multrand_batch.RData")
###############

d2 <- simfun(nblock=15)
t2_lme4 <- system.time(g2_lme4 <- glmer(y~x+(1|f)+(1|g),family="poisson",data=d2))
t2_GA <- system.time(g2_GA <- glmmadmb(y~x+(1|f)+(1|g),family="poisson",data=d2))
save.image("multrand_batch.RData")

## with nblock=15: works, but takes 880 seconds
## with nblock=20:
##   Need to increase the maximum number of separable calls allowed to at least 20001
##   Current value is 20000
}
