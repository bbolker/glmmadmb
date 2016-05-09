library(glmmADMB)


rnbinom1 <- function(n,mu,size) {
  ##    (1+mu/size2) = size
  ## -> size2 = mu/(size-1)
  size2 <- mu/(size-1)
  ## cat(size,size2,mu*size,mu*(1+mu/size2),"\n")
  rnbinom(n,mu=mu,size=size2)
}

set.seed(1002)

## vvec <- replicate(500,var(rnbinom1(500,mu=3,size=2)))
y1 <- rnbinom1(500,mu=3,size=2)
c(mean(y1),var(y1)) ## approx. 6?

if (!check_rforge()) {
g0 <- glmmadmb(y1~1,family="nbinom")
g1 <- glmmadmb(y1~1,family="nbinom1")
g0$alpha
g1$alpha

x <- runif(500,min=-1,max=4)
y2 <- rnbinom1(500,mu=exp(1+x),size=2)

g2 <- glmmadmb(y2~x,family="nbinom")
g3 <- glmmadmb(y2~x,family="nbinom1")
AIC(g2)-AIC(g3)  ## >0; nbinom1 is better

y3 <- rnbinom(500,mu=exp(1+x),size=2)
g4 <- glmmadmb(y3~x,family="nbinom")
g5 <- glmmadmb(y3~x,family="nbinom1")
AIC(g4)-AIC(g5) ## <0; nbinom is better
}
