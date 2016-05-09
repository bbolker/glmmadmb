library(glmmADMB)
library(bbmle)
set.seed(101)

n <- 1000
d <- data.frame(x=runif(n))
d <- transform(d,y=rlogis(n,1+2*x,0.5))

m1 <- mle2(y~dlogis(m,exp(logs)),
           parameters=list(m~x),
           data=d,
           start=list(m=0,logs=0))
## mle2 FIXME: need better error message for
## Error in if (vpos0 == 1) 1 else vposvals[vpos0 - 1] + 1 : 
##   argument is of length zero
##   start=list(x=0,logs=0))
coef(m1)

if (!check_rforge()) {
m2 <- glmmadmb(y~x,data=d,family="logistic")
coef(m2)

}
## glmmadmb FIXME: bogus printout in number of observations
##  should return sigma/log-sigma info

