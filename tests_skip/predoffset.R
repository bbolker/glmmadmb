## test with lm() and offset
library(glmmADMB)

set.seed(1001)
d <- data.frame(x=1:20,y=runif(20),s=rlnorm(20))
L0 <- lm(y~x,data=d)
L1 <- lm(y~x+offset(log(s)),data=d)

p0 <- unname(predict(L0))
p1 <- unname(predict(L1))

if (!check_rforge()) {
L2 <- glmmadmb(y~x,family="gaussian",data=d)
d <- transform(d,logS=log(s))
L3 <- glmmadmb(y~x+offset(logS),
               family="gaussian",data=d)
L4 <- glmmadmb(y~x+offset(log(s)),
               family="gaussian",data=d)

## no-offset, glmmADMB vs lm
stopifnot(all.equal(predict(L2),p0,tolerance=1e-5))
## offset as variable,
stopifnot(all.equal(predict(L3),predict(L3,newdata=d)))
stopifnot(all.equal(predict(L3,newdata=d),p1,tolerance=5e-6))
stopifnot(all.equal(predict(L4),predict(L4,newdata=d)))
stopifnot(all.equal(predict(L4,newdata=d),p1,tolerance=5e-6))
predict(L3,newdata=data.frame(x=1:3,logS=0:2))
predict(L4,newdata=data.frame(x=1:3,s=1:3))

}
