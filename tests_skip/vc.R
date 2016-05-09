## tests of nesting as well 
library(glmmADMB)

## try to make reporting work for (mistakenly) nested variables, without breaking anything:

##
set.seed(10)
d <- expand.grid(f1=LETTERS[1:5],f2=letters[1:3],rep=1:5)
d <- transform(d,f12=interaction(f1,f2),y=rnorm(nrow(d)),x=rnorm(nrow(d)),x2=rnorm(nrow(d)))
## library(nlme)
## lme(y~1, random=~1|f1/f2,data=d))
## lme(y~1, random=~1|f1/f12,data=d))

if (!check_rforge()) {
t1 <- system.time(g1 <- glmmadmb(y~1+(1|f1/f2),family="gaussian",data=d))
t2 <- system.time(g2 <- glmmadmb(y~1+(1|f1/f12),family="gaussian",data=d))
VarCorr(g1)
## except for reporting, we get the same answers
summary(g1)
summary(g2)


## need another test with another component:

## FIXME: does this model make sense?  feels like we have a (1|f1) term repeated ...
g3 <- glmmadmb(y~1+(1|f1/f2)+(x|f1),family="gaussian",data=d)
VarCorr(g3)

if (FALSE) {
    ## FIXME: these don't work, but because of non-pos-def errors rather
    ##  than anything essentially wrong?
    g4 <- glmmadmb(y~1+(x|f1),corStruct="full",family="gaussian",data=d)
    VarCorr(g4)

    g5 <- glmmadmb(y~1+(x+x2|f1),corStruct="full",family="gaussian",data=d)
    VarCorr(g5)
}

## m5 <- lme(y~1,random=~x+x2|f1,data=d)
## needs work.

## is doing it 'wrong' (explicit nesting) slower?
## library(rbenchmark)
## benchmark(glmmadmb(y~1+(1|f1/f2),family="gaussian",data=d),
##           glmmadmb(y~1+(1|f1/f12),family="gaussian",data=d),
##          replications=20)

## only a little bit, at least for this example:
 
##         n   elapsed relative user.self sys.self user.child sys.child
##  f1/f12 20  39.901 1.166935     0.572    2.344     30.154     6.317
##  f1/f2  20  34.193 1.000000     0.432    2.200     26.337     4.848

## }
}
