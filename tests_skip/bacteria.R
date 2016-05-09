data(bacteria,package="MASS")

bacteria$present <- as.numeric(bacteria$y)-1

## FIXME: restore when problems with new lme4 are sorted!
if (FALSE) {
    library(lme4)
    gfit <- glmer(present ~ trt + I(week > 2)+(1|ID),family = "binomial", data = bacteria,
                  optimizer="bobyqa", tolPwrss=1e-13, verbose=10)
}

library(glmmADMB)
if (!check_rforge()) {
bfit <-  glmmadmb(present ~ trt + I(week > 2), random = ~ 1 | ID,
                     family = "binomial", data = bacteria)
if (!is.null(bfit$phi)) {
stopifnot(all.equal(bfit$phi,
                    glmmADMB:::make_phi(model.matrix(~trt+I(week>2),data=bacteria)),
                                        tol=1e-4))
}
}

if (FALSE) {
    gsd <- attr(lme4::VarCorr(gfit)[[1]],"stddev")
    f <- lme4::fixef(gfit)
    u <- unlist(lme4::ranef(gfit)[[1]])/gsd
    bfit2 <-glmmadmb(present ~ trt + I(week > 2), random = ~ 1 | ID,
                     family = "binomial", data = bacteria,
                 start =list(fixed=f,RE_sd=log(gsd),u=u),
                 verbose=TRUE)
slist <- list(fixed=f,RE_sd=log(gsd),u=u)
bfit3 <-glmmadmb(present ~ trt + I(week > 2), random = ~ 1 | ID,
                     family = "binomial", data = bacteria,
                 start=slist,
                 extra.args="-phase 5",
                 verbose=TRUE)
}

## 9.6130683
