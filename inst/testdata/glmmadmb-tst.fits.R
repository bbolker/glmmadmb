#### Saved fits for glmmadmb testing
####  copied from lme4
####  ----------------------------------
fn0 <- "glmmadmb-tst-fits.rda"
fn <- system.file("testdata", fn0,
                  package="lme4", mustWork=TRUE)
if(FALSE) ### "Load" these by  load(fn)
    ## or "better"
    attach(fn)

to.save <- character()

library(glmmADMB)
str(packageDescription("glmmADMB")[c("Version", "Packaged", "Built")])
data("sleepstudy",package="lme4")

## intercept only in both fixed and random effects
fit_sleepstudy_0 <- glmmadmb(Reaction ~ 1 + (1|Subject), sleepstudy,
                             family="gaussian")
## fixed slope, intercept-only RE
fit_sleepstudy_1 <- glmmadmb(Reaction ~ Days + (1|Subject), sleepstudy,
                             family="gaussian")
## fixed slope, intercept & slope RE
fit_sleepstudy_2 <- glmmadmb(Reaction ~ Days + (Days|Subject), sleepstudy,
                         family="gaussian",
                         corStruct="full")
## fixed slope, independent intercept & slope RE
fit_sleepstudy_3 <- glmmadmb(Reaction ~ Days + (Days|Subject), sleepstudy,
                         family="gaussian",
                         corStruct="diag") ## (default anyway)

data("cbpp",package="lme4")
cbpp$obs <- factor(seq(nrow(cbpp)))
## intercept-only fixed effect
fit_cbpp_0 <- glmmadmb(cbind(incidence, size-incidence) ~ 1 + (1|herd),
                    cbpp, family="binomial")
## include fixed effect of period
fit_cbpp_1 <- update(fit_cbpp_0, . ~ . + period)  ## matrix not pos def warning
## include observation-level RE
fit_cbpp_2 <- update(fit_cbpp_1, . ~ . + (1|obs))
## specify formula by proportion/weights instead
## fit_cbpp_3 <- update(fit_cbpp_1, incidence/size ~ period + (1 | herd), weights = size)

save(list=c(to.save, ls(pattern="fit_")), file=fn0)
