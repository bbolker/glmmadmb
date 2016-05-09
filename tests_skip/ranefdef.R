## Make sure that ranef() definition coexists nicely
library(lme4)
gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                   family = binomial, data = cbpp)
ranef(gm1)

library(glmmADMB)

## ranef should now work for glmm.admb fits, should *still*
## work for lme4 models

ranef.x <- function(object,...) print("a")
x <- 1
class(x) <- "x"
try(ranef(x)) ## Error: could not find function "ranef" (not any more!)
library(nlme)
ranef(x)
## [1] "a"

