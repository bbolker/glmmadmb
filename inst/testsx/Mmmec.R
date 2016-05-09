library(mlmRev)
library(lme4)
library(glmmADMB)

(fm1 <- lmer(deaths ~ uvb + (1|region), Mmmec, poisson,
             offset = log(expected)))

Mmmec$logExpected <- log(Mmmec$expected)
gm1 <- glmmadmb(deaths ~ uvb + (1|region)+offset(logExpected), Mmmec, "poisson")

coef(gm1)
fixef(fm1)

## overdispersion?  not too bad
sum(residuals(fm1)^2)/(nrow(Mmmec)-2)

gm2 <- glmmadmb(deaths ~ uvb + (1|region)+offset(logExpected), Mmmec, "nbinom")


AIC(gm1)-AIC(gm2)
logLik(gm2)-logLik(gm1)
## large improvement in LL, but little difference in coefficients
coef(gm2)
