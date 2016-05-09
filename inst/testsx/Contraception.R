library(mlmRev)
library(glmmADMB)
library(lme4)

fm1 <- glmer(use ~ urban+age+livch+(1|district), Contraception, binomial)
fm2 <- glmer(use ~ urban+age+livch+(urban|district), Contraception, binomial)

gm1 <- glmmadmb(use ~ urban+age+livch+(1|district), Contraception, "binomial")
gm2 <- glmmadmb(use ~ urban+age+livch+(urban|district), Contraception, "binomial")

cbind(fixef(fm1),coef(gm1))
cbind(fixef(fm2),coef(gm2))

