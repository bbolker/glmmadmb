library(glmmADMB)
library(lme4)

gm1_lme4 <- glmer(cbind(incidence, size - incidence) ~
                  period + (1 | herd), data = cbpp, family = binomial)
if (!check_rforge()) {
gm1_glmmADMB <- glmmadmb(cbind(incidence, size -
                               incidence) ~ period + (1 | herd),
                         data = cbpp, family = "binomial")

## this is NECESSARY for right now
glmmADMB:::VarCorr(gm1_glmmADMB)
try(VarCorr(gm1_glmmADMB))
fixef(gm1_lme4)
}
