library(glmmADMB)
data("cbpp",package="lme4")

if (!check_rforge()) {
gm1 <- glmmadmb(cbind(incidence, size - incidence) ~ period + (1 | herd),
                   family = "binomial", data = cbpp)

ranef(gm1)
ranef(gm1,sd=TRUE)
ranef(gm1,scale=FALSE)
ranef(gm1,sd=TRUE,scale=FALSE)
}
