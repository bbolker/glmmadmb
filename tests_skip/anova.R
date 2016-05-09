library(glmmADMB)
data(bacteria,package="MASS")
bacteria$present <- as.numeric(bacteria$y)-1
if (!check_rforge()) {
    bfit <-  glmmadmb(present ~ trt + I(week > 2), random = ~ 1 | ID,
                      family = "binomial", data = bacteria)
    bfit2 <- update(bfit, .~.-trt)
    anova(bfit,bfit2)
    anova(bfit2,bfit)
}
