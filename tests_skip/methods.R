library(glmmADMB)
methods(class="glmmadmb")
getdata("OwlModel")
AIC(OwlModel)
coef(OwlModel)
## deviance(OwlModel) ## NULL
df.residual(OwlModel)
fitted(OwlModel)
fixef(OwlModel)
logLik(OwlModel)
head(model.frame(OwlModel))
nobs(OwlModel)
head(predict(OwlModel))
predict(OwlModel,newdata=data.frame(model.frame(OwlModel)[1:5,],BroodSize=5))
ranef(OwlModel)
head(residuals(OwlModel))
stdEr(OwlModel)
summary(OwlModel)
VarCorr(OwlModel)
vcov(OwlModel)

## AICc: implemented in MuMIn, AICcmodavg, bbmle
## glmmADMB objects work with MuMIn and bbmle versions, not
##   with AICcmodavg
## FIXME: write to AICcmodavg maintainer?
## coeftab: in coefplot2
