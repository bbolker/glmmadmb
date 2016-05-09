library(glmmADMB)
library(lme4)

## started out as an attempt to test residuals calculations
## between lme4 and

## would like to do a simpler example (e.g. glm only, from ?glm)
d.AD <- data.frame(counts=c(18,17,15,20,10,20,25,13,12),
                   outcome=gl(3,1,9),
                   treatment=gl(3,3))
glm.D93 <- glm(counts ~ outcome + treatment, family=poisson,
               data=d.AD)
if (!check_rforge()) {
glm.D93.admb <- glmmadmb(counts~outcome+treatment, family="poisson",
          data=d.AD, extra.args="-crit 1.e-6")
r1P <- residuals(glm.D93,type="pearson")
r2P <- residuals(glm.D93.admb,type="pearson")
stopifnot(max(abs(r1P-r2P))<2e-4)
r1R <- residuals(glm.D93,type="response")
r2R <- residuals(glm.D93.admb,type="response")
stopifnot(max(abs(r1R-r2R))<5e-4)

coef(glm.D93)
coef(glm.D93.admb)
## zero-inflated negative binomial is probably best for fitting data:
## glmer can do Poisson-lognormal, not NB, and not zero-inflation

## ... so fit Poisson with RE, even though it's not a great fit to the data

OwlModel_poiss.glmer <- glmer(SiblingNegotiation ~ FoodTreatment * SexParent +
                              (1|Nest)+offset(logBroodSize),
                              data=Owls, family=poisson)


OwlModel_poiss.admb <- glmmadmb(SiblingNegotiation~FoodTreatment*SexParent+
                                (1|Nest)+offset(logBroodSize),
                                data=Owls, family="poisson")

OwlModel_poiss.glmer2 <- glmer(SiblingNegotiation ~ FoodTreatment * SexParent +
                               (1|Nest)+offset(logBroodSize),
                               data=Owls, family=poisson,
                               start=list(fixef=coef(OwlModel_poiss.admb)))

## matrix not pos definite in sparse choleski
## "Estimated covariance matrix may not be positive definite"
logLik(OwlModel_poiss.admb)
logLik(OwlModel_poiss.glmer)

## avoid Suggests: ggplot2 requirement
## if (require(ggplot2)) {
##  ca <- fixef(OwlModel_poiss.glmer)
##  cg <- coef(OwlModel_poiss.admb)
##  d <- rbind(data.frame(par=names(ca),est=ca,pkg="lme4"),
##             data.frame(par=names(cg),est=cg,pkg="glmmADMB"))
##  levels(d$par) <- c("food","food:sex","intercept","sex")
##  qplot(est,par,data=d,colour=pkg)
## }

r_glmer <- residuals(OwlModel_poiss.glmer)
r_admb <- residuals(OwlModel_poiss.admb)

## plot(r_glmer,r_admb)                       

}
