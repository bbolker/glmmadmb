library(glmmADMB)

## Dobson (1990) Page 93: Randomized Controlled Trial :
d <- data.frame(counts=c(18,17,15,20,10,20,25,13,12),
                outcome=gl(3,1,9),
                treatment=gl(3,3))
glm.D93 <- glm(counts ~ outcome + treatment,
               data=d,
               family=poisson())

if (FALSE) {
##    library("ggplot2")
##    library("plyr")
    gmean <- function(x,...) exp(log(mean(x,...)))
    treatment_mean <- ddply(d,"treatment",summarise,counts=gmean(counts))
    outcome_mean <- merge(ddply(d,"outcome",summarise,counts=gmean(counts)),
                  expand.grid(outcome=factor(1:3),treatment=factor(1:3)))
    qplot(outcome,counts,colour=treatment,data=d)+
        geom_line(aes(group=treatment))+
            geom_hline(data=treatment_mean,lty=2,
                       aes(yintercept=counts,colour=treatment))+
           geom_line(data=outcome_mean,lty=3,
                     aes(group=treatment))
}

if (!check_rforge()) {
g2 <- glmmadmb(counts~outcome+treatment,family="poisson",data=d,
               extra.args="-crit 1.e-6")


coef(glm.D93)
coef(g2)

}
