apa <- data.frame(odist=ordered(sample(1:5,replace=TRUE,size=100)),
            sl=runif(100),abun=rnbinom(100,mu=2,size=1))

library(glmmADMB)
if (!check_rforge()) {
gg <- try(glmmadmb(abun~odist+(1|sl),data=apa,zeroInflation=TRUE,family="nbinom"),silent=TRUE)
## should throw an error from Droplevels
stopifnot(inherits(gg,"try-error"),grepl("Error in Droplevels",gg))
}
