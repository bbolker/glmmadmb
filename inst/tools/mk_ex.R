## code for exporting input/output files for glmmadmb testing
library(glmmADMB)

save_fun <- function(savename,argList) {
    dir.create(savename)
    argList <- c(argList,list(save.dir=savename))
    do.call("glmmadmb",argList)
    setwd(savename)
    for (i in c("pin","par","dat","std")) {
        file.copy(paste("glmmadmb",i,sep="."),
                  paste0("../",savename,".",i))
    }
    setwd("..")
    unlink(savename,recursive=TRUE,force=TRUE)
}

## bacteria
data(bacteria,package="MASS")
bacteria$present <- as.numeric(bacteria$y)-1
save_fun("bacteria",
         list(formula=present ~ trt + I(week > 2),
              random = ~ 1 | ID,
              family = "binomial", data = bacteria))

## beta-binomial
set.seed(1002)
nblock <- 10
nperblock <- 50
sd.u <- 0.2
ntot <- nblock*nperblock
d <- data.frame(x=runif(ntot),f=factor(rep(LETTERS[1:nblock],each=nperblock)))
r <- rnorm(nblock,mean=0,sd=sd.u)
phi <- 2
d$eta0 <- with(d,-3+5*x)
d$eta <- with(d,eta0+r[f])
d$mu0 <- plogis(d$eta0)
d$mu <- plogis(d$eta)
theta <- 2
shape1 <- theta * d$mu
shape2 <- theta * (1 - d$mu)
shape1_0 <- theta * d$mu0
shape2_0 <- theta * (1 - d$mu0)
d$y0 <- rbinom(ntot, size = 20, prob = rbeta(ntot, shape1_0, shape2_0))
d$y <- rbinom(ntot, size = 20, prob = rbeta(ntot, shape1, shape2))

save_fun("betabinom",
        list(formula=cbind(y0,20-y0)~x,data=d,family="betabinomial"))


data("bioChemists",package="pscl")
bb <- subset(bioChemists,art>0)
save_fun("hurdle_biochem1",
         list(formula=art~fem+mar+kid5+phd+ment,
              family="truncpoiss",link="log",data=bb))
bioChemists <- transform(bioChemists,nz=as.numeric(art>0))
save_fun("hurdle_biochem2",
         list(formula=nz~fem+mar+kid5+phd+ment,
               family="binomial",data=bioChemists))

Owls <- transform(Owls,
                  Nest=reorder(Nest,NegPerChick),
                  logBroodSize=log(BroodSize),
                  NCalls=SiblingNegotiation)
save_fun("zi_owls",
         list(NCalls~(FoodTreatment+ArrivalTime)*SexParent+
              offset(logBroodSize)+(1|Nest),
              data=Owls,
              zeroInflation=TRUE,
              family="poisson"))

epil2$subject <- factor(epil2$subject)
save_fun("epil2",
         list(y~Base*trt+Age+Visit+(Visit|subject),
              data=epil2, family="nbinom"))
