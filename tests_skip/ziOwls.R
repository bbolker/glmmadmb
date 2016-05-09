library(glmmADMB)
testLevel <- if (nzchar(s <- Sys.getenv("GLMMADMB_TEST_LEVEL"))) as.numeric(s) else 1

## slow, only run when requested
if (testLevel>1) {
Owls <- transform(Owls,
                  Nest=reorder(Nest,NegPerChick),
                  logBroodSize=log(BroodSize),
                  NCalls=SiblingNegotiation)
if (!check_rforge()) {
fit_zipoiss <- glmmadmb(NCalls~(FoodTreatment+ArrivalTime)*SexParent+
                        offset(logBroodSize)+(1|Nest),
                        data=Owls,
                        zeroInflation=TRUE,
                        family="poisson")
}
sessionInfo()
}
