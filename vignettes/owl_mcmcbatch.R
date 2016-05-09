library(glmmADMB)
Owls <- transform(Owls,
  Nest=reorder(Nest,NegPerChick),
  logBroodSize=log(BroodSize),
  NCalls=SiblingNegotiation)

fit_zinbinom1_bs_mcmc <- glmmadmb(NCalls~(FoodTreatment+ArrivalTime)*SexParent+
                                     BroodSize+(1|Nest),
                                     data=Owls,
                                     zeroInflation=TRUE,
                                     family="nbinom1",
                                  mcmc=TRUE,
                                  mcmc.opts=mcmcControl(mcmc=50000,thin=50))

save("fit_zinbinom1_bs_mcmc",file="owls_mcmcbatch.RData")
