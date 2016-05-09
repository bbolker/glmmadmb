library(glmmADMB)
Owls <- transform(Owls,
                  Nest=reorder(Nest,NegPerChick),
                  NCalls=SiblingNegotiation)

time0 <- system.time(OwlModel <- glmmadmb(NCalls~(FoodTreatment+ArrivalTime)*SexParent+
                                          offset(log(BroodSize))+(1|Nest),
                                          data=Owls,
                                          zeroInflation=TRUE,
                                          family="nbinom"))
save(list=c(ls(pattern="time"),ls(pattern="OwlModel")),file="OwlModel.rda")

time1 <- system.time(OwlModel_nb1_bs <- glmmadmb(NCalls~(FoodTreatment+ArrivalTime)*SexParent+
                                                 BroodSize+(1|Nest),
                                                 data=Owls,
                                                 zeroInflation=TRUE,
                                                 family="nbinom1"))
save(list=c(ls(pattern="time"),ls(pattern="OwlModel")),file="OwlModel.rda")

time2 <- system.time(OwlModel_nb1_bs_mcmc <- glmmadmb(NCalls~(FoodTreatment+ArrivalTime)*SexParent+
                                                      BroodSize+(1|Nest),
                                                      data=Owls,
                                                      zeroInflation=TRUE,
                                                      family="nbinom1",
                                                      mcmc=TRUE,
                                                      mcmc.opts=mcmcControl(mcmc=50000)))

save(list=c(ls(pattern="time"),ls(pattern="OwlModel")),file="OwlModel.rda")
file.copy("OwlModel.rda","../inst/extdata/OwlModel.rda",overwrite=TRUE)

if (FALSE) {
  load("owls_mcmcbatch.RData")
  fit0 <- glmmadmb(NCalls~(FoodTreatment+ArrivalTime)*SexParent+
                   BroodSize+(1|Nest),
                   data=Owls,
                   zeroInflation=TRUE,
                   family="nbinom1",
                   save.dir="owls_mcmcdir0",
                   mcmc=TRUE,
                   admb.opts=admbControl(run=FALSE))
  
  fit_zinbinom1_bs_mcmc <- glmmadmb(NCalls~(FoodTreatment+ArrivalTime)*SexParent+
                                    BroodSize+(1|Nest),
                                    data=Owls,
                                    zeroInflation=TRUE,
                                    family="nbinom1",
                                    save.dir="owls_mcmcdir",
                                    mcmc=TRUE,
                                    mcmc.opts=mcmcControl(mcmc=50000),
                                    admb.opts=admbControl(run=FALSE))

  
  ## In matrix(readBin(f, "double", n = (fs - isize)/dsize), byrow = TRUE,  :
  ## data length [16384] is not a sub-multiple or multiple of the number of rows [443]

  save("fit0","time0","fit_zinbinom1_bs_mcmc","mcmctime",file="owls_mcmcbatch.RData")


  fit0 <- fit_zinbinom2_bs_mcmc
  pnames <- colnames(fit0$mcmc)

  load("owls_mcmcbatch.RData")
  colnames(fit_zinbinom1_bs_mcmc$mcmc) <- pnames
  
  save("fit0","fit_zinbinom1_bs_mcmc","mcmctime",file="owls_mcmcbatch.RData")

}
