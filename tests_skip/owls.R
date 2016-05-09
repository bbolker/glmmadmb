data(Owls)
g1 <- glmmadmb(SiblingNegotiation~1,family="poisson",data=Owls)
g2 <- glmmadmb(SiblingNegotiation~1,family="poisson",data=Owls,
               admb.opts=admbControl(poiss_prob_bound=0))
c0 <- log(mean(Owls$SiblingNegotiation))  ## correct value
g3 <- glm(SiblingNegotiation~1,family="poisson",data=Owls)
all.equal(c2 <- unname(coef(g2)),c0)
