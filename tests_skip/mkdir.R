## Example from ?glm, ultimately from Dobson (1990) p. 93
d <- data.frame(counts=c(18,17,15,20,10,20,25,13,12),
                outcome=gl(3,1,9),
                treatment <- gl(3,3))

library(glmmADMB)
if (!check_rforge()) {
g2 <- glmmadmb(counts~outcome+treatment,family="poisson",data=d,
               extra.args="-crit 1.e-6")
g2 <- glmmadmb(counts~outcome+treatment,family="poisson",data=d,
               extra.args="-crit 1.e-6",
               save.dir="dirtst")
v <- unlink("dirtst",recursive=TRUE)
stopifnot(v==0)
}
