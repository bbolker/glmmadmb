library(glmmADMB)  ## testing version

ss <- system.file("testdata","multrand_batch.RData",package="glmmADMB")
load(ss)
ss2 <- system.file("testdata","multrand_lme4_sum.RData",package="glmmADMB")
load(ss2)

## kluge, for passing tests until I can get this sorted out
setMethod("VarCorr", signature(x="glmmadmb"), glmmADMB:::VarCorr.glmmadmb)
setMethod("VarCorr", signature(x="summary.glmmadmb"), glmmADMB:::VarCorr.glmmadmb)

sumfun <- function (x,times)
  UseMethod("sumfun")

tmpnamefun <- function(nn) {
  c(outer(c("min","mean","max"),nn,
          function(x,y) paste(y,x,sep=".")))
}
tmpsumfun <- function(x) c(min=min(x),mean=mean(x),max=max(x))

## will not work for 'old' glmmADMB (wrong ranef structure)
sumfun.glmmadmb <- function(x,times) {
  fixed <- coef(x)
  ransum <- unlist(lapply(ranef(x),
                   function(z)
                   apply(z,MARGIN=2,FUN=tmpsumfun)))
  LL <- logLik(x)
  rv <- unlist(lapply(VarCorr(x),diag))
  times <- round(times[3],2)
  mm <- c(fixed,c(rv),c(LL),ransum,times)
  rnames <- names(rv)
  names(mm) <- c(names(coef(x)),
                 paste("var(RE)",rnames,sep="."),
                 "logLik",
                 paste("U",tmpnamefun(rnames),sep="."),
                 "time")
  mm
}


## sumfun.mer <- sumfun.merMod <- function(x,times) {
##     ## should work for old/new lme4 ...
##   fixed <- fixef(x)
##   ransum <- unlist(lapply(ranef(x),
##                           function(z)
##                           apply(z,MARGIN=2,FUN=tmpsumfun)))
##   LL <- logLik(x)
##   rv <- sapply(VarCorr(x),c)
##   times <- round(times[3],2)
##   mm <- c(fixed,c(rv),c(LL),ransum,times)
##   rnames <- names(rv)
##   names(mm) <- c(names(fixef(x)),
##                  paste("var(RE)",rnames,sep="."),
##                  "logLik",
##                  paste("U",tmpnamefun(rnames),sep="."),
##                  "time")
##   mm
## }

## g1_lme4_sum <- sumfun(g1_lme4,t1_lme4)
## g2_lme4_sum <- sumfun(g2_lme4,t2_lme4)
## save("g1_lme4_sum","g2_lme4_sum",file="multrand_lme4_sum.RData")

## sumfun2A <- function(modlist,tlist) {
##   mapply(sumfun,modlist,tlist)
## }

## does this work for multiple grouping variables?
g1_GA_sum <- sumfun(g1_GA,t1_GA)
g2_GA_sum <- sumfun(g2_GA,t2_GA)

cmpfun <- function(x,y,tol=1e-2,ignore="time",
                   meantol=2e-2) {
    cdif <- abs(x/y-1)
    cdif <- cdif[!names(cdif) %in% ignore]
    mvec <- grepl("mean$",names(cdif))
    all(cdif[!mvec]<tol && cdif[mvec]<meantol)
}
stopifnot(cmpfun(g1_GA_sum,g1_lme4_sum))
stopifnot(cmpfun(g2_GA_sum,g2_lme4_sum,meantol=0.2))
## this is all a bit silly since these tests don't get re-run anyway...

