## basic comparison of a series of glmmADMB models
library(glmmADMB)
testLevel <- if (nzchar(s <- Sys.getenv("GLMMADMB_TEST_LEVEL"))) as.numeric(s) else 1
## slow, only run when requested
if (testLevel>2) {

    set.seed(1001)
    nblock <- 10
    N <- 500
    d <- data.frame(x=runif(N),
                    f=factor(rep(seq(nblock),each=N/nblock)))
    u <- rnorm(nblock,sd=1)
    d$y <- with(d,rnbinom(N,mu=exp(1+2*x+u[f]),size=1.2))
    d$obs <- factor(seq(nrow(d)))

    library(glmmADMB)

    m_poiss <- glmmadmb(y~x+(1|f),family="poisson",data=d)
    m_ZIP <- glmmadmb(y~x+(1|f),family="poisson",data=d,
                      zeroInflation=TRUE,
                      start=list(fixed=fixef(m_poiss)))
    m_LNPoiss <- glmmadmb(y~x+(1|f)+(1|obs),family="poisson",data=d)
    m_NB <- glmmadmb(y~x+(1|f),family="nbinom",data=d)
    m_ZINB <- glmmadmb(y~x+(1|f),family="nbinom",data=d,zeroInflation=TRUE)
    m_NB1 <- glmmadmb(y~x+(1|f),family="nbinom1",data=d)
    m_ZINB1 <- glmmadmb(y~x+(1|f),family="nbinom1",data=d,zeroInflation=TRUE)

    ndf <- function(x) {
        attr(logLik(x),"df")
    }

    nvec <- c("poiss","ZIP","LNPoiss","NB","ZINB","NB1","ZINB1")
    mlist <- lapply(paste("m_",nvec,sep=""),get)
    names(mlist) <- nvec
    sapply(mlist,ndf)
    library(bbmle)
    AICtab(mlist)
}
