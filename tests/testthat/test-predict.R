library(glmmADMB)
library(testthat)


test_that("predict with factors", {
    set.seed(101)
    dd <- data.frame(y=rnorm(20,mean=rep(1:2,each=10)),
                     f=factor(rep(c("a","b"),each=10)))
    m1 <- glmmadmb(y~f,data=dd,family="gaussian")
    expect_equal(predict(m1,newdata=data.frame(f=c("a","b"))),
                 rev(predict(m1,newdata=data.frame(f=c("b","a")))))
})
