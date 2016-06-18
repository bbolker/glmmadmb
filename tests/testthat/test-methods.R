library("testthat")
library("glmmADMB")
L <- load(system.file("testdata", "glmmadmb-tst-fits.rda",
                      package="glmmADMB", mustWork=TRUE))

## FIXME: should test for old R versions, skip reloading test data in that
## case?
fm0 <- fit_sleepstudy_0
fm1 <- fit_sleepstudy_1
fm2 <- fit_sleepstudy_2
gm1 <- fit_cbpp_1
gm2 <- fit_cbpp_2
gm3 <- fit_cbpp_3

test_that("drop1", {
    })
