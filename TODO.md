# Wish/to-do list for glmmADMB

* formula interface for zero-inflation
* more complex variance models (ASREML syntax)
* Can we add some sort of progress bar or "working ..." message in the absence of voluminous messages (difficult); better isolation of messages to set command-line flags
* Add more examples, tests, vignettes ... especially a smaller (quicker) data set for the main example
* work on predict() method; predictions incl REs
* simulate() method
* save (more) fitted models in inst/ so they can be used easily in examples without doing the fits each time
* speed up (Skaug)
* more robust first pass, weighted SS or quasi-likelihood (Fournier)
* access to profiling
* MCMC fixes: transform correlation parameter
* sort out `terms`, `formula` methods to make `step` etc. work properly (check `lme4` implementation, add `fixed.only=FALSE` option?)
* would making predictions into sdreport variables give (Wald) CIs that included uncertainty in Var-Corr?
* vignette: more on troubleshooting, R2admb intro, importance sample, AGQ, ... ? 




