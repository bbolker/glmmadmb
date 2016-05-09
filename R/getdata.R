getdata <- function(s) {
    ## FIXME: if s is missing, generate list of available files
    ## FIXME: be more flexible about extensions
    ## FIXME: can data() be adapted for this??
    p <- parent.frame()
    load(system.file("extdata",paste(s,".rda",sep=""),
                    package="glmmADMB"),envir=p)
}

