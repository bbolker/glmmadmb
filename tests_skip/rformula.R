
## code for extracting random effects
## note: for this to work we shouldn't have anything funky
##  with nested parentheses.  Also don't know what will
##  happen if fixed effects are strewn between random effects:

## multiple terms of the form
##  (<expr1>|<expr2>)
## where <expr1> describes the covariates affected by/interacting with
## the random effect and <expr2> describes the structure of the random effect

process_randformula <- function(f,data) {
  rchar <- as.character(f[2]) ## RHS
  ## drop bits before first/after last parenthesis
  rchar <- gsub("^[^()]*\\(","",
                gsub(")[^()]*$","",rchar))
  randbits <- grep("\\|",strsplit(rchar,"[()]")[[1]],value=TRUE)
  splitbits <- strsplit(randbits,"\\|")
                              
  cfun <- function(lbit,mdata) {
    model.matrix(as.formula(paste("~",lbit)),mdata)
  }

  ## here we want to expand RHS and provide a list of indices
  ## into the appropriate factor:
  rfun <- function(rbit,rdata) {
    f <- as.formula(paste("~",rbit,"-1"))
    t <- terms(f,data=rdata)
    ## ugly: "If the answer is parse() you should usually rethink the question" but ??
    labs <- attr(t,"term.labels")
    sapply(labs,
           function(lab) as.numeric(with(data,eval(parse(text=lab)))))
  }
  termnames <- gsub("\\|","_bar_",
                    gsub(":","_int_",
                         gsub("/","_nest_",
                              gsub("\\*","_cross",
                                   gsub(" ","",randbits)))))
  L <- list(mmats=lapply(lapply(splitbits,"[",1),cfun,mdata=data),
       codes=lapply(lapply(splitbits,"[",2),rfun,rdata=data))
  names(L$mmats) <- names(L$codes) <- termnames
  L
}

## print tp file (not necessarily 
write_randformula <- function(x,name) {
  require(glmmADMB)
  fn <- if(substring(name,nchar(name)-3)==".dat") {
    name
  } else paste(name,".dat",sep="")
  cat("### design matrices for random effects:\n",file=fn)
  R2admb::dat_write(name,x$mmat,append=TRUE)
  cat("### indices for random effects:\n",file=fn,append=TRUE)
  R2admb::dat_write(name,x$codes,append=TRUE)
}

set.seed(1001)
tdat <- expand.grid(f1=LETTERS[1:5],f2=letters[1:5],rep=1:3)
tdat$x <- runif(nrow(tdat))

p1 <- process_randformula(~ x + (1|f1)+(x|f2),data=tdat)
p2 <- process_randformula(~ x + (1|f1+f2),data=tdat)
p3 <- process_randformula(~ x + (1|f1)+(1|f2),data=tdat)


write_randformula(p1,"test")
