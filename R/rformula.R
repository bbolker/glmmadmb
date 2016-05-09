
## code for extracting random effects
## note: for this to work we shouldn't have anything funky
##  with nested parentheses.  Also don't know what will
##  happen if fixed effects are strewn between random effects:

## multiple terms of the form
##  (<expr1>|<expr2>)
## where <expr1> describes the covariates affected by/interacting with
## the random effect and <expr2> describes the structure of the random effect

get_fixedformula <- function(f) {
  lchar <- as.character(f[2])
  rchar <- as.character(f[3]) ## RHS
  offsetstr <- ""
  has_offset <- function(x) {grepl("offset\\(",x)}
  if (has_offset(rchar)) {
    ## protect/remove offset
    off1str <- "offset\\([^()]*\\)"
    off2str <- "offset\\([^()]*\\([^()]*\\)[^()]*\\)"
    ## try nested parens
    offsetstr <- gsub(paste(".*(\\+ *",off2str,").*",sep=""),"\\1",rchar)
    rchar <- gsub(off2str,"",rchar)
    if (has_offset(rchar)) {
      offsetstr <- gsub(paste(".*(\\+ *",off1str,").*",sep=""),"\\1",rchar)
      rchar <- gsub(off1str,"",rchar)
    }    
    if (has_offset(rchar)) stop("unable to process offset term")
  }
  rchar <- gsub("\\([^)|]+\\|[^)|]+\\)","",rchar) ## parentheses containing |
  rchar <- gsub("(\\+ *\\+ *)+","+",rchar) ## duplicated +
  rchar <- gsub(" *\\++ *$","",rchar) ## terminating + (possibly multiple)
  if (rchar=="") rchar <- "1"  ## empty fixed formula
  as.formula(paste(lchar,"~",rchar,offsetstr))
}

## test strings
## get_fixedformula(y~Base*trt+Age+Visit+(Visit|subject))
## get_fixedformula(y~Base*trt+Age+Visit+poly(a,b,c)+(Visit|subject))
## get_fixedformula(y~Base*trt+Age+Visit+poly(a,b,c)+(Visit|subject)+(1|zzz))

### needs to be outside so can be seen even if data= not specified
Droplevels <- function(x) {
  if (!is.factor(x)) stop("all grouping variables in random effects must be factors")
  droplevels(x)
}

process_randformula <- function(f,random,data) {
  rchar <- as.character(f[3]) ## RHS
  if (!missing(random)) {
    if (grepl("\\(.+\\|.+\\)",rchar))
      stop("must specify random effects *either* as part of 'formula' *or* in 'random'")
    rchar <- as.character(random[2])
  } else {
    ## drop bits before first/after last parenthesis
    rchar <- gsub("^[^()]*\\(","",
                  gsub(")[^()]*$","",rchar))
  }

  ## separate random components into a character vector
  randbits <- grep("\\|",strsplit(rchar,"[()]")[[1]],value=TRUE)
  

  ## FIXME: test whether this works with fixed effects embedded
  ##   between random effects??

  ## random components divided into structure and grouping variable
  ##  (before and after |)
  splitbits <- strsplit(randbits,"\\|")

  ## generate model matrix from component
  ## FIXME: fails with ~1 and no data ...
  cfun <- function(lbit,mdata) {
    ##
    if (!is.null(nrow(mdata))) {
      m <- model.matrix(as.formula(paste("~",lbit)),mdata)
    } else {
      ## hack: pull LHS of formula from global f
      m <- model.matrix(as.formula(paste(as.character(f[2]),"~",lbit)))
    }
    m
  }

  
  ## expand RHS and provide a list of indices
  ## into the appropriate factor:
  rfun <- function(rbit,rdata) {
    f <- as.formula(paste("~",rbit,"-1"))
    t <- terms(f,data=rdata)
    ## ugly: "If the answer is parse() you should usually rethink the question" but ??
    labs <- attr(t,"term.labels")
    sapply(labs,
    ##       function(lab) as.numeric(with(data,Droplevels(eval(parse(text=lab))))))
           function(lab) as.numeric(Droplevels(eval(parse(text=lab),data))))
  }
  
  termnames <- gsub("\\|","_bar_",
                    gsub(":","_int_",
                         gsub("/","_nest_",
                              gsub("\\*","_cross",
                                   gsub(" ","",randbits)))))
  groups <- gsub("^ +","",lapply(splitbits,"[",2))

  ## expand formulae

  ## terms list
  tL <- lapply(as.list(groups),
         function(x) {
           tt <- terms(formula(paste("~",x,sep="")),data=data)
           tt
         })

  ## list of ATOMIC factors
  groupL <- lapply(tL,
                   function(z) {
                     ## FAILED attempt to reduce variables by getting rid of unused levels:
                     ## x <- lapply(eval(attr(z,"variables"),data),Droplevels)
                     x <- eval(attr(z,"variables"),data)
                     names(x) <-rownames(attr(z,"factors"))
                     x
                   })

  nterms <- sapply(tL,function(x) length(attr(x,"term.labels")))
                 

  dgroupL <- lapply(tL,
                   function(z) {
                     ff <- as.list(colnames(attr(z,"factors")))
                     ## FAILED attempt to reduce variables by getting rid of unused levels:
                     lapply(ff,function(x) Droplevels(eval(parse(text=x),data)))
                     ## DANGER WILL ROBINSON: eval(parse(...)) !
                     ##lapply(ff,function(x) eval(parse(text=x),data))
                   })
  
  ## flat group list
  fgroupL <- unlist(groupL,recursive=FALSE)

  
  ## nonfactors <- groupL[!sapply(data[groups],inherits,"factor")]
  ## FIXME: want to be able to have valid formulas
  ##  (e.g. a/b/c) on RHS of | ...
  if (any(!sapply(fgroupL,inherits,"factor")))
    stop("all grouping variables must be factors")
  ## FIXME: identify which variables are non-factors

  ## n.b. this is LHS of the RE chunk, not of the overall formula
  LHS <- gsub("^ +","",lapply(splitbits,"[",1))

  ngroupfac <- sapply(groupL,length)
  ## if (any(ngroupfac>1)) {
  ## LHS <- rep(LHS,ngroupfac)
  ## }

  L <- list(mmats=lapply(LHS,cfun,mdata=data),
            codes=lapply(groups,rfun,rdata=data),
            levels=lapply(unlist(dgroupL,recursive=FALSE),levels),
            nterms=nterms,   ## num. derived terms per RE mat
            rstring=rchar) 
  ## for (i in seq_along(L$mmats)) {
  ##   attr(L$mmats[[i]],"levels") <- levels(fgroupL[[1]])
  ## }
  rnames <- unlist(lapply(L$codes,colnames)) 
  names(L$mmats)  <- groups
  names(L$levels) <- rnames
  names(L$codes) <- termnames
  ## ??? junk ???
  ## attr(names,"groupnames1") <- groups
  ## attr(names,"groupnames2") <- rnames
  ## attr(names,"levelnames") <- termnames
  L
}

## print to file
write_randformula <- function(x,name) {
  fn <- if(substring(name,nchar(name)-3)==".dat") {
    name
  } else paste(name,".dat",sep="")
  cat("### design matrices for random effects:\n",file=fn)
  R2admb::dat_write(name,x$mmat,append=TRUE)
  cat("### indices for random effects:\n",file=fn,append=TRUE)
  R2admb::dat_write(name,x$codes,append=TRUE)
}

