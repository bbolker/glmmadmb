admbControl <- function(impSamp=0,
                        maxfn=500,
                        imaxfn=500,
                        maxph=5,
                        noinit=TRUE,
                        shess=TRUE,
                        run=TRUE,
                        ZI_kluge=FALSE,
                        poiss_prob_bound=TRUE) {
  list(impSamp=impSamp,maxfn=maxfn,imaxfn=imaxfn,noinit=noinit,shess=shess,
       ZI_kluge=ZI_kluge,run=run,maxph=maxph,poiss_prob_bound=poiss_prob_bound)
  ## FIXME: do something clever with formals/match.call() ?
}

## FIXME: var-cov parameter names
## FIXME: correlation terms
mcmc_transform <- function(m,phi,fnames) {
    ## zero-inflation
    t_pz <- NULL
    if ("pz" %in% colnames(m)) {
        pz <- m[,"pz",drop=FALSE]
        t_pz <- pz  ## (not transformed)
    }
    ## fixed effects
    fixed <- m[,grep("^beta",colnames(m)),drop=FALSE]
    t_fixed <- fixed %*% phi
    colnames(t_fixed) <- fnames
    ## variance parameters: log std dev
    theta <- m[,grep("^tmpL",colnames(m)),drop=FALSE]
    t_theta <- exp(theta)
    ## corr parameters ("offdiagonal elements of cholesky-factor of correlation matrix")
    t_corr <- NULL
    if (any(corrcols <- grepl("^tmpL1",colnames(m)))) {
        ## warning("correlation parameters in MCMC not transformed")
        corr <- m[,corrcols,drop=FALSE]
        t_corr <- corr
    }
    ## scale/overdispersion parameter
    t_alpha <- NULL
    if ("log_alpha" %in% colnames(m)) {
        logalpha <- m[,grep("^log_alpha",colnames(m)),drop=FALSE]
        t_alpha <- matrix(exp(logalpha),dimnames=list(NULL,"alpha"))
    }
    ## random effects
    ## FIXME: need to scale u by appropriate sd (sigh)
    t_re <- NULL
    if (any(recols <- grepl("^u\\.[0-9]+",colnames(m)))) {
        ## warning("random effects in MCMC not transformed")
        re <- m[,recols,drop=FALSE]
        t_re <- re
    }
    m2 <- cbind(t_pz,t_fixed,t_theta,t_corr,t_alpha,t_re) 
    m2 <- mcmc(m2)
    ## FIXME: set start, end, thin explicitly?
    m2
}

mcmcControl <- function(mcmc=1000,
                         mcmc2=0,
                         mcsave,
                         mcnoscale=FALSE,
                         mcgrope=FALSE,
                         mcmult=1) {
  if (missing(mcsave)) mcsave <- pmax(1,floor(mcmc/1000))
  if (mcmc>0 && mcmc2>0) stop("may not specify both mcmc and mcmc2>0")
  r <- list(mcsave=mcsave,mcnoscale=mcnoscale,mcgrope=mcgrope,mcmult=mcmult)
  if (mcmc>0) c(list(mcmc=mcmc),r) else c(list(mcmc2=mcmc2,r))
}

mcmcArgs <- function(L) {
  argstr <- mapply(function(n,val) {
    if (is.numeric(val)) {
      paste("-",n," ",val,sep="")
    } else {
      if (isTRUE(val)) paste("-",val,sep="")
    }
  },names(L),L)
  paste(unlist(argstr),collapse=" ")
}

## FIXME: use .Deprecated?
glmm.admb <- function(...) {
  warning("'glmm.admb' has been renamed to 'glmmadmb' in the new version; please modify your code")
  glmmadmb(...)
}

get_bin_version <- function(file_name="glmmadmb") {
  res <- get_bin_loc(file_name)
  platform <- res$platform
  bin_loc <- res$bin_loc
  r <- run_bin(platform,bin_loc,file_name,cmdoptions="-version")
  message(r,sep="\n")
  invisible(r)
}

get_bin_loc <- function(file_name="glmmadmb",debug=FALSE) {
   nbits <- 8 * .Machine$sizeof.pointer
   if (.Platform$OS.type == "windows") {
     platform <- "windows"
   } else {
     ## MacOS detection OK?
     ## http://tolstoy.newcastle.edu.au/R/e2/help/07/01/8497.html
     if (substr(R.version$os,1,6)=="darwin") {
       platform <- "macos"
       unameres <- system("uname -v",intern=TRUE)
       if (grepl("Version 9",unameres)) {
         stop("glmmADMB binaries are not available for Mac OS 10.5 (Leopard). ",
              "Please see http://glmmadmb.r-forge.r-project.org for other options")
       }
     } else {
       if (R.version$os=="linux-gnu") {
         platform <- "linux"
       } else {
         stop("glmmadmb binary not available for OS ",R.version$os)
         ## FIXME: allow user to supply their own binary?
       }
     }
   }
   if (debug) cat("platform:",platform,nbits,"\n")
   execname <- if (platform=="windows") paste(file_name,"exe",sep=".") else file_name
   if (debug) cat("executable name:",execname,"\n")
   bin_loc <- system.file("bin",paste(platform,nbits,sep=""),
                          execname,
                          package="glmmADMB")
   if (debug) cat("bin_loc:",bin_loc,"\n")
   if (nchar(bin_loc)==0) stop(
              sprintf("glmmadmb binary should be available, but isn't (%s, %d bits) ",platform,nbits))
   list(bin_loc=bin_loc,platform=platform)
 }

run_bin <- function(platform,bin_loc,file_name,
                    cmdoptions,run=TRUE,rm_binary=TRUE,
                    debug=FALSE,verbose=FALSE,copy.exec=NULL) {
    ## copy executable even if not running code (i.e. make complete copy needed
    ##   to run ADMB outside of R)
    ## FIXME: for what platforms do we really need to copy the binary?
    ##  can't we just run it in place?
    ## Or does it do something silly and produce
    ##  output in the directory in which the binary lives rather than the
    ##  current working directory (feel like I struggled with this earlier
    ##  but have now forgotten -- ADMB mailing list archives??)
    if (is.null(copy.exec)) {
        copy.exec <- (platform!="windows")
    }
  if (!copy.exec) {
      cmd <- paste("\"",bin_loc, "\"", " ", cmdoptions, sep="")
  } else {
      pref <- if (platform!="windows") "./" else ""
      cmd <- paste0(pref, file_name, " ", cmdoptions)
      file.copy(bin_loc,".")
      ## file.copy seems to strip executable permissions ...
      Sys.chmod(file_name,mode="0755") 
  }
  if (debug) cat("Command line:",cmd,"\n")
  if (run) {
    if (platform=="windows") {
      sys.result <- shell(cmd, invisible=TRUE,intern=!verbose)
    } else  {
      sys.result <- system(cmd,intern=!verbose)
      if (rm_binary) unlink(file_name)
    }
  } else sys.result <- NULL
  sys.result
}

glmmadmb <- function(formula, data, family="poisson", link,start,
                     random,
                     corStruct="diag",  easyFlag=TRUE,
                     zeroInflation=FALSE,
                     admb.opts=admbControl(),
                     mcmc=FALSE,
                     mcmc.opts=mcmcControl(),
                     save.dir, verbose=FALSE,
                     extra.args="",
                     bin_loc=NULL,
                     copy.exec=TRUE,
                     debug=FALSE)
{

  file_name <- "glmmadmb"

  if (!missing(easyFlag)) warning("easyFlag argument ignored")

  ## find platform, binary location information
  res <- get_bin_loc(file_name,debug=debug)
  platform <- res$platform
  if (is.null(bin_loc)) bin_loc <- res$bin_loc

  ## extract admb.opts contents
  impSamp <- admb.opts$impSamp
  maxfn <- admb.opts$maxfn
  imaxfn <- admb.opts$imaxfn
  run <- admb.opts$run
  
  if (use_tmp_dir <- missing(save.dir)) {
    repeat {
      save.dir <- paste(tempfile(pattern="glmmADMB"))
      if (!file.exists(save.dir)) {
        if (debug) cat("using temp directory",save.dir,"\n")
        break
      }
    }
  }
  if (debug) cat("saving files into directory",save.dir,"\n")
  if (newdir <- !file_test("-d",save.dir)) {
    if (debug) cat("creating temp directory\n")
    dir.create(save.dir)
  }
  owd <- setwd(save.dir)
  if (debug) cat("changed working directory to",getwd(),"\n")
  on.exit({
    setwd(owd)
    if (debug) cat("restored working directory to",getwd(),"\n")
  })
  if (use_tmp_dir) {
    on.exit({
      unlink(save.dir,recursive=TRUE)
      if (debug) cat("removed temp directory",save.dir,"\n")
    },
            add=TRUE)
  }

  zi_model_flag <- alpha_model_flag <- FALSE
  if (is.list(formula)) {
      if (!is.null(zi_model <- formula[["zi"]])) zi_model_flag <- TRUE
      if (!is.null(alpha_model <- formula[["alpha"]])) alpha_model_flag <- TRUE
      formula <- if (names(formula)[1]=="") formula[[1]] else formula[["eta"]]
  }

  if (zi_model_flag || alpha_model_flag)
      stop("models for zero-inflation and dispersion parameters not yet implemented")
  
  call <- match.call()
  has_rand <- !missing(random) || length(grep("\\|",as.character(formula)[3]))>0
  ## FIXME: we can no longer tell easily if impSamp was specified
  ##  (because it is specified via admb.control())
  if (!has_rand && (!missing(corStruct))) 
      stop("No random effects specified: \"corStruct\" does not make sense")

  if (!is.character(family)) stop("must specify 'family' as a character string")
  family <- tolower(family)
  family <- gsub("binomial","binom",family)
  nbinom1_flag <- 0
  if (family %in% c("nbinom1","truncnbinom1") ) {
    nbinom1_flag <- 1
  }
  ## [trunc]nbinom2 is a synonym for [trunc]nbinom
  if (length(grep("nbinom",family))>0)
    family <- gsub("[12]$","",family)

  like_type_flag <- switch(family,poisson=0,binom=1,nbinom=2,gamma=3,beta=4,gaussian=5,
                           truncpoiss=6,truncnbinom=7,logistic=8,betabinom=9,
                           stop("unknown family"))
  has_alpha <- family %in% c("nbinom","gamma","beta","gaussian","truncnbinom","logistic",
                             "betabinom")

  if (missing(link)) {
      link <- switch(family,
                     binom=, beta=, betabinom="logit",
                     nbinom=, poisson=, truncpoiss=, truncnbinom=, gamma="log",
                     gaussian=, logistic="identity")
  }
  linkfun <- switch(link,log=log,logit=qlogis,probit=qnorm,inverse=function(x) {1/x},
                    cloglog=function(x) {log(-log(1-x))},
                    identity=identity,
                    stop("unknown link function"))
  ilinkfun <- switch(link,log=exp,logit=plogis,probit=pnorm, inverse=function(x) {1/x},
                     cloglog=function(x) {
                       pmax(pmin(-expm1(-exp(x)),
                                 1 - .Machine$double.eps), .Machine$double.eps)
                       },
                       identity=identity
                     )

  link_type_flag <- switch(link,log=0,logit=1,probit=2,inverse=3,cloglog=4,identity=5)
  
  ## from glm()
  ## extract x, y, etc from the model formula and frame
  if (missing(data)) 
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  
  m <- match(c("formula", "data", "subset", "na.action", 
               "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  ## FIXME/WARNING!!! added parent.frame().  Does this always work???
  mf$formula <- get_fixedformula(eval(mf$formula,parent.frame()))

  comb_formula <- function(old,new) {
      if (is.null(new)) return(old)
      update.formula(old,paste0(".~.+",as.character(new)[2]))
  }

  ## FIXME: use Reduce()?
  if (alpha_model_flag) {
      mf$formula <- comb_formula(mf$formula,alpha_model)
  }
  if (zi_model_flag) {
      mf$formula <- comb_formula(mf$formula,zi_model)
  }

  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")

  ## FIXME: warning on evaluation with nested factors?
  mf_orig <- mf
  mf <- eval(mf, parent.frame())

  mt <- attr(mf, "terms") # allow model.frame to have updated it

  if (is.data.frame(data) && nrow(mf)<nrow(data))
    warning("NAs removed in constructing fixed-effect model frame: ",
            "you should probably remove them manually, e.g. with na.omit()")
  
  y <- model.response(mf, "any") # e.g. factors are allowed
  if (inherits(y,"factor")) {
    if (family!="binom") stop("factors as response variables only allowed for binomial models")
    ## slightly odd, but this is how glm() defines it: first level is failure, all others success
    y <- 1-as.numeric(as.numeric(y)==1)
  }
  y <- as.matrix(y)
  ## FIXME: add tolerance to this test ...
  if (family %in% c("nbinom","truncnbinom","poisson","truncpoiss","binom","betabinom") &&
      any(y!=floor(y))) {
    warning("non-integer response values in discrete family")
  }
  if (substr(family,1,5)=="trunc" && any(y==0)) {
    warning("zero response values in truncated family")
  }
  
  n <- nrow(y)
  p_y <- ncol(y)

  nyobs <- if (family %in% c("binom","betabinom")) rowSums(y) else 1

  ## BMB: nobs() method? model.matrix, model.frame, terms methods?

  offset <- as.vector(model.offset(mf))
  has_offset <- !is.null(offset)
  if (!has_offset) offset <- rep(0,n)

  X <- model.matrix(mt,mf)  ## BMB: allow contrasts?
  p <- ncol(X)

  if ((rankX <- rankMatrix(X))<p)
      stop(gettextf("rank of X = %d < ncol(X) = %d", rankX, p))

  zi_p <- 1
  zi_X <- matrix(0,1,1)
  if (zi_model_flag) {
      zi_X <- model.matrix(zi_model,mf)
      zi_p <- ncol(zi_X)
  }
  
  alpha_p <- 1
  alpha_X <- matrix(0,1,1)
  if (alpha_model_flag) {
      alpha_X <- model.matrix(alpha_model,mf)
      alpha_p <- ncol(alpha_X)
  }

  if (has_rand) {
    REmat <- process_randformula(formula,random,data=data)

    ## kluge: replicate elements
    if (any(REmat$nterms>1)) {
      tmpind <- rep(seq_along(REmat$mmats),REmat$nterms)
      REmat$mmats <- REmat$mmats[tmpind]
      ## fix names: ???
      names(REmat$mmats) <- names(REmat$levels)
      REmat$codes <-
        unlist(lapply(REmat$codes,
                      function(x) lapply(as.list(as.data.frame(x)),matrix)),
               recursive=FALSE)
      ## split into a list of single matrices
    }
    
    ## Number of random effects within each crossed term
    m <- sapply(REmat$mmats, ncol)

    ## Number of crossed terms
    M <- length(m)

    ## FIXME: allow variable corStructs
    ## if (length(corStruct)>1) {
    ## stop("varying corStructs for different random effects not yet implemented")
    ## }
    corStruct <- rep(corStruct,length.out=M)
    ## if (!all(corStruct %in% c("diag","full"))) stop("corStruct must be either 'diag' or 'full'")

    m2 <- ifelse(corStruct=="diag",m,m*(m+1)/2)

    ## Number of levels of each crossed term
    q <- sapply(REmat$levels,length)

    ## Construct design matrix (Z) for random effects
    Z <- do.call(cbind,REmat$mmats)
    colnames(Z) <- paste(rep(names(REmat$mmat),m),
                         unlist(sapply(REmat$mmat,colnames)),sep=".")

    ## Make matrix of pointers into ADMB vector of random effects "u".
    ## Rows in II contain pointers for a given data point 
    ## First offset the columns of II to point into different segments of u
    II <- matrix(rep(q,m),nrow=n,ncol=sum(m),byrow=TRUE)
    if(sum(m)>1)
    for(i in 2:sum(m))
      II[,i] = II[,i] + II[,i-1]
    II <- cbind(0,II[,-ncol(II)])
    ##    tmpf <- function(x) c(0,cumsum(x)[-length(x)])
    ##    II <- t(apply(II,1,tmpf))
    ii <- 1
    for(i in 1:M)		# Fill in the actual random effects codes
      for(j in 1:m[i])
        {
          II[,ii] = II[,ii] + REmat$codes[[i]]
          ii = ii+1
        }
    
    ## Splits u into "correlation blocks" (no way yet of specifying uncorrelated random effects)
    cor_block_start = cumsum(m)-m[1]+1
    cor_block_stop = cumsum(m)
    numb_cor_params = max(sum(m*(m-1)/2),1)	# For technical reasons we need at least 1 parameter: see glmmadmb.tpl
    
  } else {
    ## no random factor: fill in with dummies
    m <- M <- q <- 1
    Z <- matrix(0,ncol=1,nrow=n)
    II <- matrix(rep(q,m),nrow=n,ncol=sum(m),byrow=TRUE)
    cor_block_start <- cor_block_stop <- 1
    numb_cor_params <- 1
  }

  
  cmdoptions <- paste("-maxfn",maxfn)
  if (!is.null(admb.opts$maxph) && !is.na(admb.opts$maxph)) {
      cmdoptions <- paste(cmdoptions,"-maxph",admb.opts$maxph)
  }
  if (admb.opts$noinit) cmdoptions <- paste(cmdoptions,"-noinit")
  if (admb.opts$shess) cmdoptions <- paste(cmdoptions,"-shess")
  if (has_rand && impSamp>0) cmdoptions <- paste(cmdoptions,"-is",impSamp)
  if (mcmc) cmdoptions <- paste(cmdoptions,mcmcArgs(mcmc.opts))
  if (!missing(extra.args)) cmdoptions <- paste(cmdoptions,extra.args)
  cor_flag <- as.numeric(corStruct=="full")
  dat_list <- list(n=n, p_y=p_y,
                   y=y, p=p, X=X, M=M, q=q, m=m, ncolZ=ncol(Z),Z=Z, II=II, 
                   cor_flag=cor_flag,
                   cor_block_start=cor_block_start,
                   cor_block_stop=cor_block_stop,
                   numb_cor_params=numb_cor_params,
                   like_type_flag=like_type_flag,
                   link_type_flag=link_type_flag,
                   rlinkflag=1L,   ## always robust (for now)
                   ## as.numeric(family=="poisson"||(family=="binomial"&&link=="logit")),
                   no_rand_flag=as.numeric(!has_rand),
                   zi_flag=as.numeric(zeroInflation),
                   zi_kluge=as.numeric(admb.opts$ZI_kluge),
                   poiss_prob_bound=as.numeric(admb.opts$poiss_prob_bound),
                   nbinom1_flag=as.numeric(nbinom1_flag),
                   intermediate_maxfn=imaxfn, 
                   has_offset=as.numeric(has_offset), 
                   offset=offset)

  ## TO DO: implement models for zero-inflation and dispersion parameters
  ## zi_model_flag=as.numeric(zi_model_flag),
  ## zi_X=zi_X,
  ## zi_p=zi_p,
  ## alpha_model_flag=as.numeric(alpha_model_flag),
  ## alpha_X=alpha_X,
  ## alpha_p=alpha_p)


  ## FIXME: pz=0.0001 should be clearly specified,
  ##    possibly made into a control parameter
  pz_start <- c(if(zeroInflation) 0.02 else 0.0001,rep(0,zi_p-1))
  log_alpha_start <- c(1.0,rep(0,alpha_p-1))
  pin_list <- list(pz=pz_start,
                   b=numeric(p),                        ## fixed effects
                   tmpL=0.25+numeric(sum(m)),           ## log-std dev of RE
                   tmpL1=1e-4+numeric(numb_cor_params), ## off-diag of cholesky factor of corr matrix
                   log_alpha=log_alpha_start,           ## dispersion param
                   u=rep(0,sum(m*q)))                   ## random effects/conditional modes

  if (!missing(start) && !is.null(start$fixed)) {
      ## message("converting fixed to orthogonalized variant ...")
      phi <- make_phi(X)
      start$fixed <- start$fixed %*% solve(phi)
  }

  rnames <- list(c("fixed","b"),
                 c("RE_sd","tmpL"),
                 c("RE_cor","tmpL1"),
                 c("pz","pz"),
                 c("log_alpha","log_alpha"))
  if (!missing(start)) {
      for (i in seq_along(rnames)) {
          names(start)[names(start)==rnames[[i]][1]] <- rnames[[i]][2]
      }
      ns <- names(start)
      for (i in seq_along(start)) {
          pp <- pin_list[[ns[i]]]
          if (is.null(pp)) {
              stop("unmatched start component ",ns[i])
          }
          if (length(pp) != length(start[[i]])) {
              stop("length mismatch in start component ",ns[i],
                   ": ",length(pp),"!=",length(start[[i]]))
          }
          pin_list[[ns[i]]] <- start[[i]]
      }
  }

  if (debug) {
      cat("writing .dat and .pin files\n")
      cat("working directory:",getwd(),"\n")
  }

  dat_write(file_name, dat_list)
  pin_write(file_name, pin_list)
  par_file <- paste0(file_name, ".par")
  std_file <- paste0(file_name, ".std")
  if(file.exists(std_file) && run) warning("file ",std_file," exists: overwriting")
  ## file.remove(std_file)

  sys.result <- run_bin(platform,bin_loc,file_name,cmdoptions,run,
                        rm_binary=!use_tmp_dir,debug=debug,
                        copy.exec=copy.exec,verbose=verbose)
  
  if (!file.exists(std_file)) {
    if (run) {
        if (file.exists(par_file)) {
            message("Parameters were estimated, but standard errors were not: ",
                    "the most likely problem is that the curvature at MLE was zero or negative")
            ##
            ## would like to read parameter files, but not exactly in standard format ...
            ## cat("Parameter file:\n")
            ## file.show(par_file)
        }
        if (debug) cat("run failed:",sys.result,"\n")
        stop("The function maximizer failed (couldn't find parameter file)",
             " Troubleshooting steps include (1) run with 'save.dir' set ",
             "and inspect output files; (2) change run parameters: see '?admbControl';",
             "(3) re-run with debug=TRUE for more information on failure mode")
    }
    message("'run=FALSE' specified, STD file not found: stopping")
    return(NULL)
  }
  ## parameter order:
  ## pz (1)
  ## beta      ## nfixpar
  ## real_beta ## nfixpar
  ## tmpL (log standard dev)   ## nvar
  ## tmplL1 (offdiag elements) ## ncor
  ## log_alpha                 ## alpha
  ##
  tmp <- read.table(paste(file_name,"std",sep="."), skip=1,as.is=TRUE)
  ## FIXME: could we change the TPL file to write everything out in full precision??
  ##   ... otherwise to read .par file or binary versions ...
  ## BMB: if 'std dev' were written without a space we could use header=TRUE
  tmpindex <- tmp[,2]
  out <- list(n=n, q=q, formula=formula,
              fixed = get_fixedformula(formula),
              call=call,
              frame=mf,
              terms=mt,
              family=family, corStruct=corStruct, impSamp=impSamp,
              easyFlag=easyFlag, zeroInflation=zeroInflation)
  if(zeroInflation) {
    out$pz <- as.numeric(tmp[tmpindex=="pz", 3])
    out$sd_pz <- as.numeric(tmp[tmpindex=="pz", 4])
  }
  out$b <- as.numeric(tmp[tmpindex=="real_beta", 3])
  out$stdbeta <- as.numeric(tmp[tmpindex=="real_beta", 4])
  names(out$stdbeta) <- names(out$b) <- colnames(X)

  if (has_alpha)
    {
      out$alpha <- as.numeric(tmp[tmpindex=="alpha", 3])
      out$sd_alpha <- as.numeric(tmp[tmpindex=="alpha", 4])
    } 


  out$offset <- offset
  
  ## if(!missing(link))
  out$link <- link
  out$linkfun <- linkfun
  out$ilinkfun <- ilinkfun
  
  if(has_rand) {
      ## BMB: fixme! make sure this works for multiple random effects
      ## reconstruct variance-covariance matrices
      Svec <- tmp[tmpindex=="S",3]  ## lower triangle??
      parnames <- lapply(REmat$mmat,colnames)
      groupnames <- names(REmat$mmats)
      dn <- mapply(list,
                   parnames,
                   parnames,SIMPLIFY=FALSE)
      ##
      ## construct matrix from diagonal or lower triangle
      mkmat <- function(x,n,dimnames) {
          if (length(x)==1) { ## single element (circumvent diag(x) weirdness)
              m <- matrix(x)
          } else if (length(x)==n) {
              m <- diag(x)
          } else {
              m <- matrix(NA,nrow=n,ncol=n)
              m[lower.tri(m,diag=TRUE)] <- x
              m[upper.tri(m)] <- t(m)[upper.tri(m)] ## symmetrize
          }
          dimnames(m) <- dimnames
          m
      }
      out$S <- mapply(mkmat,split(Svec,rep(1:length(m),m2)),
                      n=m,
                      dimnames=dn,
                      SIMPLIFY=FALSE)
      sd_Svec <- tmp[tmpindex=="S",4]
      out$sd_S <- mapply(mkmat,split(sd_Svec,rep(1:length(m),m2)),
                         n=m,
                         dimnames=dn,
                         SIMPLIFY=FALSE)
      names(out$S) <- names(out$sd_S) <- groupnames
      ## if(corStruct == "diag") {
      ##     replace_offdiag <- function(m,r=0) {
      ##         m[upper.tri(m) | lower.tri(m)] <- r
      ##         m
      ##     }
      ##     out$S <- lapply(out$S,replace_offdiag)
      ##     out$sd_S <- lapply(out$sd_S,replace_offdiag)
      ##     ## FIXME: replace with NA for sd_S?
      ## }
      ## out$random <- random
      ## FIXME: check dimnames for a wider range of cases ...
      uvec <- tmp[tmpindex=="u",3]
      ulist <- split(uvec,rep(1:length(q),q*m))
      dn <- mapply(list,
                   REmat$levels, ## lapply(REmat$mmat,attr,which="levels"),
                   parnames,
                   SIMPLIFY=FALSE)
      out$U <- mapply(matrix,
                      ulist,
                      ncol=m,
                      nrow=q,
                      dimnames=dn,
                      SIMPLIFY=FALSE)
      names(out$U) <- groupnames
      ## BMB/FIXME: should this be byrow or not? check ...
      ii <- 1
      allU <- matrix(nrow=n,ncol=sum(m))
      for (i in 1:length(m)) {
          allU[,ii:(ii+m[i]-1)] <- out$U[[i]][c(REmat$codes[[i]]),]
          ii <- ii+m[i]
      }
    ## U <- matrix(as.numeric(tmp[tmpindex=="u",3]), ncol=m, byrow=TRUE,
    ## ## dimnames=list(data[,group][group_d[-1]-1],colnames(Z)))
    ## dimnames=list(NULL,colnames(Z)))
    sd_uvec <- tmp[tmpindex=="u",4]
    sd_ulist <- split(sd_uvec,rep(1:length(q),q*m))
    out$sd_U <- mapply(matrix,
                       sd_ulist,
                       ncol=m,
                       nrow=q,
                       dimnames=dn,
                       SIMPLIFY=FALSE)
    ## FIXME: check byrow = TRUE or FALSE
    ## matrix(as.numeric(tmp[tmpindex=="u",4]), ncol=m, byrow=TRUE,
    ## dimnames=list(data[,group][group_d[-1]-1],colnames(Z)))
    ##    dimnames=list(NULL,colnames(Z)))
  } else {
    out$U <- out$sd_U <- matrix(rep(0,q), ncol=1, byrow=TRUE)
  } ## !has_rand

  mu <- as.numeric(X %*% out$b) + offset

  if (has_rand) {
      sdvals <- unlist(lapply(out$S,function(z)sqrt(diag(z))))
      eta <- mu + rowSums(Z * sweep(allU,2,sdvals,"*"))
  } else eta <- mu
  
  lambda <- ilinkfun(eta)

  out$fitted <- lambda
  if(zeroInflation) {
      out$fitted <- out$fitted * (1-out$pz)
  }

  out$sd.est <- with(out,switch(family,
                                poisson=sqrt(lambda),
                                nbinom=sqrt(lambda*(1+lambda/alpha)),
                                nbinom1=sqrt(lambda*alpha),
                                betabinom=, binom=sqrt(fitted*(1-fitted)/nyobs),
                                beta=sqrt(lambda*(1-lambda)/(1+alpha)),
                                ## gamma: parameterized as scale, mean (=shape*scale)
                                gamma=sqrt(lambda*alpha),
                                    {warning("sd.est not defined for this family"); rep(NA,length(lambda))}))
  ##  stop("sd.est not defined for family",family))

  if (family %in% c("binom","betabinom") && p_y>1) {
      out$residuals <- y[,1]/nyobs-lambda
  } else  out$residuals <- y-lambda
  tmp <- par_read(file_name)
  out$npar <- tmp$npar   ## BMB: should this be total number of parameters or number of fixed parameters?
  bpar <- bar_read(file_name,n=tmp$npar+1)[-1] ## drop ZI parameter (first)
  if (file.exists("phi.rep")) {
      out$phi <- matrix(scan("phi.rep",quiet=TRUE),nrow=length(out$b),byrow=TRUE)
      ## tmp$beta <- bpar
      newb <- bpar[seq_along(out$b)] %*% out$phi
      ## FIXME: find a way to check consistency between .bar and .par that isn't subject
      ##        to false positives (i.e., replicate the rounding that ADMB does?)
      ## if (!isTRUE(all.equal(c(newb),unname(out$b),tol=5e-5))) {
      ## if (!all(round(newb,3)==round(out$b,3))) {
      ## warning("apparent discrepancy between .bar and .par values: using .par values")
      out$b[] <- newb ## replace values
  }  ## FIXME: should figure out why phi.rep is not getting created on MacOS -- binary created from old TPL?
  out$loglik <- tmp$loglik
  out$gradloglik <- tmp$gradient
  nfixpar <- length(out$b)
  ## drop cors that don't correspond to fixed-effect parameters
  out$corMat <- tmp$cor[1:nfixpar,1:nfixpar]

  out$conv <- 0
  out$convmsg <- ""
  
  if (abs(out$gradloglik) >= 0.001) {
    out$convmsg <- paste("log-likelihood of gradient=",out$gradloglik)
    out$conv <- 1
    warning("Convergence failed:",out$convmsg)
  }

  ## BMB: warning/convergence code for bad hessian?

  ## FIXME: figure out names of parameters (other than fixed)
  if (mcmc) {
      out$mcmc <- as.matrix(R2admb::read_psv(file_name))
      colnames(out$mcmc) <- R2admb:::rep_pars(tmpindex)[1:ncol(out$mcmc)]
      out$mcmc <- mcmc_transform(out$mcmc,out$phi,fnames=names(out$b))
  }

  if (has_rand) out$random <- REmat$rstring
  
  class(out) <- "glmmadmb"

  return(out)
}  
