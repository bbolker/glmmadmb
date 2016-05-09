par_read <- function(name)
{
  filename <- if(substring(name,nchar(name)-3)==".par") name else paste(name,".par",sep="")
  ## filename <- if(tools::file_ext(name)=="par") name else paste(name,".par",sep="") # R 2.11 and later

  tmp <- scan(filename, what="", quiet=TRUE)
  tmp2 <- split(tmp, cumsum(tmp=="#"))
  x <- tmp2[-1]

  for(i in 1:length(x))
  {
    y <- x[[i]]
    n <- nchar(y[2])
    x[[i]] <- as.numeric(y[-(1:2)])
    names(x)[i] <- substring(y[2], 1, n-1)
  }

  x$loglik <- -as.numeric(tmp2[[1]][11])
  x$gradient <- -as.numeric(tmp2[[1]][16])
  x$npar <- as.numeric(tmp2[[1]][[6]])

  ## stuff stolen from R2ADMB
  rt <- function(f,ext,...) {
    fn <- paste(f,ext,sep=".")
    if (file.exists(fn)) read.table(fn,...) else NA
  }

  ## FIXME:: better R2admb integration
  ncorpar <- length(readLines(paste(name,"cor",sep=".")))-2
  cor_dat <- rt(name,"cor", skip = 2, fill=TRUE, 
                as.is=TRUE,col.names=paste("X",1:(4+ncorpar),sep=""))
  start_pos <- which(cor_dat[,2]=="real_beta")[1]-1
  cormat <- as.matrix(cor_dat[start_pos+(1:x$npar),start_pos+4+(1:x$npar)])
  cormat[upper.tri(cormat)] <- t(cormat)[upper.tri(cormat)]
  x$cormat <- cormat
  dimnames(x$cormat) <- list(names(x$b),names(x$b))
  
  return(x)
}

bar_read <- function(name,n)
{
  filename <- if(substring(name,nchar(name)-3)==".bar") name else paste(name,"bar",sep=".")
  ## filename <- if(tools::file_ext(name)=="par") name else paste(name,".par",sep="") # R 2.11 and later

  f <- file(filename,open="rb")
  r <- readBin(f,what="numeric",n=n)
  close(f)
  r
}

