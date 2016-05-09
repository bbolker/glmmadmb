pin_write <- function(name, L)
{
  filename <- if(substring(name,nchar(name)-3)==".pin") name else paste(name,".pin",sep="")
  ## filename <- if(tools::file_ext(name)=="pin") name else paste(name,".pin",sep="") # R 2.11 and later

  cat("# \"", name, ".pin\" produced by pin_write() from ADMButils; ", date(), "\n", file=filename, sep="")

  for(i in 1:length(L))
  {
    x <- L[[i]]
    if(data.class(x) == "numeric")
      cat("#", names(L)[i], "\n", L[[i]], "\n\n", file=filename, append=TRUE)
    if(data.class(x) == "matrix")
    {
      cat("#", names(L)[i], "\n", file=filename, append=TRUE)
      write.table(L[[i]], col.names=FALSE, row.names=FALSE, quote=FALSE, file=filename, append=TRUE)
      cat("\n", file=filename, append=TRUE)
    }
  }

  invisible(NULL)
}
