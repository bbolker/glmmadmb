## functions to download glmmADMB binaries from buildbot

## buildbot_base <- "http://www.admb-project.org/buildbot/glmmadmb/glmmadmb-admb-"
buildbot_base <- "http://www.admb-project.org/buildbot/glmmadmb/glmmadmb-"
## assume we are starting from the root of the package directory
## 27 Dec: change gcc 4.5 -> 4.6
platform_str <- c(linux32="fedora18-gcc4.7-32bit",
                  linux64="ubuntu12-gcc4.6-64bit",
                  macos32="macos10.8-xcode4.6-32bit",
                  macos64="macos10.8-xcode4.6-64bit",
                  windows32="windows7-vc10-32bit",
                  windows64="windows7-vc10-64bit")
                  ## FIXME: 
                  ## windows64="windows7-vc10-64bit-64bit-"
                  ## windows64="windows8-mingw")
  ## c(linux32="linux-gcc4.5.2-32bit",
  ##   linux64="linux-gcc4.5.2-64bit",
  ##   macos32="macos10.6.7-xcode3.2.6-32bit",
  ##   macos64="macos10.6.7-xcode3.2.6-64bit")

## assumed that we are starting in the head package directory ...
if (!file.exists("inst/bin")) {
  setwd("../..")  ## try moving up
  if (!file.exists("inst/bin")) {
    stop("can't find bin directory -- are you in the package root?")
  }
}
setwd("inst/bin") ## move to bin directory

## get new versions of all binaries
get_allbin <- function(release) {
  if (missing(release)) {
    g <- suppressWarnings(get_bbot_versions())
    ## suppress "incomplete final line" warning
    m <- sapply(platform_str,grep,g$desc)
    OK <- sapply(m,length)>0
    m <- unlist(m)
    release <- paste("r",g[m,]$ver,sep="")
  }
  cdir <- getwd()
  on.exit(setwd(cdir))
  plist <- names(platform_str)[OK]
  status <- rep(NA,length(plist))
  names(status) <- plist
  for (i in seq_along(plist)) {
    p <- plist[i]
    cat(p,"\n")
    setwd(p)
    srcext <- if (grepl("windows",p)) ".exe" else ".bin"
    destext <- if (grepl("windows",p)) ".exe" else ""
    fn <- paste(buildbot_base,release[i],"-",
                platform_str[p],
                srcext,sep="")
    tt <- try(download.file(fn,
                            destfile=paste("glmmadmb",destext,sep="")))
    if (!inherits(tt,"try-error")) status[i] <- release[i]
    setwd("..")
  }
  status
}
  
## glmmadmb-r107-bcc5.5-32bit 	0B 	[text/html] 	
## glmmadmb-r107-linux-gcc4.4.3-32bit 	1M 	[text/html] 	
## glmmadmb-r107-linux-gcc4.5.2-32bit 	1M 	[text/html] 	
## glmmadmb-r107-linux-gcc4.5.2-64bit 	1M 	[text/html] 	
## glmmadmb-r107-linux-solarisstudio12-32bit 	2M 	[text/html] 	
## glmmadmb-r107-macos10.6.7-xcode3.2.6-32bit 	2M 	[text/html] 	
## glmmadmb-r107-macos10.6.7-xcode3.2.6-64bit 	2M 	[text/html] 

bburl <- "http://admb-project.org/buildbot/glmmadmb/"

get_bbot_versions <- function(os="all",rev="latest",bits="all",
                              bburl="http://admb-project.org/buildbot/glmmadmb/") {
    require("plyr")
    require("stringr")
    OSvals <- c("fedora","macos","ubuntu","windows")
    z <- suppressWarnings(readLines(url(bburl)))
    desc <- z[grep("glmmadmb",z)]
    sizes <- z[grep("[0-9][BKM]",z)]
    desc <- desc[-1]  ## header line
    desc <- gsub(".+\"(([a-z0-9.]|-)+).*","\\1",desc)
    ## assume 0B = 0 -- everything else will be in M
    sizes <- as.numeric(gsub(".+>([0-9.]+).*","\\1",sizes))
    d <- data.frame(desc,sizes)
    d <- subset(d,sizes>0)
    if (os!="all") d <- d[grep(os,d$desc),]
    if (bits!="all") d <- d[grep(paste(d$bits,"bit",sep="")),]
    if (rev=="latest") {
        ## d$type <- gsub("-r[0-9]+.*$","",d$desc)
        d$type <- gsub("\\.bin$","",gsub("(glmmadmb-|r[0-9]+-)","",d$desc))
        d$ver <- as.numeric(gsub(".*-r([0-9]+).*$","\\1",d$desc))
        d <- droplevels(ddply(d, .(type), function(x) x[which.max(x$ver),]))
    } else if (ver!="all") {
        d <- d[d$ver==ver,]
    }
    OSopts <- paste(OSvals,collapse="|")
    d <- mutate(d,
                desc=as.character(desc),
                type=as.character(type),
                sdesc=gsub("(glmmadmb|r[0-9]+|\\.bin|\\.exe)","",desc),
                glmmadmbver=as.numeric(gsub("[-r]","",
                str_extract(desc,"-r[0-9]+-"))),
                fullOS=str_extract(sdesc,paste0("(",OSopts,")[.0-9]+")),
                OSver=str_extract(fullOS,"[0-9.]+"),
                OSstr=str_extract(fullOS,"[[:alpha:]]+"))
    d <- mutate(d,sdesc=gsub("^-+","",gsub("--","-",sdesc)))
    return(d)
}

gg <- get_bbot_versions()
download_bin <- function(x,quiet=FALSE) {
  download.file(paste(bburl,x,sep=""),x)
}

test_bin <- function(x) {
  Sys.chmod(x)
  v <- try(system(paste("./",x," 1> tmp.log 2>&1",sep=""),intern=TRUE))
  readLines("tmp.log")
}

test_OK <- function(x) length(x)==3 && x[3]==" This is usual caused by a missing DAT file "

test_allbin <- function(...) {
  d <- get_bbot_versions(...)
  z <- lapply(d$desc,function(x) {
    download_bin(x)
    test_bin(x)
  })
  invisible(sapply(d$desc,unlink))
  unlink(c("b1","b2","s1","s2","glmmadmb*.log",
           "tmp.log","variance","eigv.rpt","\001",
           "fmin.log"))
  data.frame(type=d$type,ver=d$ver,OK=sapply(z,test_OK))
}

## test_allbin()
## get_allbin()

if (FALSE) {
    ## add to get_bbot_versions(); add latest.only
    gg <- get_bbot_versions()
    gg$sdesc <- gsub("(glmmadmb|-r[0-9]+|\\.bin|\\.exe)","",gg$desc)
    gg$sdesc <- gsub("^-","",gg$sdesc)
    table(gg$sdesc)
    library("stringr")
    glmmadmbver <- as.numeric(gsub("[-r]","",str_extract(gg$desc,"-r[0-9]+-")))
    fullOS <- str_extract(gg$sdesc,"(fedora|macos|ubuntu|macos|windows)[.0-9]+")
    OSver <- str_extract(fullOS,"[0-9.]+")
    OSstr <- str_extract(fullOS,"[[:alpha:]]+")
    gg2 <- with(gg,data.frame(OSstr,OSver,bits,ADMBver=ver))
}
## ISSUES
## Johnoel: adopt consistent naming format?  (maybe we can just be clever about parsing ...)
## back-compatibility; add MacOS < 10.9 builds?
## should we change the system to allow downloads of different binaries?  (maybe there's
##  a default binary but you can download a different one?)  If we have to include binaries
## for all OS versions/compilers the package is going to get really really big ...


## INSTRUCTIONS for manual replacement of binaries:
## 1. go to http://admb-project.org/buildbot/glmmadmb/
## 2. find the latest version (or the version you want) that's compatible with your OS
## 3. download it somewhere
## 4. from within R, glmmADMB:::get_bin_loc()
## 

if (FALSE){
    library("glmmADMB"
    bb <- glmmADMB:::get_bin_loc()[["bin_loc"]]
    bpath <- gsub("glmmadmb$","",bb)
    file.copy(bb,paste0(bpath,"glmmadmb.bak"))
    bburl <- "http://admb-project.org/buildbot/glmmadmb/"
    download.file(paste0(bburl,"glmmadmb-mingw64-r2885-windows8-mingw64.exe"),
                  dest=bb)
}

